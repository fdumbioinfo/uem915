#' @title Quality Controls
#' @description
#' do descriptive statistics : histogram, boxplot , hierarchical clustering and ACP for column
#' @param dat matrix numeric
#' @param sif data.frame
#' @param inputdata logical to export input data or not
#' @param histo logical do histogram if TRUE by defaut
#' @param boxplot logical do boxplot if TRUE by defaut
#' @param hc logical do hierarchical clustering if TRUE by defaut
#' @param acp logical do principal component analysis if TRUE by defaut
#' @param dirname character
#' @param path character
#' @details return pval for each factor of anova model
#' function use doparalle
#' @return data.frame
#' @examples
#' # not run
#' # qc(dat, metadata)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @import gplots
#' @import ggplot2
#' @importFrom foreach foreach %dopar%
#' @importFrom graphics title hist
#' @export
qc <- function(
    dat , sif = NULL, inputdata = F,
    histo = TRUE, boxplot = TRUE, hc = TRUE, acp = TRUE,
    dirname = NULL , path = "." )
{
  dat[,-1] -> Dat1
  if( is.null( dirname ) ) { paste("QC_",ncol(Dat1),"_",nrow(Dat1),sep="") -> DirName }else{ dirname -> DirName  }
  path %>% file.path( DirName ) -> Path
  Path %>% dir.create
  #
  if( inputdata )
  {
    paste("data_",ncol(Dat1),"_",nrow(Dat1),".tsv",sep="") -> FileName
    Path %>% file.path( FileName ) -> FileName
    data.frame( "rowID" = dat[,1] %>% as.character, Dat1, stringsAsFactors=T ) -> Dat2
    Dat2 %>% output( FileName )
    # metadata
    if( !is.null( sif ) )
    {
      for( i in 1:dim(sif)[2] )
      {
        if( ( length( table( as.character(sif[,i]) ) ) < length( as.character( sif[,i] ) ) ) &
            length( table( as.character(sif[,i]) ) ) < 28 )
        { sif[,i]  <- factor( sif[,i] ) }
      }
      paste("metadata_",nrow(sif),".tsv",sep="") -> FileName
      Path %>% file.path( FileName ) -> FileName
      sif %>% output( FileName )
    }else
      {
        paste("metadata_",ncol(Dat1),".tsv",sep="") -> FileName
        Path %>% file.path( FileName ) -> FileName
        data.frame( SampleID = paste("s",1:ncol(Dat1),sep=""), FileName=colnames(Dat1) ) %>% output(FileName)
      }
  }
  # histogram
  if( histo )
  {
    paste("Histo_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep="") -> FileName
    Path %>% file.path(FileName) -> FileName
    paste("Histogram",ncol(Dat1)," x ",nrow(Dat1)) -> Title
    pdf( FileName )
    Dat1 %>% as.matrix %>% hist( breaks=70, main=Title )
    graphics.off()
  }
  # without factor
  if( is.null(sif) )
  {
    # boxplot
    if( boxplot )
    {
      paste( "boxplot_data_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep="" ) -> FileName
      Path %>% file.path( FileName ) -> FileName
      paste("Boxplot data ",ncol(Dat1)," x ",nrow(Dat1), sep = ""  ) -> Title
      # plot
      pdf( FileName )
      Dat1 %>% graphics::boxplot( outline = F , las = 2 )
      title( main = Title  )
      graphics.off()
    }
    # hc
    if( hc )
    {
      paste("hc_data_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
      Path %>% file.path( FileName ) -> FileName
      paste("Hierarchical clustering data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
      # plot
      pdf( FileName )
      Dat1 %>% select(-1) %>% hc( title = Title )
      graphics.off()
    }
    # acp
    if( acp )
    {
      paste("acp_data_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep = "") -> FileName
      Path %>% file.path( FileName ) -> FileName
      paste( "Principal Component Analysis data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
      # plot
      Dat1 %>% acp( title = Title )
      ggsave( FileName )
    }
  }
  # with factor
  if( !is.null(sif) )
  {
    sif %>% sapply( is.factor ) %>% which -> sel
    # boxplot
    if( boxplot )
    {
      Path %>% file.path( "boxplot" ) %>% dir.create
      foreach( i=1:length(sel) ) %do%
        {
          paste("boxplot_data_",colnames(sif)[sel[i] ],"_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "boxplot", FileName ) -> FileName
          paste("Boxplot data ",ncol(Dat1)," x ",nrow(Dat1), sep = ""  ) -> Title
          # plot
          pdf( FileName )
          Dat1 %>% boxplot( factor=sif[,sel[i]], outline=F, title=Title, legendtitle=colnames(sif)[sel[i]] )
          graphics.off()
        }
    }
    # hc
    if( hc )
    {
      Path %>% file.path( "hc") %>% dir.create
      foreach( i=1:length(sel) ) %do%
        {
          paste("hc_data_",colnames(sif)[ sel[i] ],"_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "hc", FileName ) -> FileName
          paste("Hierarchical clustering data ",ncol(Dat1)," x ",nrow(Dat1),sep = "") -> Title
          #
          pdf( FileName )
          Dat1 %>% hc( factor=sif[,sel[i]], title=Title, legendtitle=colnames( sif )[sel[i]] )
          graphics.off()
        }
    }
    # ACP
    if( acp )
    {
      Path %>% file.path( "acp") %>% dir.create
      # pc 1 2
      foreach( i=1:length(sel) ) %do%
        {
          paste("acp_data_",colnames(sif)[sel[i] ],"_pc12_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "acp", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
          #
          Dat1 %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=2, title=Title, legendtitle=colnames(sif)[sel[i]] )
          ggsave( FileName )
        }
      # pc 1 3
      foreach( i=1:length(sel) ) %do%
        {
          paste("acp_data_",colnames(sif)[ sel[i] ],"_pc13_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "acp", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
          #
          Dat1 %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=3, title=Title, legendtitle=colnames(sif)[sel[i]] )
          ggsave( FileName )
        }
      # pc 1 4
      foreach( i=1:length(sel) ) %do%
        {
          paste("acp_data_",colnames(sif)[ sel[i] ],"_pc14_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "acp", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(Dat1)," x ",nrow(Dat1), sep = ""  ) -> Title
          #
          Dat1 %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=4, title=Title, legendtitle=colnames(sif)[ sel[i] ] )
          ggsave( FileName )
        }
    }
  }
}