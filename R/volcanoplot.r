#' @title Volcanoplot
#' @description To make a volcanoplot
#' @param dat data.frame
#' @param pval numeric p-value threshold
#' @param fc numeric fold-change threshold
#' @param GeneName logical displya GeneName or not. TRUE by default
#' @param GeneNameN logical the number of gene to be displayed. 50 by default
#' @param GeneNameList character vector of gene Symbol. see details.
#' @param GeneNameSize numeric
#' @param title character
#' @details
#' dat argument must have 4 column: featureid , pvalue , foldchange , Symbol.
#' Only Symbol column must have this name.
#' If GeneNameList not supplied, first 50 up and down features will be display.
#' @return no returned value
#' @examples
#' # not run
#' # volcanoplot()
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr slice
#' @import ggplot2
#' @noRd
volcanoplot <- function(
  dat, pval = 0.05, fc = 1.5,
  GeneName = FALSE , GeneNameN = 50 , GeneNameList = NULL, GeneNameSize = 2 ,
  title = "Volcanoplot" )
{
  dat[,2] %>% log10 %>% "*"( . , -1) -> dat[,2]
  dat[,3] %>% fctolog2ratio -> dat[,3]
  # up
  dat %>% dplyr::slice( which( dat[,2] %>% ">"(.,-log10(pval)) & dat[,3] > log2(fc))) -> Up
  # down
  dat %>% dplyr::slice( which( dat[,2] %>% ">"(.,-log10(pval)) & dat[,3] < -log2(fc))) -> Down
  ggplot() -> p
  p + geom_point(aes(x=dat[,3], y=dat[,2]), colour="black") -> p
  # up
  p + geom_point(aes(x=Up[,3], y=Up[,2]), colour="red") -> p
  # down
  p + geom_point(aes(x=Down[,3], y=Down[,2]), colour="green4") -> p
  if(GeneName)
  {
    if(is.null(GeneNameList))
    {
      # display GeneName up
      dat %>% dplyr::slice(which(dat[,3] > 1)) -> Dat1
      Dat1[order(-Dat1[,2]),] %>% dplyr::slice(1:GeneNameN)-> GeneNameUp
      p + geom_text(aes(x=GeneNameUp[,3], y=GeneNameUp[,2]), label=GeneNameUp[,4], size=GeneNameSize) -> p
      # display GeneName down
      dat %>% dplyr::slice(which(dat[,3] < -1)) -> Dat1
      Dat1[order(-Dat1[,2]),] %>% dplyr::slice(1:GeneNameN) -> GeneNameDown
      p + geom_text(aes(x=GeneNameDown[,3], y=GeneNameDown[,2]), label=GeneNameDown[,4], size=GeneNameSize) -> p
    }else
      {
        GeneNameList %>% as.character %>% unique %>% data.frame(Symbol=.) %>% inner_join(dat, by="Symbol") %>% select(c(2,3,4,1)) -> GeneName
        # up
        GeneName[order(-GeneName[,2]),] %>% dplyr::slice(1:GeneNameN) -> GeneNameUp
        GeneNameUp %>% slice(which(GeneNameUp[,3] > 1)) -> GeneNameUp
        p + aes(x=GeneNameUp[,3], y=GeneNameUp[,2]) -> p
        p + geom_point(aes(x=GeneNameUp[,3], y=GeneNameUp[,2]), colour="orange") -> p
        p + geom_text(aes(x=GeneNameUp[,3], y=GeneNameUp[,2]), label=GeneNameUp[,4], size=GeneNameSize, fontface=2, hjust=-0.2) -> p
        # down
        GeneName[order(-GeneName[,2]),] %>% dplyr::slice(1:GeneNameN )-> GeneNameDown
        GeneNameDown %>% dplyr::slice(which(GeneNameDown[,3] < -1)) -> GeneNameDown
        p + geom_point(aes(x=GeneNameDown[,3], y=GeneNameDown[,2]), colour="orange") -> p
        p + geom_text(aes(x=GeneNameDown[,3], y=GeneNameDown[,2]), label=GeneNameDown[,4], size=GeneNameSize, fontface=2, hjust=1) -> p
      }
  }
  p + theme_bw() -> p
  p + theme(plot.title=element_text(size=10, face="bold", hjust=0.5), legend.position="none") -> p
  p + geom_hline(yintercept=-log10(pval), linetype="solid", colour="orange", size=0.1) -> p
  p + geom_vline(xintercept=log2(fc), linetype="solid", colour="orange", size=0.1) -> p
  p + geom_vline(xintercept=-log2(fc), linetype="solid", colour="orange", size=0.1 ) -> p
  p + labs(x=paste("|log2ratio| > ",fc,sep=""), y=paste("-log10(",pval,")",sep="")) -> p
  p + ggtitle(paste(title, sep="")) -> p
  p %>% return()
}