#' @title MSigDB enrichment analysis
#' @description MSigDB enrichment analysis
#' @param SymbolList character Symbol or NCBI gene ID
#' @param geneannot data.frame 
#' @param species character hs mm rn dr ss 
#' @param bg numeric 
#' @param filtergeneset regexp to filter geneset database
#' @param overlapmin numeric for minimum overlap between geneset and list
#' @param enaScoremin numeric for minimum ratio ena
#' @param top numeric top features to plot
#' @param labsize numeric size of function in barplot
#' @param dpibarplot character barplot resolution
#' @param path character for relative path of output directory
#' @param dirname character name for output
#' @return file with enrichment analysis results
#' @examples
#' # not run
#' # ena( Symbollist , filtergeneset = "reactome")
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select
#' @importFrom rlang .data
#' @importFrom stats binom.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @export
ena <- function(
    SymbolList = NULL, geneannot = NULL, species = "hs",
    bg = 25000 , filtergeneset = "all",
    overlapmin = 2 , enaScoremin = 1,
    top = 80, labsize = 11, dpibarplot = "screen", path = "." , dirname = NULL )
{
  i=j=1
  ifelse( species != "hs" , Ortholog <- T , Ortholog <- FALSE )
  if( is.null(geneannot) )
  {
    SymbolList %>% "!="( . , "" ) %>% which %>% SymbolList[.] %>% unique -> SymbolList
    SymbolList %>% annot( species = species, ortholog = Ortholog ) -> a0
  }else{ geneannot -> a0 }
  #
  if( nrow(a0) > 0 )
  {
    a0$GeneID -> GeneidList
    # enainfo : match symbol
    enainfo <- list()
    a0$Symbol %>% length -> enainfo[[1]]
    # enainfo : total genesets db
    moalannotgene::genesetdb -> GenesetDb0
    GenesetDb0 %>% unlist(recursive = F) %>% length -> enainfo[[2]]
    ### genesetdb filtering
    GenesetDb0 %>% names %>% grep( filtergeneset , . , ignore.case = T  ) %>% GenesetDb0[.] -> GenesetDb1
    if( filtergeneset == "all" ){ GenesetDb0 -> GenesetDb1 }
    # enainfo : filtergeneset
    filtergeneset -> enainfo[[3]]
    # enainfo : genesets number after filtering
    GenesetDb1 %>% sapply( length ) %>% sum -> enainfo[[4]]
    # enainfo : background value
    bg -> enainfo[[5]]
    # ena computing
    ena0 <- foreach( i=1:length( GenesetDb1 )) %:%
        foreach( j=1:length( GenesetDb1[[i]] ), .combine = "rbind" ) %do%
      {
        GenesetDb1[[i]][[j]][[1]] -> GenesetName
        GeneidList %in% GenesetDb1[[i]][[j]][[4]] %>% which %>% a0$GeneID[.] -> OverlapGeneIDList
        GeneidList %in% GenesetDb1[[i]][[j]][[4]] %>% which %>% a0$Symbol[.] -> OverlapSymbolList
        OverlapSymbolList %>% length -> OverlapSize
        GenesetDb1[[i]][[j]][[4]] %>% length -> GenesetSize
        ( OverlapSize / GenesetSize ) %>% round( 2 ) -> OverlapRatio
        ( ( OverlapSize / a0$Symbol %>% length ) / (  GenesetSize / bg ) ) %>% round( 2 ) -> EnaScore
        if(  OverlapSize > overlapmin & EnaScore > enaScoremin )
        {
          binom.test( OverlapSize, a0$Symbol %>% length, ( GenesetSize / bg ) )[3] %>% unlist -> pval
          c(
            GenesetName,
            paste0( OverlapSymbolList, collapse = "|" ),
            paste0( OverlapGeneIDList, collapse = "|" ),
            OverlapSize,GenesetSize,OverlapRatio,EnaScore,pval)
        }
      }
    ena0 %>% setNames( GenesetDb1 %>% names ) -> ena0
    ena0 %>% sapply(is.null) %>% "!"(.) %>% which %>% ena0[.] -> ena0
    #
    # ena output
    #
    HTML <- F
    if( length(ena0) > 0 )
    {
      HTML <- T
      # output file
      ifelse( is.null(dirname) , "ena_ListName" -> DirName , paste("ena_",dirname,sep="") -> DirName )
      path %>% file.path( DirName ) -> Path
      Path %>% dir.create
      Path %>% file.path( "files" ) %>% dir.create
      #
      foreach( i=1:length(ena0) ) %do%
        {
          #
          # table
          #
          ena0[[i]] %>% data.frame( stringsAsFactors = F ) %>%
            stats::setNames(
              c("Name","SymbolList","GeneIDList","OverlapSize","GenesetSize",
                "OverlapRatio","ENAScore","pval" ) ) -> ena1
          ena1$OverlapSize %>% as.numeric -> ena1$OverlapSize
          ena1$GenesetSize %>% as.numeric -> ena1$GenesetSize
          ena1$OverlapRatio %>% as.numeric -> ena1$OverlapRatio
          ena1$ENAScore %>% as.numeric -> ena1$ENAScore
          ena1$pval %>% as.numeric -> ena1$pval
          ena1 %>% mutate( pvalFDR = .data$pval %>% p.adjust(method = "fdr") ) %>%
            mutate( log10pvalFDR = .data$pvalFDR %>% log10 %>% '*'( . , -1 ) %>% round( . , 4) ) %>%
            arrange( .data$pvalFDR  ) -> ena1
          # output
          paste(
            "ena_" ,ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1), "_",
            ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2),"_",
            dim(ena1)[1] , ".tsv" , sep = "" ) -> FileName
          Path %>% file.path( "files", FileName  ) %>% output( ena1 , . )
          #
          # barplot
          #
          ena1 %>% dplyr::select( c(1,10,7) ) -> ena2
          if( nrow(ena2) > 0 )
          {
            paste( ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1),
                   ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2) , sep = " "  ) -> title
            
            enabarplot( dat = ena2 , top = top , labsize = labsize , title = title ) -> p
            # output
            ifelse( nrow(ena2) < top , topfilename <- nrow(ena2) , topfilename <- top )
            paste(
              "ena_barplot_",top, "_" ,ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1), "_",
              ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2),"_",
              topfilename, ".jpeg" , sep = "" ) -> FileName
            Path %>% file.path( "files" , FileName ) -> FileName1
            ggsave( filename=FileName1 , plot = p , width = 12, height = 15 , dpi = dpibarplot )
          }
        }
      }
    #
    # html
    #
    if( HTML )
    {
      #
      # ena to html
      #
      Path %>% list.files( full.names = T, recursive = T) -> l0
      l0 %>% basename %>% grep(".jpeg$" , . ) %>% l0[.]  -> l1 # summary barplot 
      l0 %>% basename %>% grep("^ena_.*.tsv" , . ) %>% l0[.] -> l0
      foreach( i = 1:length( l0 ),.packages = c("magrittr") ) %do% { enatohtmldo( file.path = l0[i] ) }
      # 
      # summary barplot 
      #
      Path %>% file.path( "files" , "barplot_summary.html" ) %>% file( "w") -> f
      paste('<html>') %>% writeLines(f)
      paste('<head>') %>% writeLines(f)
      paste('<meta charset="utf-8"/>') %>% writeLines(f)
      paste('<title> Barplot collection summary </title>') %>% writeLines(f)
      paste('</head>') %>% writeLines(f)
      paste('<body>') %>% writeLines(f)
      foreach( i=1:length(l1) ) %do%
        {  
          paste("<pr>") %>%
            paste("<img src=",l1[i] %>% basename,"style='width:150px;height:150px;'>") %>%
            paste("</pr>") %>% writeLines(f)
        }
      paste('</body>') %>% writeLines(f)
      paste('</html>') %>% writeLines(f)
      close(f)
      #
      # summary html
      #
      # file
      Path %>% file.path( "enrichment_analysis.html" ) %>% file( "w") -> f
      paste('<html>') %>% writeLines(f)
      # head
      paste('<head>') %>% writeLines(f)
      paste('<meta charset="utf-8"/>') %>% writeLines(f)
      paste('<title> Results analysis </title>') %>% writeLines(f)
      paste('</head>') %>% writeLines(f)
      paste('<body>') %>% writeLines(f)
      paste('<p><b><FONT color = "red" size = "5pt"> Enrichment analysis results </FONT></b></p>') %>% writeLines(f)
      # ena info
      paste('<p><b> Match Symbols : ', enainfo[[1]] ,'</b></p>') %>% writeLines(f)
      paste('<p><b> Total genesets : ', enainfo[[2]] ,'</b></p>') %>% writeLines(f)
      paste('<p><b> geneset filter : ', enainfo[[3]] ,'</b></p>') %>% writeLines(f)
      paste('<p><b> genesets used : ', enainfo[[4]] ,'</b></p>') %>% writeLines(f)
      paste('<p><b> Background value : ', enainfo[[5]] ,'</b></p>') %>% writeLines(f)
      paste('<p><b><a href = "https://www.gsea-msigdb.org/gsea/msigdb/index.jsp"> Molecular Signatures Database v7.5.1: ',
            '</a></b></p>' , sep = "") %>% writeLines(f)
      # ena list info
      paste('<table border = 1 >') %>%
        paste('<thead>') %>%
        paste('<th>Collection</th>') %>%
        paste('<th>Database</th>') %>%
        paste('<th>All</th>') %>%
        paste('<th>Results</th>') %>%
        paste('</thead>') %>% writeLines(f)
      foreach( i=1:length(ena0) ) %do%
        {
          l0[i] %>% basename %>% sub("^ena_(.*_.*)_.*_.*.tsv","\\1",.) -> Collection0
          l0[i] %>% basename %>% sub("^ena_.*_.*_(.*)_.*.tsv","\\1",.) -> Database0
          GenesetDb0 %>% names %>% grep(Database0,.,value=T) %>% strsplit("\\|") %>% unlist %>% "["(3) -> GeneSetNumber
          paste( '<tr>') %>%
            paste( '<td>' ,  Collection0 , '</td>' ) %>%
            paste( '<td>' , Database0 , '</td>' ) %>%
            paste( '<td>' ,  GeneSetNumber , '</td>' ) %>% writeLines(f)
          l0[i] %>% sub( "\\.tsv" ,".html" , . ) %>% basename %>% file.path("files" , . ) -> link
          paste('<td><a href = "', link  ,'">' , sep = "") %>%
              paste( l0[i] %>% basename %>% sub( ".*_(.*).tsv$" ,"\\1" , . ) ) %>%
            paste('</a></td>') %>%
          paste( '</tr>') %>% writeLines(f)
        }
      # Path %>% file.path( DirName , "files", "barplot_summary.html" ) -> link 
      file.path( "files", "barplot_summary.html" ) -> link 
      paste('<p><b><a href = ',link,'> Barplot summary</a></b></p>' , sep = "") %>% writeLines(f)
      paste('</body>') %>% writeLines(f)
      paste('</html>') %>% writeLines(f)
      close(f)
    }
  }
}