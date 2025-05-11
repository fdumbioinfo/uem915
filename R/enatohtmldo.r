#' @title Convert tsv ena function result file into html file
#' @description Convert tsv ena function result file into html file.
#' @param file.path character path to file to convert
#' @return file in html format
#' @details
#' Convert tsv ena function result file into html file
#' @examples
#' # not run
#' # filename %>% enatohtml
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach foreach %do%
#' @noRd
enatohtmldo <- function( file.path )
{
  file.path %>% input -> ena
  file.path %>% basename %>% gsub("tsv" , "html" , . ) -> FileName
  file.path %>% dirname %>% file.path( FileName ) %>% file( "w") -> f
  # html
  paste('<html>') %>% writeLines(f)
  # head
  paste('<head>') %>% writeLines(f)
  paste('<meta charset="utf-8"/>') %>% writeLines(f)
  paste('<title> Results analysis </title>') %>% writeLines(f)
  paste('</head>') %>% writeLines(f)
  # body
  paste('<body>') %>% writeLines(f)
  # table
  paste('<table border = 1 >') %>% writeLines(f)
  paste('<thead>') %>% writeLines(f)
  paste('<th>  </th>') %>% writeLines(f)
  paste('<th> Name </th>') %>% writeLines(f)
  paste('<th> Overlap Symbol list </th>') %>% writeLines(f)
  paste('<th> Overlap Size </th>') %>% writeLines(f)
  paste('<th> Geneset Size </th>') %>% writeLines(f)
  paste('<th> ENAScore </th>') %>% writeLines(f)
  paste('<th> log10(pFDR) </th>') %>% writeLines(f)
  paste('</thead>') %>% writeLines(f)
  # table body
  paste('<tbody>') %>% writeLines(f)
  i=1
  foreach( i=1:dim(ena)[1] ) %do%
  {
    paste('<tr>') %>% writeLines(f)
    # row number
    paste('<td>') %>% paste(i) %>%
      paste('</td>') %>% writeLines(f)
    # Name
    paste('<td><b>') %>% writeLines(f)
    paste('<a href = "https://www.gsea-msigdb.org/gsea/msigdb/cards/',ena[ i , 1 ],'.html">' , sep = "") %>%
      paste(ena[ i , 1 ]) %>%
      paste('</a>') %>%
      paste('</b></td>') %>% writeLines(f)
    # Overlap Symbols list
    ena$GeneIDList[i] %>% strsplit("\\|") %>% unlist -> geneidList
    ena$SymbolList[i] %>% strsplit("\\|") %>% unlist -> SymboldList
    geneidList -> geneidList1
    SymboldList -> SymboldList1
    if( length(geneidList) > 30 )
    {
      geneidList[1:30] -> geneidList1
      SymboldList[1:30] -> SymboldList1
    }
    paste('<td>') %>% writeLines(f)
    paste('<a href = "https://www.ncbi.nlm.nih.gov/gene/?term=',geneidList1,'">', sep = "") %>%
      paste('<FONT size = "2pt">', SymboldList1 , '</FONT>') %>%
      paste('</a>') %>%
      paste("|") %>% writeLines(f)
    if( length(geneidList) > 30 ) { paste("...") %>% writeLines(f) }
    paste('</td>') %>% writeLines(f)
    # Overlap Size
    paste('<td>') %>%
      paste(ena$OverlapSize[i]) %>%
      paste('</td>') %>% writeLines(f)
    # geneset Size
    paste('<td>') %>%
      paste(ena$GenesetSize[i]) %>%
      paste('</td>') %>% writeLines(f)
    # ENAScore
    paste('<td>') %>%
      paste( ena$ENAScore[i] %>% round(.,1) ) %>%
      paste('</td>') %>% writeLines(f)
    # Overlap pval
    paste('<td>') %>%
      paste( ena$log10pvalFDR[i] %>% round(.,1) ) %>%
      paste('</td>') %>% writeLines(f)
    paste('</tr>') %>% writeLines(f)#end row
  }
  paste('</tbody>') %>% writeLines(f)
  paste('</body>') %>% writeLines(f)
  paste('</html>') %>% writeLines(f)
  close(f)
}