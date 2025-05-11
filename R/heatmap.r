#' @title Heatmap
#' @description To make a heatmap
#' @param dat matrix numeric
#' @param factor factor
#' @param method character
#' @param dendrogram character to display 'none', 'row', 'column' or 'both' (by default) dendrograms
#' @param labCol character
#' @param cexCol numeric
#' @param labRow Character
#' @param cexRow numeric
#' @param cexlegend numeric
#' @param k numeric number of clusters to colorize for rows
#' @param keysize numeric
#' @param keycolor character of 3 for low mid high value of the key
#' @param parmar numeric 4 values for margin sizes
#' @details To make a heatmap from a matrix or a data.frame
#' @return no returned value
#' @examples
#' # not run
#' # library(magrittr)
#' # data(sif1)
#' # data(mat1)
#' # mat1 %>% heatmap(sif1$F3)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom graphics par legend
#' @importFrom gplots heatmap.2 colorpanel
#' @importFrom dendextend set
#' @export
heatmap <- function(
  dat , factor , method = "complete" , dendrogram = "both", k = NULL,
  labCol = "" , cexCol = 0.85 , labRow = "" , cexRow = 0.35,
  cexlegend = 0.65 , keysize = 0.9, keycolor = c( "darkgreen", "orange" , "darkred") , parmar = c(5,4,5,6))
{
  par(xpd=T, mar=parmar)
  if(dendrogram == "both")
  {
    dat %>% hc(factor=factor, plot=F , method = method ) -> coldend
    dat %>% t %>% hc(plot=F, method=method) -> rowdend
  }
  if( dendrogram == "none" )
  {
    coldend <- F
    rowdend <- F
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if( dendrogram == "row" )
  {
    coldend <- F
    rowdend <- T
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if(dendrogram == "col")
  {
    coldend <- T
    rowdend <- F
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if( !is.null(k) ){ rowdend %>% dendextend::set( "branches_k_color", value=palette0[1:k], k=k ) -> rowdend }
  dat %>% as.matrix %>%
    gplots::heatmap.2(
      ColSideColors=factortocolor(factor), labCol=labCol, cexCol=cexCol, srtCol=90, labRow=labRow, cexRow=cexRow, srtRow = 0,
      trace="none", key=T, keysize=keysize, density.info="none", key.xlab="", key.title=NA,
      dendrogram=dendrogram, Colv=coldend, Rowv=rowdend,
      col=colorpanel(50, low=keycolor[1], mid=keycolor[2], high=keycolor[3]), scale="row", distfun=dist2 )
  legend(
    "topright", inset = c(-0.22,-0.1), legend = levels( factor( factor ) ),
    col = palette0[1:length( levels( factor( factor ) ) ) ],
    lty = 1, lwd = 5, cex = cexlegend )
}