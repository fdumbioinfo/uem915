#' @title loading regular libraries
#' @description load magrittr, dplyr, gplots, ggplot2, foreach, parallel, doParallel
#' @export
env <- function()
{
  c( "magrittr","dplyr","gplots","ggplot2","foreach","parallel","doParallel") %>%
    sapply( library , character.only = T )
}