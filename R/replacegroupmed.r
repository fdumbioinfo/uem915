#' @title replace value by group median
#' @param dat character vector of mzxml file
#' @param value numeric value to substitute
#' @param factor character
#' @details
#' replace value by column median if group size is one and replace by row group median if group size > 1
#' @return data.frame
#' @examples
#' # not run
#' # replacegroupmed(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach  %do% foreach
#' @importFrom stats median setNames
#' @export
replacegroupmed <- function( dat , value = 0 , factor )
{
  i=j=1
  if( dat %>% unlist %>% "=="(0) %>% which %>% length %>% ">"(0) )
  {
    Dat1 <- foreach( i=1:length(factor %>% levels) , .combine = "cbind") %do%
      {
        factor %>% levels %>% "["(.,i) %>% "=="(.,factor) %>% which -> sel
        dat[,sel] %>% as.data.frame -> Dat0
        Dat0 %>% unlist %>% median -> MedGroupCol
        Dat0 %>% "=="( ., value ) %>% which( arr.ind=T ) -> selvalue
        if( nrow(selvalue) > 0 )
        {
          foreach( j=1:nrow(selvalue) ) %do%
            {
              Dat0[selvalue[j,1],] %>% as.numeric %>% median -> MedGroupRow
              if( MedGroupRow %>% "=="(.,0) ){ Dat0[ selvalue[j,1] , selvalue[j,2] ] <- MedGroupCol }else
              { Dat0[ selvalue[j,1] , selvalue[j,2] ] <- MedGroupRow }
            }
        }
        Dat0 %>% as.data.frame
      }
  }
  Dat1 %>% setNames( colnames(dat) ) %>% return()
}