#' @export
#' @author Nick McKay
#' @title Read varve thicknesses from shape file
#' @description Pulls varve thicknesses and error codes from a shapefile counted in GIS
#' @importFrom sf st_read
#' @import dplyr
#' @import magrittr
#' @param filename the path to a shape file
#' @return A data.frame with the varve data
readVarveShapefile <- function(filename,codeCol = "conf",varveDir = "vertical"){
shape <-  sf::st_read(filename,quiet = TRUE)

varveCoords <- data.frame(x1 = unlist(lapply(shape$geometry,"[[",1)),
                          x2 = unlist(lapply(shape$geometry,"[[",2)),
                          y1 = unlist(lapply(shape$geometry,"[[",3)),
                          y2 = unlist(lapply(shape$geometry,"[[",4)))

if(varveDir == "vertical"){
varves <- varveCoords %>%
  dplyr::arrange(desc(y1)) %>%
  dplyr::mutate(count = seq_along(y1)) %>%
  dplyr::mutate(thick = abs(y1-y2))
}else if(varveDir == "horizontal"){
  varves <- varveCoords %>%
    dplyr::arrange(desc(x1)) %>%
    dplyr::mutate(count = seq_along(x1)) %>%
    dplyr::mutate(thick = abs(x1-x2))
}

out <- dplyr::select(varves,count,thick)

if(any(names(shape) == codeCol)){
  out$varveCode = shape[[codeCol]]
}



return(out)
}
