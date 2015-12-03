#' Polygon grid spanning the Indian subcontinent
#' 
#' A grid with 490 cells each \eqn{1^0} by \eqn{1^0} in size spanning the indian 
#' subcontinent. Each cell is assigned a hypothetical assemblage of varying
#' size from a regional pool of 100 species. The attribute table of the shape
#' file has one uniqe id field and three explanatory variables. 
#' \itemize{
#'    \item NUMMER a unique id
#'    \item ELEV_MEAN mean elevation of each grid
#'    \item ELEV_MIN minimum elevation of each grid
#'    \item ELEV_MAX maximum elevation of each grid}
#'    
#'@format an ESRI shapefile imported with \code{\link[maptools]{readShapeSpatial}} with 
#'        four columns in attribute table
#'@name shp
NULL