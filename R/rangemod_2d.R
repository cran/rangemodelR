#' Range cohesion models for spatial polygon grids
#' @description species richness model for two dimentional grid with geometric
#'              constraints and range cohesion
#' @param spmat a site by species matrix or data frame with species in columns
#' @param shp shapefile of sites where species occurences are recorded
#' @param reps number of replicates
#' @param nb a neighbour object similar to one generated from 
#'        \code{\link[spdep]{poly2nb}} of \pkg{\link[spdep]{spdep}}
#'        If 'nb' is NA then result is range scatter
#' @param field a number or character vector indicating which column in the dbf
#'        of shapefile is the unique id
#' @param var an optional vector containing explanatory variable for
#'        constraining the randomization
#' @param first If true, 'var' is used while choosing the first occurence as
#'        well.if 'var' is null, first is always set 'FALSE'
#' @param degen If true, each randomized site by species matrix is saved and
#'        provided in output
#' @param rsize which rangesizes to use for simulation, can be an integer vector
#'        of same length as number of species(collumns) 
#'        or either 'observed' or'unif'. See details for explanations
#' @details rangemod.2d impliments simulations used by Rahbeck et.al. (2007) to
#'          species distribution data on a continuous grid. In 'spmat' the sites
#'          (rows) represent each cell in the grid.The species occurences across
#'          sites are randomly spread maintaining strict range cohesion.
#'          A neighbour object is used to limit the choice
#'          of cells during random selctions to immidiate neighbours. Options 
#'          for creating four cell (rook) or eight cell (queen) neighbours can
#'          be accessed while creating the 'nb' object, (typically from package 
#'          \pkg{\link[spdep]{spdep}}). The randomisation proceeds by selecting a
#'          single site,(weighted by 'var' if provided, and first is TRUE), and
#'          then continues selecting one site at a time from a vector of
#'          available neighbours taken from 'nb' and weighted by 'var' if
#'          provided. The vector of available sites is updated after each site
#'          is selected. 
#' @return  A list containing following elements:
#' \itemize{
#'  \item{"out.df"}{ a data frame with two columns 'mod.rich' and 'sd.rich' 
#'                    with mean and standard deviation respectively of the
#'                    predicted species richness }
#'  \item{"out.shp"}{ same as the input shapefile with the two colums of
#'                    'out.df' added in attribute table}
#'  \item{"degenerate.matrices"} {a list of all the randomized matrices(only 
#'                                present if 'degen' is TRUE)}                  
#'  }
#' @references Rahbek, C., Gotelli, N., Colwell, R., Entsminger, G., Rangel,
#'              T. & Graves, G. (2007) Predicting continental-scale patterns of
#'              bird species richness with spatially explicit
#'              models. Proceedings of the Royal Society B: Biological
#'              Sciences, 274, 165.
#'              
#'              Gotelli, N.J., Anderson, M.J., Arita, H.T., Chao, A., Colwell,
#'              R.K., Connolly, S.R., Currie, D.J., Dunn, R.R., Graves, G.R. &
#'              Green, J.L. (2009) Patterns and causes of species richness:
#'              a general simulation model for macroecology. Ecology Letters,
#'              12, 873-886.  
#' @examples \dontrun{
#' data(shp)
#' data(neigh_ob)
#' data(spmat)
#' library(ggplot2)
#' library(plyr)
#' mod.out <- rangemod.2d(spmat,shp,"NUMMER",nb = neigh_ob,rsize = "observed",
#'                        var = NULL,reps = 5,degen = F,first = F)
#' shp.out <- mod.out$out.shp
#' shp.out$id <- shp.out@@data[,"NUMMER"]
#' shp.out.df <- as.data.frame(shp.out)
#' shp.out.fort <- fortify(shp.out,region = "id")
#' shp.out.gg <- join(shp.out.fort,shp.out.df,by = "id")
#' ggplot(shp.out.gg)+
#'   geom_map(map=shp.out.gg,aes_string("long","lat",map_id="NUMMER",
#'                                      fill = "mod.rich"))+
#'   geom_path(aes(x = long,y = lat,group = group),colour = "white")+
#'   coord_equal()
#' }
#' @export
rangemod.2d <- function(spmat,shp,field,nb,rsize = c("observed","unif"),var,
                        reps,degen,first){
  if (!requireNamespace("maptools", quietly = TRUE)) {
    stop("'maptools' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  ####sanity check of arguments####
  if(any(!(shp@data[,field]%in%rownames(spmat)))){
    stop("all unique ids in 'shp' should appear in 'spmat'")
  }
  
  if(!is.na(nb)&&!length(nb) == nrow(spmat)){
    stop("length of 'nb' should be same as number of sites: ",
         length(nb)," and ", nrow(spmat))
  }
  
  if(!is.null(var)&& !length(var) == nrow(spmat)){
    stop("'var' should be of same length as number of sites: ",
         length(var),"and",nrow(spmat),".")
  }
  
  if(is.vector(rsize,mode = "numeric")&& !length(rsize) == ncol(spmat)){
    stop("rsize should be of same length as number of species: ",
         length(rsize)," and ",ncol(spmat))
  }

  ####chunk1 - prepare objects for input and output####
  spmat[spmat>0] <- 1
  if(is.vector(rsize,mode = "numeric")){
    range.size <- rsize
  }else{
    rsize <- match.arg(rsize)
    range.size <- switch(rsize,observed = {colSums(spmat)},
                         unif = {sample(1:nrow(spmat),ncol(spmat),replace = T)})
  }
  mat.temp <- as.matrix(spmat)
  mat.out <- matrix(nrow = nrow(spmat),ncol = reps,
                    dimnames = list(rownames(spmat),1:reps))
  uid <- shp@data[,field]
  degen.mats <- list()
  
  ####chunk2 - use 'random.range' to spread ranges on the matrix####
  if(degen == TRUE){
    for(i in 1:reps){
      mat.temp[which(mat.temp > 0)] <- 0
      for(j in 1:length(range.size)){
        temp.vec1 <- random.range(uid = uid,nb=nb,
                                  range.size = range.size[j],var = var,
                                  first = first)
        mat.temp[which(rownames(mat.temp)%in%as.character(temp.vec1)),j] <- 1
      }
      mat.out[,i] <- rowSums(mat.temp)
      degen.mats[[i]] <- mat.temp
    }
  }else{
    for(i in 1:reps){
      mat.temp[which(mat.temp > 0)] <- 0
      for(j in 1:length(range.size)){
        temp.vec1 <- random.range(uid = uid,nb=nb,
                                  range.size = range.size[j],var = var,
                                  first = first)
        mat.temp[which(rownames(mat.temp)%in%as.character(temp.vec1)),j] <- 1
      }
      mat.out[,i] <- rowSums(mat.temp)
      
    }
  }
  
  out.df <- data.frame(mod.rich = apply(mat.out,1,mean),
                       mod.sd =  apply(mat.out,1,stats::sd),
                       q2.5 = apply(mat.out,1,stats::quantile,probs = 0.025),
                       q95.5 = apply(mat.out,1,stats::quantile,probs = 0.975))
  
  ####chunk3 - add simulated data to shp and prepare for plotting in ggplot2
  shp@data <- data.frame(shp@data,out.df)
  
  if(degen==TRUE){
    outlist <- list(out.df = out.df,out.shp = shp,
                    degenerate.matrices=degen.mats)
  }else{
    outlist <- list(out.df = out.df,out.shp = shp)
  }
  outlist
}