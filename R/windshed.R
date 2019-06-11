

# dispersal function to compute transition likelihood
diffuse <- function(x){
      p <- x[c(1,3,5,7,9,11,13,15)] # SW ...cw... S
      g <- x[17:20]
      if(g[1]<g[2] & g[4]<g[3]) return(p[1]/sqrt(2)/sqrt(2)) #SW
      if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
      if(g[2]<g[1] & g[4]<g[3]) return(p[3]/sqrt(2)/sqrt(2)) #NW
      if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
      if(g[2]<g[1] & g[3]<g[4]) return(p[5]/sqrt(2)/sqrt(2)) #NE
      if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
      if(g[1]<g[2] & g[3]<g[4]) return(p[7]/sqrt(2)/sqrt(2)) #SE
      if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
      # dividing by 2 root 2 prevents the geocorrection from distorting the probability field
}


transition_stack <- function(x, transitionFunction, directions, symm, ...){

      brk <- x
      x <- brk[[1]]

      tr <- new("TransitionLayer",
                nrows=as.integer(nrow(x)),
                ncols=as.integer(ncol(x)),
                extent=extent(x),
                crs=projection(x, asText=FALSE),
                transitionMatrix = Matrix(0,ncell(x),ncell(x)),
                transitionCells = 1:ncell(x))
      transitionMatr <- transitionMatrix(tr)
      Cells <- which(!is.na(getValues(x)))
      adj <- adjacent(x, cells=Cells, pairs=TRUE, target=Cells, directions=directions)

      ##### start modificiations #####
      # format raster data layers to feed to transitionFunction
      # col order is x[[1]][from, to] ... x[[n]][from,to]
      dataVals <- lapply(1:nlayers(brk),
                         function(i) cbind(values(brk[[i]])[adj[,1]],
                                           values(brk[[i]])[adj[,2]]))
      dataVals <- do.call("cbind", dataVals)
      ##### end modifications #####

      transition.values <- apply(dataVals, 1, transitionFunction, ...)
      transitionMatr[adj] <- as.vector(transition.values)
      transitionMatrix(tr) <- transitionMatr
      matrixValues(tr) <- "resistance"
      return(tr)
}



# create a transition object from a windrose raster stack
wind_trans <- function(windrose, correction="c"){
      windrose <- add_coords(windrose)
      trans <- transition_stack(windrose, diffuse, directions=8, symm=F)
      geoCorrection(trans, type=correction)
}


windshed <- function(trans, # transition object created by wind_trans()
                     coords, # lenghth-2 vector: lon-lat of focal location
                     upwind = F){
      coords <- matrix(coords, ncol=2)
      if(upwind){
            cost <- raster(trans)
            cost[] <- costDistance(trans, coordinates(cost), coords)
      } else{
            cost <- accCost(trans, coords)
      }
      return(cost)
}
