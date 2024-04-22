setClass("wind_graph",
        contains = "TransitionLayer",
        slots = c(p = "numeric",
                  direction = "character"))

#' Build a wind connectivity graph
#'
#' This function constructs a landscape connectivity graph from a windrose raster stack.
#'
#' @param x An object of class `wind_rose`.
#' @param direction Either "downwind" (the default) or "upwind", indicating whether outbound or inbound wind conductance should be computed, respectively.
#' @return a `wind_graph` object consisting of a gdistance \link[gdistance]{Transition-class} object and additional metadata
#' @export
#' @aliases build_wind_graph
wind_graph <- function(x, direction = "downwind"){
      g <- transition_stack(add_coords(x), windflow, directions = 8, symm = F, direction = direction)
      g <- as(g, "wind_graph")
      g@p <- g@p
      g@direction <- direction
      g
}


#' A modified version of `gdistance::transition`,
#'
#' This version accepts a raster stack and a custom transition function; the
#' original does not support arbitrary transition functions for multi-layer
#' transition data.
transition_stack <- function(x, transitionFunction, directions, symm, ...){

      # require(raster)
      brk <- raster::stack(x)
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

      ##### start modifications #####
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

