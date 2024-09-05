#' An S4 object class representing a wind graph
#'
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
#' @param wrap Should the left and right edges of the raster be connected? Default is FALSE. Set to TRUE if, for example, `x` is a global raster where -180 and 180 are equivalent longitudes.
#' @return a `wind_graph` object consisting of a gdistance \link[gdistance]{Transition-class} object and additional metadata
#' @export
#' @aliases build_wind_graph
wind_graph <- function(x, direction = "downwind", wrap = FALSE){
      g <- transition_stack(add_coords(x), flow, directions = 8, symm = F, wrap = wrap, direction = direction)
      g <- as(g, "wind_graph")
      g@p <- g@p
      g@direction <- direction
      g
}


#' A wind transition function to supply to `transition_stack`.
#'
#' @param x as used internally within `transition_stack`, this is a vector of
#'   collated windrose values for FROM and TO cells
#' @param direction either "downwind" (the default) or "upwind", indicating
#'   whether outbound or inbound wind conductance should be computed.
flow <- function(x, direction = "downwind"){

      # node row and column indices
      g <- x[17:20]

      # edge wind loadings (clockwise from southwest)
      p <- c(x[1]+x[2],
             x[3]+x[4],
             x[5]+x[6],
             x[7]+x[8],
             x[9]+x[10],
             x[11]+x[12],
             x[13]+x[14],
             x[15]+x[16]) / 2

      if(direction=="upwind") p <- p[c(5:8, 1:4)]

      if(abs(g[3] - g[4]) == 1){ # standard non-suture edges
            if(g[1]<g[2] & g[4]<g[3]) return(p[1]) #SW
            if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
            if(g[2]<g[1] & g[4]<g[3]) return(p[3]) #NW
            if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
            if(g[2]<g[1] & g[3]<g[4]) return(p[5]) #NE
            if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
            if(g[1]<g[2] & g[3]<g[4]) return(p[7]) #SE
            if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
      }else{ # suture edges (only relevant if wrap == TRUE)
            if(g[1]<g[2] & g[4]>g[3]) return(p[1]) #SW
            if(g[1]==g[2] & g[4]>g[3]) return(p[2]) #W
            if(g[2]<g[1] & g[4]>g[3]) return(p[3]) #NW
            if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
            if(g[2]<g[1] & g[3]>g[4]) return(p[5]) #NE
            if(g[1]==g[2] & g[3]>g[4]) return(p[6]) #E
            if(g[1]<g[2] & g[3]>g[4]) return(p[7]) #SE
            if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
      }
}



#' A modified version of `gdistance::transition`,
#'
#' This version accepts a raster stack and a custom transition function; the
#' original does not support arbitrary transition functions for multi-layer
#' transition data.
transition_stack <- function(x, transitionFunction, directions, symm, wrap = FALSE, ...){

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


      ##### start modifications to gdistance::transition #####

      # add adjacencies between left and right edges
      if(wrap){
            i <- raster(x)
            i[] <- Cells
            i <- terra::as.matrix(i)
            nc <- ncol(i)
            nr <- nrow(i)
            s <- rbind(cbind(i[, 1], i[, nc]), # horizontal
                       cbind(i[2:nr, 1], i[1:(nr-1), nc]),
                       cbind(i[1:(nr-1), 1], i[2:nr, nc]))
            adj <- rbind(adj, s, s[,2:1])
      }

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





