#' A wind transition function to supply to `transition_stack`.
#'
#' @param x as used internally within `transition_stack`, this is a vector of
#'   collated windrose values for FROM and TO cells
#' @param direction either "downwind" (the default) or "upwind", indicating
#'   whether outbound or inbound wind conductance should be computed.
windflow <- function(x, direction="downwind"){

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

      if(g[1]<g[2] & g[4]<g[3]) return(p[1]) #SW
      if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
      if(g[2]<g[1] & g[4]<g[3]) return(p[3]) #NW
      if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
      if(g[2]<g[1] & g[3]<g[4]) return(p[5]) #NE
      if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
      if(g[1]<g[2] & g[3]<g[4]) return(p[7]) #SE
      if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
}
