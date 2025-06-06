#' conn_gnp
#' 
#' A function to simulate multilayered slowly changing Erdos-Renyi network
#' @param p edge probability
#' @param n numebr of nodes per layer
#' @param max.it maximum number of iterations to try to get a connected network
#' @export
conn_gnp = function(n, p, max.it = 10){
  bool = FALSE
  iter = 1
  while (!bool){
    g=sample_gnp(n,p)
    bool = is_connected(g)
    if (iter==max.it){
      stop('A connected graph could not be produced for max iterations.')
    }else if(!bool){
      iter = iter+1
    }
  }
  return(g)
}