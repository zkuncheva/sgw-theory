#' get supra adjacency
#' 
#' A function to get supra adjacency from list of adjacency matrices
#' @param G list of graphs
#' @param periodic whether to have periodic supra adj
#' @param omega inter-layer weight
#' @export
get_supra_adj = function(G, periodic = TRUE, omega = 1){
  # TODO : add LART option for weights
  # get block diagonal matrix
  n = dim(G[[1]])[1]
  adjG = bdiag(G)
  # get inter-layer weights on diagonals
  adjG[row(adjG) + n == col(adjG)] = omega
  adjG[col(adjG) + n == row(adjG)] = omega
  if (periodic){
    tt=length(G)
    adjG[row(adjG) + (tt-1)*n == col(adjG)] = omega
    adjG[col(adjG) + (tt-1)*n == row(adjG)] = omega
  }
  return(adjG)
}
