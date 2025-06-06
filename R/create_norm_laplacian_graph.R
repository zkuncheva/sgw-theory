#' create_norm_laplacian_graph
#' 
#' A function to obtain normalized laplacian
#' @param A adjacency/supra adjacency matrix
#' @export

create_norm_laplacian_graph = function(A){
  if(!isSymmetric(A)){
    stop('Provided matrix A should be symmetric.')
  }
  n = dim(A)[1]
  D = matrix(0, n, n)
  diag(D) = A %*% rep(1,n)
  diag(D) = 1/sqrt(diag(D))
  Q = D %*% A %*% D  
  Ln = diag(1,n)-Q
  Ln[is.nan(Ln)] = 0
  return(Ln)
  
}