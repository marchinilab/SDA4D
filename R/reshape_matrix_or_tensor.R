

#' matrixToTensor converts 2D tensor to 4D  array suitable for SDA4D.
#' 
#' matrixToTensor converts a matrix of dimension N by L to a tensor
#'  of dimension N by L by 1 by 1, as is required to run SDA4D. 
#' 
#' @param input_matrix numeric matrix
#' @returns four dimensional tensor with entries from input_matrix.
#' @examples
#' matrixToTensor(matrix(rnorm(50),25,2))
#' @export
matrixToTensor<- function(input_matrix){
  d_vec = dim(input_matrix)
  stopifnot(length(d_vec)==2)
  Y = array(c(input_matrix),dim = c(d_vec,1,1))
  return(Y)
}