


#' Run SDA4D Bayesian Sparse Tensor Decomposition with default parameters.
#'
#' RunSDA4D performs a tensor decomposition of a 4D tensor with the default
#' parameters as used in the accompanying paper.
#'
#' @param data_tensor four dimensional numeric tensor with dimensions $(N,L,M,T)$
#'  where $N$ is the number of samples, $L$ is the number of genes.
#' @param dimn_vector numeric vector of length 4 which determines the dimension
#'  of the tensor.
#' @param num_components integer number of components to run the decomposition 
#' with. Default is 4, minimum is 3.
#' @param max_iters integer maximum number of iterations to run. Defaults to
#'  2000.
#' @param stopping binary TRUE/FALSE. Defaults to TRUE. If TRUE then method will stop either
#' after \code{maxiters} or when the average number of PIPs passing the
#' threshold of 0.5 over \code{track} iterations drops below 1.
#' @param track integer number of iterations to keep track of if implementing
#' the stopping criterion as described in \code{stopping} description.
#'   Defaults to 10.
#'
#' @return a list of inferred parameters and ELBO trace.
#'
#'
#' @examples 
#' set.seed(531)
#' ex_data <- generate_example_data()$dataTensor
#' result <- RunSDA4D(ex_data, dim(ex_data))
#' 
#' @export
RunSDA4D<-function(data_tensor,
                   dimn_vector,
                   num_components=4,
                   max_iters=2000,
                   stopping=TRUE,
                   track=10){
    stopifnot(length(dim(data_tensor))==4)
    stopifnot(all(dimn_vector>0))
    stopifnot(num_components>=3)
    if(any(is.na(data_tensor))){
      stop('ERROR: the data tensor contains missing data')
    }
    if(any(is.infinite(data_tensor))){
      stop('ERROR: the data tensor contains infinite data')
    }
    if(!is.numeric(data_tensor)){
      stop('ERROR: the data tensor contains non-numeric data')
    }
    stopifnot(length(dimn_vector)==length(dim(data_tensor)))
    if(!all(dim(data_tensor)==dimn_vector)){
        stop('ERROR: dimension of tensor does not match provided dimensions,
             stopping')
    }
    
    if(dimn_vector[3]>1 & dimn_vector[4]==1){
      stop("ERROR: The 3D model applies to tensors of dimensions N,L,1,T, 
           check the dimensions of your tensor")
    }
    if(dimn_vector[1]<=1 || dimn_vector[2]<=1){
      stop("ERROR: The first two dimensions N,L of the data tensor must be
           larger than 1.")
    }
    
    stopifnot(num_components>1)
    stopifnot(max_iters>1)
    if(stopping & track<1){
        stop('ERROR: stopping is set to TRUE but track is < 1 , set track to at
             least 1.')
    }

    ## next, setup the default hyperparameters and dimension vectors.
    ## These are labelled identically to the details in the paper.
    params_to_run=list(N=dimn_vector[1],
                       L=dimn_vector[2],
                       M=dimn_vector[3],
                       T=dimn_vector[4],
                       C=num_components,
                       e=1e-6,
                       f=1e6,
                       g=0,
                       h=0,
                       u=1e-6,
                       v=1e6,
                       r=1,
                       z=1)

    #check and report dimensions
    if(params_to_run$M==1 && params_to_run$T > 1){
        print('Dimensions M = 1, T>1 so running as 3D tensor decomposition')
    }else if(params_to_run$M==1 && params_to_run$T == 1){
      print('Dimensions M = 1, T = 1 so running as 2D factor analysis decomposition')
    }else{
        print('Running as 4D tensor decomposition')
        print(paste('with dimensions ',dimn_vector))
    }

    ## run the method and return variables
    res<-SparsePARAFAC(params_to_run,data_tensor,maxiter=max_iters,
                       stopping=stopping,track=track,debugging=FALSE)
    if(res$Error==1){
        print('Finished but Negative Free Energy shows signs of a decrease,
              this should be checked.')
    }else{
        print('Finished')
    }
    res[[1]]<-NULL #removes the error indicator from the output.
    if(stopping){
        print(paste(res$maximumiteration,' Iterations were carried out.'))
    }
    names(res)[1] = "ELBO" # renames Neg_FE to ELBO for the vignette.
    if(res$ELBO[1]!=0){
        res$ELBO = res$ELBO[res$ELBO!=0]
    }
    print('returning output')
    return(res)
}
