\name{expokit_dmexpv_Qmat}
\alias{expokit_dmexpv_Qmat}
\title{EXPOKIT dmexpv matrix exponentiation on Q matrix}
\usage{
  expokit_dmexpv_Qmat(Qmat = NULL, t = 2.1,
    transpose_needed = TRUE, anorm = NULL)
}
\arguments{
  \item{Qmat}{an input Q transition matrix}

  \item{t}{one or more time values to exponentiate by}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal)}

  \item{anorm}{\code{dmexpv} requires an initial guess at
  the norm of the matrix. Using the R function
  \code{\link{norm}} might get slow with large matrices. If
  so, the user can input a guess manually (\code{Lagrange}
  seems to just use 1 or 0, if I recall correctly).}
}
\value{
  \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_}
  produces additional output relating to accuracy of the
  output matrix etc.; these can be by a direct call of
  dmexpv.
}
\description{
  This function converts a matrix to COO format and
  exponentiates it via the EXPOKIT dmexpv function
  (designed for sparse matrices) and wrapper functions
  \code{wrapalldmexpv_} around dmexpv. \cr From EXPOKIT:\cr
  * The method used is based on Krylov subspace
  projection\cr * techniques and the matrix under
  consideration interacts only\cr * via the external
  routine `matvec' performing the matrix-vector \cr *
  product (matrix-free method).\cr *\cr * This is a
  customised version for Markov Chains. This means that
  a\cr * check is done within this code to ensure that the
  resulting vector \cr * w is a probability vector, i.e., w
  must have all its components \cr * in [0,1], with sum
  equal to 1. This check is done at some expense\cr * and
  the user may try DGEXPV which is cheaper since it ignores
  \cr * probability constraints.\cr \cr COO (coordinated
  list) format is a compressed format that is\cr required
  for EXPOKIT's sparse-matrix functions (like dmexpv and\cr
  unlike EXPOKIT's padm-related functions.\cr \cr COO
  (coordinated list) format is described here:\cr \cr
  \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
  \cr If \code{Qmat} is NULL (default), a default matrix is
  input.\cr
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{mat2coo}}

  \code{\link{expokit_dmexpv_wrapper}}

  \code{\link{expokit_wrapalldmexpv_tvals}}
}

