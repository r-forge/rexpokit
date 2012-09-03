\name{expokit_dgexpv_Qmat}
\alias{expokit_dgexpv_Qmat}
\title{EXPOKIT dgexpv matrix exponentiation on Q matrix}
\usage{
  expokit_dgexpv_Qmat(Qmat = NULL, t = 2.1,
    transpose_needed = TRUE, anorm = NULL)
}
\arguments{
  \item{Qmat}{an input Q transition matrix}

  \item{t}{a time value to exponentiate by}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal)}

  \item{anorm}{\code{dgexpv} requires an initial guess at
  the norm of the matrix. Using the R function
  \code{\link{norm}} might get slow with large matrices. If
  so, the user can input a guess manually (\code{Lagrange}
  seems to just use 1 or 0, if I recall correctly).}
}
\value{
  \code{tmpoutmat} the output matrix. \code{wrapalldgexpv_}
  produces additional output relating to accuracy of the
  output matrix etc.; these can be by a direct call of
  dgexpv.
}
\description{
  This function converts a matrix to COO format and
  exponentiates it via the EXPOKIT dgexpv function
  (designed for sparse matrices) and wrapper functions
  \code{wrapalldgexpv_} around dgexpv.\cr \cr NOTE: DGEXPV
  vs. DMEXPV. According to the EXPOKIT documentation,
  DGEXPV should be faster than DMEXPV, however DMEXPV runs
  an accuracy check appropriate for Markov chains, which is
  not done in DGEXPV.\cr \cr From EXPOKIT:\cr \cr * The
  method used is based on Krylov subspace projection\cr *
  techniques and the matrix under consideration interacts
  only\cr * via the external routine `matvec' performing
  the matrix-vector \cr * product (matrix-free method).\cr
  *\cr * This [DMEXPV, not DGEXPV -- NJM] is a customised
  version for Markov Chains. This means that a\cr * check
  is done within this code to ensure that the resulting
  vector \cr * w is a probability vector, i.e., w must have
  all its components \cr * in [0,1], with sum equal to 1.
  This check is done at some expense\cr * and the user may
  try DGEXPV which is cheaper since it ignores \cr *
  probability constraints.\cr \cr I (NJM) have not noticed
  a difference between the outputs of these two functions,
  but it might occur with large matrices.
}
\details{
  COO (coordinated list) format is a compressed format that
  is required for EXPOKIT's sparse-matrix functions (like
  dgexpv and unlike EXPOKIT's padm-related functions. COO
  format is described here:\cr \cr
  \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
  \cr If \code{Qmat} is NULL (default), a default matrix is
  input.\cr
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{mat2coo}}

  \code{\link{expokit_dgexpv_wrapper}}

  \code{\link{expokit_wrapalldgexpv_tvals}}
}

