\name{expokit_dgpadm_Qmat}
\alias{expokit_dgpadm_Qmat}
\title{EXPOKIT dgpadm matrix exponentiation on Q matrix}
\usage{
  expokit_dgpadm_Qmat(Qmat = NULL, t = 2.1,
    transpose_needed = TRUE)
}
\arguments{
  \item{Qmat}{an input Q transition matrix}

  \item{t}{one or more time values to exponentiate by}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal)}
}
\value{
  \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_}
  produces additional output relating to accuracy of the
  output matrix etc.; these can be obtained by a direct
  call of wrapalldmexpv_.
}
\description{
  This function exponentiates a matrix via the EXPOKIT padm
  function (designed for small dense matrices) and wrapper
  function \code{wrapalldgpadm_} around dmexpv.\cr \cr From
  EXPOKIT:\cr \cr * Computes exp(t*H), the matrix
  exponential of a general matrix in \cr * full, using the
  irreducible rational Pade approximation to the \cr *
  exponential function exp(x) = r(x) = (+/-)( I +
  2*(q(x)/p(x)) ), \cr * combined with
  scaling-and-squaring.  \cr \cr If \code{Qmat} is NULL
  (default), a default matrix is input.\cr
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{mat2coo}}
}

