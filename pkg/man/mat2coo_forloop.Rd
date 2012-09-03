\name{mat2coo_forloop}
\alias{mat2coo_forloop}
\title{Convert matrix to COO format using nested for-loops}
\usage{
  mat2coo_forloop(tmpmat)
}
\arguments{
  \item{tmpmat}{A square matrix}
}
\value{
  tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia},
  \code{ja}, and \code{a}
}
\description{
  Converts a matrix to COO format. This version of the
  function uses for-loops, which is presumably less
  efficient than \code{\link{mat2coo}}.
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{mat2coo}}
}

