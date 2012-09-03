\name{SparseM_coo_to_REXPOKIT_coo}
\alias{SparseM_coo_to_REXPOKIT_coo}
\title{Convert a SparseM COO matrix to a plain matrix}
\usage{
  SparseM_coo_to_REXPOKIT_coo(tmpmat_in_SparseMcoo_fmt)
}
\arguments{
  \item{tmpmat_in_SparseMcoo_fmt}{A square matrix S4 object
  derived from SparseM's as.matrix.coo}
}
\value{
  tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia},
  \code{ja}, and \code{a}
}
\description{
  Converts a SparseM COO-formatted matrix (an S4 object) to
  a plain matrix, with \cr column #1 = ia = i index\cr
  column #2 = ja = j index\cr column #3 = a = nonzero
  values of the matrix\cr
}
\details{
  Background: COO (coordinated list) format, is described
  here:\cr \cr
  \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
  \cr In \code{EXPOKIT} and its wrapper functions, a
  COO-formated matrix is input as 3 vectors (first two
  integer, the third double):\cr \cr ia = row number\cr ja
  = column number\cr a = value of that cell in the matrix
  (skipping 0 cells)\cr
}
\examples{
# Example use:
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{mat2coo_forloop}}
}

