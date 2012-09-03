\name{expokit_wrapalldmexpv_tvals}
\alias{expokit_wrapalldmexpv_tvals}
\title{Run EXPOKIT's dmexpv}
\usage{
  expokit_wrapalldmexpv_tvals(Qmat = NULL, tvals = c(2.1),
    transpose_needed = TRUE, COO_needed = TRUE,
    coo_n = NULL, force_list_if_1_tval = FALSE)
}
\arguments{
  \item{Qmat}{an input Q transition matrix}

  \item{tvals}{one or more time values to exponentiate by
  (doesn't have to literally be a time value, obviously)}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal)}

  \item{COO_needed}{Should the matrix be tranposed to COO?
  COO format is required for EXPOKIT's sparse-matrix
  functions (like dmexpv and unlike the padm-related
  functions. Default TRUE; if false, user must put a
  COO-formated matrix in \code{Qmat}.}

  \item{coo_n}{If a COO matrix is input, \code{coo_n}
  specified the order of the matrix.}

  \item{force_list_if_1_tval}{Default FALSE, but set to
  TRUE if you want a single matrix to be returned inside a
  list}
}
\value{
  \code{tmpoutmat} the output matrix, if 1 t-value is
  input; \code{list_of_matrices_output}, if more than 1
  t-value is input; to get a single output matrix in a
  list, set \code{force_list_if_1_tval=TRUE}
}
\description{
  The function runs EXPOKIT's \code{dmexpv} function on a Q
  matrix and \emph{one or more} time values.  If
  \code{Qmat} is NULL (default), a default matrix is input.
}
\examples{
# Make a square instantaneous rate matrix (Q matrix)
# This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
# to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
# Unleashed" at:
# http://www.bioinf.org/molsys/data/idiots.pdf
#
# The Q matrix includes the stationary base freqencies, which Pmat
# converges to as t becomes large.
Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)

# Make a series of t values
tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)

# DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
# DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.

# DGEXPV, single t-value
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)

# This function runs the for-loop itself (sadly, we could not get mapply() to work
# on a function that calls dmexpv/dgexpv), returning a list of probability matrices.

# DGEXPV functions
list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals, transpose_needed=TRUE)
list_of_P_matrices_dgexpv
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_dmexpv_wrapper}}

  \code{\link{expokit_dmexpv_Qmat}}
}

