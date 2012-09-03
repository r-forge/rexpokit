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
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_dmexpv_wrapper}}

  \code{\link{expokit_dmexpv_Qmat}}
}

