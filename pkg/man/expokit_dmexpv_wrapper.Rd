\name{expokit_dmexpv_wrapper}
\alias{expokit_dmexpv_wrapper}
\title{EXPOKIT dmexpv wrapper function}
\usage{
  expokit_dmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp,
    lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
}
\arguments{
  \item{n}{number of rows in Q matrix}

  \item{m}{n-1}

  \item{t}{the value to exponentiate the rate matrix by
  (often e.g. a time value)}

  \item{v}{variable to store some results in; should have n
  elements (and perhaps start with 1)}

  \item{w}{same length as v}

  \item{tol}{tolerance for approximations; usually set to
  0.01}

  \item{anorm}{the norm of the Q matrix}

  \item{lwsp}{length of workspace (wsp); for dmexpv,
  lwsp=n*(m+2)+5*(m+2)^2+ideg+1}

  \item{wsp}{workspace to store some results in; should be
  a double with lwsp elements}

  \item{liwsp}{length of integer workspace; for dmexpv,
  liwsp=m+2}

  \item{iwsp}{integer workspace}

  \item{itrace}{option, set to 0}

  \item{iflag}{option, set to 0}

  \item{ia}{i indices of Qmat nonzero values}

  \item{ja}{j indices of Qmat nonzero values}

  \item{a}{nonzero values of Qmat (ia, ja, a are columns of
  a COO-formatted Q matrix)}

  \item{nz}{number of non-zeros in Qmat}

  \item{res}{space for output probability matrix (n x n)

  EXPOKIT needs the input matrix to be transposed compared
  to normal) COO format is required for EXPOKIT.}
}
\value{
  \code{tmpoutmat} the output matrix for the (first) input
  t-value
}
\description{
  This function wraps the .C call to EXPOKIT for the dmexpv
  function.  Only the output probability matrix is
  returned.
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_dmexpv_Qmat}}

  \code{\link{expokit_wrapalldmexpv_tvals}}
}

