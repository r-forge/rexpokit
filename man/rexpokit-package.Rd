\docType{package}
\name{rexpokit-package}
\alias{rexpokit}
\alias{rexpokit-package}
\title{Matrix exponentiation with EXPOKIT in R}
\description{
  Matrix exponentiation with EXPOKIT in R
}
\details{
  \tabular{ll}{ Package: \tab rexpokit\cr Type: \tab
  Package\cr Version: \tab 0.21\cr Date: \tab 2013-01-30\cr
  License: \tab GPL (>= 2)\cr LazyLoad: \tab yes\cr }

  This package wraps some of the matrix exponentiation
  utilities from EXPOKIT
  (\url{http://www.maths.uq.edu.au/expokit/}), a FORTRAN
  library that is widely recommended for fast matrix
  exponentiation (Sidje RB, 1998. "Expokit: A Software
  Package for Computing Matrix Exponentials." \emph{ACM
  Trans. Math. Softw.} 24(1): 130-156).

  EXPOKIT includes functions for exponentiating both small,
  dense matrices, and large, sparse matrices (in sparse
  matrices, most of the cells have value 0).  Rapid matrix
  exponentiation is useful in phylogenetics when we have a
  large number of states (as we do when we are inferring
  the history of transitions between the possible
  geographic ranges of a species), but is probably useful
  in other ways as well.

  Various messages on discussion boards have asked whether
  or not there is an R package that uses EXPOKIT.  There
  are only two as of this writing (January 2013) --
  \code{\link[diversitree:find.mle]{diversitree}} and
  \code{\link[ctarma]{ctarma}}. However, diversitree's
  usage is nested deeply in a series of dynamic functions
  and integrated with additional libraries (e.g. deSolve)
  and so is very difficult to extract for general usage,
  and \code{\link[ctarma]{ctarma}} implements only ZEXPM
  via
  \code{\link[ctarma]{ctarma}}::\code{\link[ctarma:zexpm]{zexpm}}.
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk} is
  also working on an implementation of certain EXPOKIT
  functions.

  (See the additional notes file
  EXPOKIT_For_Dummies_notes_v1.txt for additional notes on
  wrappers for EXPOKIT in Python etc.)

  As it turns out, the EXPOKIT documentation and code is
  far from trivial to figure out, since the code as
  published does not run "out of the box" -- in particular,
  the Q transition matrix ("matvec"), which is the major
  input into an exponentiation algorithm, is not input
  directly, but rather via another function, which requires
  the user to hack together some FORTRAN code to do this
  and make a wrapper for the core function.  I couldn't
  figure it out in a short amount of time, but Stephen
  Smith did for his "Lagrange" biogeography package, so I
  essentially copied this chunk of his code to get started.
  \cr \cr \bold{Installation hints}\cr \cr Installing
  \code{rexpokit} from source will require a gfortran
  compiler to convert the FORTRAN code files in /src (*.f)
  to object files (*.o), and g++ to compile and link the
  C++ wrapper. \code{rexpokit} was developed on an Intel
  Mac running OS X 10.7.  I (NJM) successfully compiled it
  using g++ and gfortran from (gcc version 4.2.1).\cr \cr
  \cr \bold{Acknowledgements/sources}\cr \cr \bold{1.}
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}
  helped greatly with the initial setup of the package. \cr
  \bold{2.} EXPOKIT, original FORTRAN package, by Roger B.
  Sidje \email{rbs@maths.uq.edu.au}, Department of
  Mathematics, University of Queensland, Brisbane,
  QLD-4072, Australia, (c) 1996-2006 All Rights Reserved\cr
  \cr Sidje has given permission for us to include EXPOKIT
  code in this R package under the usual GPL license for
  CRAN R packages. For the full EXPOKIT copyright and
  license, see expokit_copyright.txt \cr

  EXPOKIT was published by Sidje in: Sidje RB (1998).
  "Expokit. A Software Package for Computing Matrix
  Exponentials." \emph{ACM-Transactions on Mathematical
  Software}, 24(1):130-156.
  \url{http://tinyurl.com/bwa87rq}\cr \cr \bold{3.} C++
  code wrapping EXPOKIT was copied in part from a file in
  LAGRANGE, C++ version by Stephen Smith:\cr
  \url{http://code.google.com/p/lagrange/}\cr
  \url{https://github.com/blackrim/lagrange}\cr \cr
  Specifically:\cr \cr \code{ * RateMatrix.cpp}\cr \code{ *
  }\cr \code{ * Created on: Aug 14, 2009}\cr \code{ *
  Author: smitty}\cr \code{ *}\cr \cr ...and the my_*.f
  wrappers for the EXPOKIT *.f code files.\cr \cr \bold{4.}
  Also copied in part (to get the .h file) from:\cr \cr
  Python package "Pyprop":\cr
  \url{http://code.google.com/p/pyprop/}\cr
  \url{http://pyprop.googlecode.com/svn/trunk/core/krylov/expokit/expokitpropagator.cpp}\cr
  \url{http://www.koders.com/python/fidCA95B5A4B2FB77455A72B8A361CF684FFE48F4DC.aspx?s=fourier+transform}\cr
  \cr Specifically:\cr
  pyprop/core/krylov/expokit/f2c/expokit.h \cr \cr
  \bold{5.} The EXPOKIT FORTRAN package is available at:\cr
  \url{http://www.maths.uq.edu.au/expokit/}\cr \cr
  Copyright:\cr
  \url{http://www.maths.uq.edu.au/expokit/copyright}\cr
  ...or...\cr expokit_copyright.txt in this install\cr
}
\examples{
# Example code
# For background, see EXPOKIT_For_Dummies_notes_v1.txt

library(rexpokit)

# Make a square instantaneous rate matrix (Q matrix)
# This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
# to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
# Unleashed" at:
# \\url{http://www.bioinf.org/molsys/data/idiots.pdf}
#
# The Q matrix includes the stationary base freqencies, which Pmat
# converges to as t becomes large.
Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)

# Make a series of t values
tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)

# Exponentiate each with EXPOKIT's dgpadm (good for small dense matrices)
for (t in tvals)
	{
	Pmat = expokit_dgpadm_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
	cat("\\n\\nTime=", t, "\\n", sep="")
	print(Pmat)
	}

# Exponentiate each with EXPOKIT's dmexpv (should be fast for large sparse matrices)
for (t in tvals)
	{
	Pmat = expokit_dmexpv_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
	cat("\\n\\nTime=", t, "\\n", sep="")
	print(Pmat)
	}

# DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
# DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.

# DMEXPV, single t-value
expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)

# DGEXPV, single t-value
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)

# These functions runs the for-loop itself (sadly, we could not get mapply() to work
# on a function that calls dmexpv/dgexpv), returning a list of probability matrices.

# DMEXPV functions
list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals, transpose_needed=TRUE)
list_of_P_matrices_dmexpv

# DGEXPV functions
list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals, transpose_needed=TRUE)
list_of_P_matrices_dgexpv

# Check if there are differences in the results (might only happen for large problems)
cat("\\n")
cat("Differences between dmexpv and dgexpv\\n")

for (i in 1:length(list_of_P_matrices_dmexpv))
	{
	diffs = list_of_P_matrices_dmexpv[[i]] - list_of_P_matrices_dgexpv[[i]]
	print(diffs)
	cat("\\n")
	}
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://www.maths.uq.edu.au/expokit/} \cr
  \url{http://www.maths.uq.edu.au/expokit/copyright}

  Sidje RB (1998). "Expokit. A Software Package for
  Computing Matrix Exponentials." _ACM Trans. Math.
  Softw._, *24*(1), pp. 130-156.

  Moler C and Loan CV (2003). "Nineteen dubious ways to
  compute the exponential of a matrix, twenty-five years
  later." _SIAM review_, pp. 3-49.

  Foster PG (2001). "The Idiot's Guide to the Zen of
  Likelihood in a Nutshell in Seven Days for Dummies,
  Unleashed." Online PDF, widely copied, <URL:
  \url{http://www.bioinf.org/molsys/data/idiots.pdf}>.
}
\seealso{
  \code{\link{expokit_wrapalldmexpv_tvals}}
}
\keyword{expokit}
\keyword{exponentiation,}
\keyword{matrix}
\keyword{matrix,}
\keyword{package,}
\keyword{phylogenetics,}
\keyword{transition}

