#' Matrix exponentiation with EXPOKIT in R
#'
#' \tabular{ll}{
#' Package: \tab rexpokit\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2012-06-13\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' This package wraps some of the matrix exponentiation
#' utilities from EXPOKIT
#' (\url{http://www.maths.uq.edu.au/expokit/}), a FORTRAN
#' library that is widely recommended as the best/fastest way
#' to exponentiate a matrix, especially a large, sparse matrix.
#' 
#' This is useful in phylogenetics when we have a large number
#' of states (as we do when we are inferring the history of
#' transitions between the possible geographic ranges of a
#' species), but is probably useful in other ways as well.
#' 
#' Various messages on discussion boards have asked whether or
#' not there is an R package that uses EXPOKIT, but I haven't
#' found one.  (But see EXPOKIT_For_Dummies_notes_v1.txt for
#' the various tidbits I did find.)
#' 
#' As it turns out, the EXPOKIT documentation and code is far
#' from trivial to figure out, since the code as published does
#' not run "out of the box" -- in particular, the Q transition
#' matrix ("matvec"), which is the major input into an
#' exponentiation algorithm, is not input directly, but rather
#' via another function, which requires the user to hack
#' together some FORTRAN code to do this and make a wrapper for
#' the core function.  I couldn't figure it out in a short
#' amount of time, but Stephen Smith did for his "Lagrange"
#' biogeography package, so I essentially copied this chunk of
#' his code.
#' 
#' The package also contains some example functions to simulate
#' dense and sparse transition matrices to test the speed of
#' the various algorithms, so that users can pick the best
#' choice for their purposes.  This will of course require
#' users to install the relevant R packages (I won't make those
#' required dependencies, however).\cr
#' \cr
#' Acknowledgements/sources:\cr
#' \cr
#' 1. Copied in part from a file in Lagrange, C++ version by Stephen Smith:\cr
#' \url{http://code.google.com/p/lagrange/}\cr
#' \url{https://github.com/blackrim/lagrange}\cr
#' \cr
#' Specifically:\cr
#'  * RateMatrix.cpp\cr
#'  * \cr
#'  *  Created on: Aug 14, 2009\cr
#'  *      Author: smitty\cr
#'  *\cr
#' ...and the my_*.f wrappers for the EXPOKIT *.f code files.\cr
#' \cr
#' 2. Also copied in part (to get the .h file) from:\cr
#' \cr
#' Python package "Pyprop":\cr
#' \url{http://code.google.com/p/pyprop/}\cr
#' \url{http://pyprop.googlecode.com/svn/trunk/core/krylov/expokit/expokitpropagator.cpp}\cr
#' \url{http://www.koders.com/python/fidCA95B5A4B2FB77455A72B8A361CF684FFE48F4DC.aspx?s=fourier+transform}\cr
#' \cr
#' Specifically:\cr
#' pyprop/core/krylov/expokit/f2c/expokit.h \cr
#' \cr
#' 3. EXPOKIT package is available at:\cr
#' \url{http://www.maths.uq.edu.au/expokit/}\cr
#' \cr
#' Copyright:\cr
#' \url{http://www.maths.uq.edu.au/expokit/copyright}\cr
#' ...or...\cr
#' expokit_copyright.txt in this install\cr
#' 
#'
#' @name rexpokit-package
#' @aliases rexpokit
#' @docType package
#' @title Matrix exponentiation with EXPOKIT in R
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @references
#' \url{http://www.maths.uq.edu.au/expokit/}
#' \url{http://www.maths.uq.edu.au/expokit/copyright}
#' @bibliography /Dropbox/_njm/__packages/rexpokit_setup/rexpokit_refs.bib
#'   @cite FosterIdiots
#'   @cite moler2003nineteen
#'   @cite Sidje1998
#' @keywords package, matrix, matrix exponentiation, phylogenetics, transition matrix, expokit
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @examples # Example code
#' # For background, see EXPOKIT_For_Dummies_notes_v1.txt
#' # For installation hints, see notes (?)
#' 
#' library(rexpokit)
#' 
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # http://www.bioinf.org/molsys/data/idiots.pdf
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # Exponentiate each with EXPOKIT's dgpadm (good for small dense matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dgpadm_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#' 
#' # Exponentiate each with EXPOKIT's dmexpv (should be fast for large sparse matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dmexpv_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#'
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # These functions runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DMEXPV functions
#' list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dmexpv
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#' 
#' # Check if there are differences in the results (might only happen for large problems)
#' cat("\n")
#' cat("Differences between dmexpv and dgexpv\n")
#' 
#' for (i in 1:length(list_of_P_matrices_dmexpv))
#' 	{
#' 	diffs = list_of_P_matrices_dmexpv[[i]] - list_of_P_matrices_dgexpv[[i]]
#' 	print(diffs)
#' 	cat("\n")
#' 	}
#' 

require(roxygen2)