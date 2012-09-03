#' @include fermat.R
#' @include rexpokit-package.R
require(roxygen2)
#sourcedir = '/Dropbox/_njm/'
#source3 = '_genericR_v1.R'
#source(paste(sourcedir, source3, sep=""))
#roxygenize()

# Re-source this R code after editing, without reinstalling from scratch:
# sourcedir = "/Dropbox/_njm/__packages/rexpokit_setup/"
# source8 = 'rexpokit_v1.R'
# source(paste(sourcedir, source8, sep=""))



# Original source:
# 
#sourcedir = "/Dropbox/_njm/"
#source8 = '_matrix_utils_v1.R'
#source(paste(sourcedir, source8, sep=""))

# for e.g. calc_loglike
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))


#######################################################
# EXPOKIT-RELATED FUNCTIONS
#######################################################

#' EXPOKIT dgpadm matrix exponentiation on Q matrix
#'
#' This function exponentiates a matrix via the EXPOKIT padm function
#' (designed for small dense matrices) and wrapper function 
#' \code{wrapalldgpadm_} around dmexpv.\cr
#'\cr
#' From EXPOKIT:\cr
#' \cr
#' *     Computes exp(t*H), the matrix exponential of a general matrix in \cr
#' *     full, using the irreducible rational Pade approximation to the   \cr
#' *     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), \cr
#' *     combined with scaling-and-squaring.                              \cr
#' \cr
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t one or more time values to exponentiate by
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' obtained by a direct call of wrapalldmexpv_.
#' @seealso \code{\link{mat2coo}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_dgpadm_Qmat <- function(Qmat=NULL, t=2.1, transpose_needed=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	

	
	# FOR DGPADM
	# ideg = 
	# "(input) the degree of the diagonal Pade to be used.
	# a value of 6 is generally satisfactory."
	ideg = as.integer(6)

	# Order (numrows/numcols) of the matrix
	# "(input) order of H."
	m = as.integer(nrow(Qmat))

	# output matrix
	res = double(length=m*m)

	# Prepare input matrix
	matvec = Qmat
	if (transpose_needed == TRUE)
		{
		tmatvec = t(matvec)
		H = as.numeric(tmatvec)
		} else {
		H = as.numeric(matvec)
		}
	
	# "H(ldh,m)  : (input) argument matrix."
	# (ldh = numrows and m is numcols, or something)
	ldh = m
	
	# No tolin PADM
	# tol or t?  should be t
	# tol = as.double(1)
	
	# lwsp = length of wsp, the workspace
	# "wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1."
	lwsp = as.integer(4*m*m+ideg+1)
	wsp = double(length=lwsp)
	
	# "ipiv(m)   : (workspace)"
	ipiv = integer(length=m)
	
	#  "iexph     : (output) number such that wsp(iexph) points to exp(tH)"
	#  "            i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)"
	iexph = as.integer(0)
	
	# "ns        : (output) number of scaling-squaring used."
	ns = as.integer(0)
	
	# "iflag     : (output) exit flag."
	# "	*                      0 - no problem"
	# "	*                     <0 - problem"
	iflag = as.integer(0)
	
	# Run the function:
	res <- .C("wrapdgpadm_", as.integer(ideg), as.integer(m), as.double(t), as.double(H), as.integer(ldh), as.double(wsp), as.integer(lwsp), as.integer(ipiv), as.integer(iexph), as.integer(ns), as.integer(iflag))
	
	output = res[[6]]
	output_Pmat_is = seq(res[[9]], res[[9]]+m*m-1, by=1)
	output_Pmat = output[output_Pmat_is]
	output_Pmat = matrix(output_Pmat, nrow=m, byrow=TRUE)
	#print(output_Pmat)
	
	return(output_Pmat)
	}







#' EXPOKIT dmexpv matrix exponentiation on Q matrix
#'
#' This function converts a matrix to COO format and exponentiates
#' it via the EXPOKIT dmexpv function (designed for sparse matrices)
#' and wrapper functions \code{wrapalldmexpv_} around dmexpv.
#'\cr
#' From EXPOKIT:\cr
#' *     The method used is based on Krylov subspace projection\cr
#' *     techniques and the matrix under consideration interacts only\cr
#' *     via the external routine `matvec' performing the matrix-vector \cr
#' *     product (matrix-free method).\cr
#' *\cr
#' *     This is a customised version for Markov Chains. This means that a\cr
#' *     check is done within this code to ensure that the resulting vector \cr
#' *     w is a probability vector, i.e., w must have all its components \cr
#' *     in [0,1], with sum equal to 1. This check is done at some expense\cr
#' *     and the user may try DGEXPV which is cheaper since it ignores \cr
#' *     probability constraints.\cr
#'\cr
#' COO (coordinated list) format is a compressed format that is\cr
#' required for EXPOKIT's sparse-matrix functions (like dmexpv and\cr
#' unlike EXPOKIT's padm-related functions.\cr
#'\cr
#' COO (coordinated list) format is described here:\cr
#'\cr
#' \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#' \cr
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t one or more time values to exponentiate by
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param anorm \code{dmexpv} requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. If so, the user
#' can input a guess manually (\code{Lagrange} seems to just use 1 or 0, if I
#' recall correctly).
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' by a direct call of dmexpv.
#' @seealso \code{\link{mat2coo}}
#' @seealso \code{\link{expokit_dmexpv_wrapper}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_dmexpv_Qmat <- function(Qmat=NULL, t=2.1, transpose_needed=TRUE, anorm=NULL)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	# number of non-zeros
	nz  = sum(Qmat != 0)

	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	n=nrow(Qmat)
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	v=double(n)
	v[1] = 1
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# lwsp = length of wsp
	# wsp = workspace to hold various variables, cells of the matrix, etc.
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = m+2
	iwsp = integer(length=liwsp)
	
	res = double(length=n*n)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = Qmat
	
	if (transpose_needed == TRUE)
		{
		tmatvec = t(matvec)
		}
	#rowSums(tmatvec)
	#colSums(tmatvec)
	
	# This might (?) get slow with large matrices -- doesn't seem to
	if (is.null(anorm))
		{
		# Use the 1-norm or one-norm
		anorm = as.numeric(norm(matvec, type="O"))
		}
	
	# The itrace flag, if set to 1, results in dmexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0	
	
	# Make the input COO matrix
	# COO = coordinate list format, useful for sparse matrices with lots of zeros:
	# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
	# ia = rownum in the matrix
	# ja = colnum in the matrix
	# a  = value of that cell
	tmpmat = tmatvec
	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	# Run the wrapper function	
	res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
		
	tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	#print(tmpoutmat)
	
	return(tmpoutmat)
	}



#' EXPOKIT dmexpv wrapper function
#'
#' This function wraps the .C call to EXPOKIT for the dmexpv function.  Only the output probability
#' matrix is returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dmexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dmexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#' @param res space for output probability matrix (n x n)
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' COO format is required for EXPOKIT.
#' @return \code{tmpoutmat} the output matrix for the (first) input t-value
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_dmexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	{
	res2 = NULL
	
	res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
	tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
	return(tmpoutmat)
	}

	


#' Run EXPOKIT's dmexpv
#'
#' The function runs EXPOKIT's \code{dmexpv} function on a Q matrix and \emph{one or more} time values.  If \code{Qmat} is NULL (default), a default matrix is input.
#'
#' @param Qmat an input Q transition matrix
#' @param tvals one or more time values to exponentiate by (doesn't have to literally be a time value, obviously)
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param COO_needed Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dmexpv and unlike the padm-related 
#' functions. Default TRUE; if false, user must put a COO-formated matrix in \code{Qmat}.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order of the matrix.
#' @param force_list_if_1_tval Default FALSE, but set to TRUE if you want a single matrix to be returned
#' inside a list
#' @return \code{tmpoutmat} the output matrix, if 1 t-value is input; \code{list_of_matrices_output},
#' if more than 1 t-value is input; to get a single output matrix in a list, set \code{force_list_if_1_tval=TRUE}
#' @seealso \code{\link{expokit_dmexpv_wrapper}}
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_wrapalldmexpv_tvals <- function(Qmat=NULL, tvals=c(2.1), transpose_needed=TRUE, COO_needed=TRUE, coo_n=NULL, force_list_if_1_tval=FALSE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	COO_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	if (is.null(tvals))
		{
		warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
		}

	
	if (COO_needed == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(Qmat != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)		
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		nz  = sum(Qmat$a != 0)

		# here
		}
	

	#######################################################
	ideg = as.integer(6)
	#######################################################
	n=nrow(Qmat)
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	v=double(n)
	v[1] = 1
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# length of wsp
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = m+2
	iwsp = integer(length=liwsp)
	
	res = double(length=n*n)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = Qmat
	tmatvec = t(matvec)
	rowSums(tmatvec)
	colSums(tmatvec)
	
	# type="O" is being used here, this is supposed to be the
	# default for norm(), although it throws an error if not
	# specified
	# 
	# From the help:
	# type - character string, specifying the type of matrix norm to be
	# computed. A character indicating the type of norm desired. 
	# 	"O", "o" or "1"
	# 		specifies the one norm, (maximum absolute column sum);
	anorm = as.numeric(norm(matvec, type="O"))
	#anorm = 1
	
	
	itrace = 0
	iflag = 0	
	
	
	#a = as.numeric(tmatvec)
	#a = as.numeric(matvec)
	tmpmat = tmatvec
	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	num_tvals = length(tvals)
	
	# If there is more than 1 t-value, or if the user desires a list even for a single
	# t-value, return a list
	if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
		{
		# Loop through the list of tvals, get the prob. matrix for each
		# sadly, mapply() etc. crash when tried on expokit_dmexpv_wrapper
		
		# Set up empty matrix
		NA_matrix = matrix(NA, nrow=n, ncol=n)
		
		# Set up list of empty matrices
		list_of_matrices_output = replicate(NA_matrix, n=num_tvals, simplify=FALSE)
		
		for (i in 1:num_tvals)
			{
			t = tvals[i]
			list_of_matrices_output[[i]] = expokit_dmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)

			} # end forloop
		
		return(list_of_matrices_output)
		
		} else {
		
		# If there is only 1 t value, just return 1 matrix
		#res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
		
		#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)

		t=tvals
		tmpoutmat = expokit_dmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)		
		
		#print(tmpoutmat)
		return(tmpoutmat)
		} # End if statment
	
	}



#' Convert matrix to COO format using SparseM function
#'
#' Converts a matrix to COO format using the SparseM function, presumably this
#' is faster than using a for-loop.\cr
#'\cr
#' \code{EXPOKIT}'s \code{dmexp}-type functions deal with sparse matrices.
#' These have a lot of zeros, and thus can be compressed
#' into COO (coordinated list) format, which is described here:\cr
#' \cr
#' \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'\cr
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'\cr
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#' 
#' @param tmpmat A square matrix
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo_forloop}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
mat2coo <- function(tmpmat)
	{
	defaults = '
	tmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
	'
	
	numrows = nrow(tmpmat)
	numcells = numrows ^2
	
	# require(sfsmisc)
	# xy.grid
	
	x = 1:numrows
	y = 1:numrows

	# Seems to be slow when numrow > 1000
	# as.vector(tmpmat) appends col1vals, then col2vals, etc., so ji = xy
	# cells_ij = expand.grid(x, y)
	# tmpa = as.vector(tmpmat)
	# 
	# # Remove 0s
	# TF = tmpa != 0
	# ia = cells_ij[,1][TF]	
	# ja = cells_ij[,2][TF]	
	# a = tmpa[TF]

	require(SparseM)	# required for the as.matrix.coo function
	
	# This produces a matrix in coo format
	# (this is an S4 object)
	tmpmat_in_SparseMcoo_fmt = as.matrix.coo(tmpmat)
	tmpmat_in_REXPOKIT_coo_fmt = SparseM_coo_to_REXPOKIT_coo(tmpmat_in_SparseMcoo_fmt)
	
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}



#' Convert a SparseM COO matrix to a plain matrix
#'
#' Converts a SparseM COO-formatted matrix (an S4 object) to a plain matrix, with \cr
#' column #1 = ia = i index\cr
#' column #2 = ja = j index\cr
#' column #3 = a = nonzero values of the matrix\cr
#'
#' Background: COO (coordinated list) format, is described here:\cr
#' \cr
#' \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'\cr
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'\cr
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#' 
#' @param tmpmat_in_SparseMcoo_fmt A square matrix S4 object derived from SparseM's as.matrix.coo
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo_forloop}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
SparseM_coo_to_REXPOKIT_coo <- function(tmpmat_in_SparseMcoo_fmt)
	{
	tmpcoo = tmpmat_in_SparseMcoo_fmt
	
	# We just need the 3 columns: i index, j index, and nonzero values
	tmpmat_in_REXPOKIT_coo_fmt = cbind(tmpcoo@ia, tmpcoo@ja, tmpcoo@ra)
	
	# Apply appropriate column names
	colnames(tmpmat_in_REXPOKIT_coo_fmt) = c("ia", "ja", "a")
	
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}




#' Convert a COO-formated matrix to standard square format
#'
#' \code{EXPOKIT}'s \code{dmexp}-type functions deal with sparse matrices.
#' These have a lot of zeros, and thus can be compressed
#' into COO (coordinated list) format, which is described here:\cr
#' \cr
#' \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'\cr
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'\cr
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#'\cr
#' This function takes a 3-column matrix or data.frame (basically \code{cbind(ia, ja, a)})
#' and the order of the matrix, \code{n} (n = the order of the matrix, i.e. number of
#' rows/columns) and converts back to standard square format.\cr
#'
#' @param coomat a 3-column matrix or data.frame (basically \code{cbind(ia, ja, a)})
#' @param n the order of the matrix
#' @return outmat
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
#' ia = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
#' ja = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)
#' a  = c(-1.218, 0.126, 0.168, 0.126, 0.504, -0.882, 0.504, 0.672, 0.336, 0.252, -1.050, 0.252, 0.378, 0.504, 0.378, -1.050)
#' coomat = cbind(ia, ja, a)
#' print(coomat)
#' n = 4
#' Qmat = coo2mat(coomat, n)
#' print(Qmat)
coo2mat <- function(coomat, n)
	{
	# Make an empty matrix of 0s
	outmat = matrix(double(length=n*n), nrow=n)
	
	# go through each row of coomat
	ia = coomat[,1]
	ja = coomat[,2]
	a = coomat[,3]
	
	for (k in 1:length(ia))
		{
		#cat(ia[k], ja[k], a[k], "\n")
		outmat[ia[k], ja[k]] = a[k]
		}
	
	return(outmat)
	}



#' Convert matrix to COO format using nested for-loops
#'
#' Converts a matrix to COO format. This version of the function uses
#' for-loops, which is presumably less efficient than \code{\link{mat2coo}}.
#'
#' @param tmpmat A square matrix
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
mat2coo_forloop <- function(tmpmat)
	{
	# Number of non-zeros
	nz = sum(tmpmat != 0)
	
	# Blank columns for COO matrix
	ia = integer(length=nz)
	ja = integer(length=nz)
	a = integer(length=nz)
	
	count = 0
	for (i in 1:nrow(tmpmat))
		{
		for (j in 1:ncol(tmpmat))
			{
			if (tmatvec[i,j] != 0)
				{
				count = count+1
				ia[count] = i
				ja[count] = j
				a[count] = tmpmat[i,j]
				}
			}
		}
	tmpmat_in_REXPOKIT_coo_fmt = cbind(ia, ja, a)
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}







#######################################################
# 
# NOTE: DGEXPV section.  Same code as dmexpv, but EXPOKIT's DGEXPV should be
# faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
# Markov chains, which is not done in DGEXPV.
#
#######################################################




#' EXPOKIT dgexpv matrix exponentiation on Q matrix
#'
#' This function converts a matrix to COO format and exponentiates
#' it via the EXPOKIT dgexpv function (designed for sparse matrices)
#' and wrapper functions \code{wrapalldgexpv_} around dgexpv.\cr
#'\cr
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.\cr
#' \cr
#' From EXPOKIT:\cr
#'\cr
#' *     The method used is based on Krylov subspace projection\cr
#' *     techniques and the matrix under consideration interacts only\cr
#' *     via the external routine `matvec' performing the matrix-vector \cr
#' *     product (matrix-free method).\cr
#' *\cr
#' *     This [DMEXPV, not DGEXPV -- NJM] is a customised version for Markov Chains. This means that a\cr
#' *     check is done within this code to ensure that the resulting vector \cr
#' *     w is a probability vector, i.e., w must have all its components \cr
#' *     in [0,1], with sum equal to 1. This check is done at some expense\cr
#' *     and the user may try DGEXPV which is cheaper since it ignores \cr
#' *     probability constraints.\cr
#'\cr
#' I (NJM) have not noticed a difference between the outputs of these two functions, but it might
#' occur with large matrices.
#'
#' COO (coordinated list) format is a compressed format that is
#' required for EXPOKIT's sparse-matrix functions (like dgexpv and
#' unlike EXPOKIT's padm-related functions. COO format is described here:\cr
#'\cr
#' \link{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#' \cr
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t a time value to exponentiate by
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param anorm \code{dgexpv} requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. If so, the user
#' can input a guess manually (\code{Lagrange} seems to just use 1 or 0, if I
#' recall correctly).
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldgexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' by a direct call of dgexpv.
#' @seealso \code{\link{mat2coo}}
#' @seealso \code{\link{expokit_dgexpv_wrapper}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_dgexpv_Qmat <- function(Qmat=NULL, t=2.1, transpose_needed=TRUE, anorm=NULL)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	# number of non-zeros
	nz  = sum(Qmat != 0)

	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	n=nrow(Qmat)
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	v=double(n)
	v[1] = 1
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# lwsp = length of wsp
	# wsp = workspace to hold various variables, cells of the matrix, etc.
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = m+2
	iwsp = integer(length=liwsp)
	
	res = double(length=n*n)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = Qmat
	
	if (transpose_needed == TRUE)
		{
		tmatvec = t(matvec)
		}
	#rowSums(tmatvec)
	#colSums(tmatvec)
	
	# This might (?) get slow with large matrices -- doesn't seem to
	if (is.null(anorm))
		{
		# Use the 1-norm or one-norm
		anorm = as.numeric(norm(matvec, type="O"))
		}
	
	# The itrace flag, if set to 1, results in dgexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0	
	
	# Make the input COO matrix
	# COO = coordinate list format, useful for sparse matrices with lots of zeros:
	# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
	# ia = rownum in the matrix
	# ja = colnum in the matrix
	# a  = value of that cell
	tmpmat = tmatvec
	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	# Run the wrapper function	
	res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
		
	tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	#print(tmpoutmat)
	
	return(tmpoutmat)
	}



#' EXPOKIT dgexpv wrapper function
#'
#' This function wraps the .C call to EXPOKIT for the dgexpv function.  Only the output probability
#' matrix is returned.
#'
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.
#' 
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dgexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dgexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#' @param res space for output probability matrix (n x n)
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' COO format is required for EXPOKIT.
#' @return \code{tmpoutmat} the output matrix for the (first) input t-value
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_dgexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	{
	res2 = NULL
	
	res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
	tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
	return(tmpoutmat)
	}

	


#' Run EXPOKIT's dgexpv on one or more t-values
#'
#' The function runs EXPOKIT's \code{dgexpv} function on a Q matrix and \emph{one or more} time values.  If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'\cr
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.\cr
#' 
#' @param Qmat an input Q transition matrix
#' @param tvals one or more time values to exponentiate by (doesn't have to literally be a time value, obviously)
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param COO_needed Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dgexpv and unlike the padm-related 
#' functions. Default TRUE; if false, user must put a COO-formated matrix in \code{Qmat}.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order of the matrix.
#' @param force_list_if_1_tval Default FALSE, but set to TRUE if you want a single matrix to be returned
#' inside a list
#' @return \code{tmpoutmat} the output matrix, if 1 t-value is input; \code{list_of_matrices_output},
#' if more than 1 t-value is input; to get a single output matrix in a list, set \code{force_list_if_1_tval=TRUE}
#' @seealso \code{\link{expokit_dgexpv_wrapper}}
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
expokit_wrapalldgexpv_tvals <- function(Qmat=NULL, tvals=c(2.1), transpose_needed=TRUE, COO_needed=TRUE, coo_n=NULL, force_list_if_1_tval=FALSE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	COO_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	if (is.null(tvals))
		{
		warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
		}

	
	if (COO_needed == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(Qmat != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)		
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		nz  = sum(Qmat$a != 0)

		# here
		}
	

	#######################################################
        ideg = as.integer(6)
	#######################################################
	n=nrow(Qmat)
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	v=double(n)
	v[1] = 1
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# length of wsp
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = m+2
	iwsp = integer(length=liwsp)
	
	res = double(length=n*n)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = Qmat
	tmatvec = t(matvec)
	rowSums(tmatvec)
	colSums(tmatvec)
	
	# type="O" is being used here, this is supposed to be the
	# default for norm(), although it throws an error if not
	# specified
	# 
	# From the help:
	# type - character string, specifying the type of matrix norm to be
	# computed. A character indicating the type of norm desired. 
	# 	"O", "o" or "1"
	# 		specifies the one norm, (maximum absolute column sum);
	anorm = as.numeric(norm(matvec, type="O"))
	#anorm = 1
	
	
	itrace = 0
	iflag = 0	
	
	
	#a = as.numeric(tmatvec)
	#a = as.numeric(matvec)
	tmpmat = tmatvec
	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	num_tvals = length(tvals)
	
	# If there is more than 1 t-value, or if the user desires a list even for a single
	# t-value, return a list
	if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
		{
		# Loop through the list of tvals, get the prob. matrix for each
		# sadly, mapply() etc. crash when tried on expokit_dgexpv_wrapper
		
		# Set up empty matrix
		NA_matrix = matrix(NA, nrow=n, ncol=n)
		
		# Set up list of empty matrices
		list_of_matrices_output = replicate(NA_matrix, n=num_tvals, simplify=FALSE)
		
		for (i in 1:num_tvals)
			{
			t = tvals[i]
			list_of_matrices_output[[i]] = expokit_dgexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)

			} # end forloop
		
		return(list_of_matrices_output)
		
		} else {
		
		# If there is only 1 t value, just return 1 matrix
		#res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
		
		#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)

		t=tvals
		tmpoutmat = expokit_dgexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)		
		
		#print(tmpoutmat)
		return(tmpoutmat)
		} # End if statment
	
	}























