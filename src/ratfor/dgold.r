
#:::::::::::
#   dgold
#:::::::::::

subroutine dgold (vmu, s, lds, qraux, nobs, nnull, tol, jpvt, M,
       *ldq, nq, q, z, y, low, upp, nlaht, score,
       *varht,info, work, twk, twk2, qwork)

#  Purpose:  To evaluate GCV/GML function based on tridiagonal form and to
#      search minimum on an interval by golden section search.

character*1       vmu
integer           ldq, n, info, jpvt(*), lds,_
                  nnull, nobs, nq, n
double precision  q(ldq,*), M(ldq,*), tol, z(*), y(*), low,_
                  upp, nlaht, score, varht, twk(2,*), twk2(*), work(*),_
                  qwork(ldq,*), qraux(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      q          the reproducing kernels, of size (ldqr,ldqc,nq).
#      ldq
#      nq         number of Q's.
#      n          size of the matrix.
#      z          U^{T} F_{2}^{T} y.
#      low        lower limit of log10(n*lambda).
#      upp        upper limit of log10(n*lambda).
#      varht      known variance if vmu=='u'.

#  On exit:
#      nlaht      the estimated log(n*lambda).
#      score      the GCV/GML/URE score at the estimated lambda.
#      varht      the variance estimate at the estimated lambda.
#      info        0: normal termination.
#                 -1: dimension error.
#                 -2: tridiagonal form is not non-negative definite.
#                 -3: vmu is none of 'v', 'm', or 'u'.

#  Work arrays:
#      twk        of size at least (2,n).
#      work       of size at least (n).

#  Routines called directly:
#      Fortran -- dsqrt
#      Blas    -- daxpy, dcopy
#      Rkpack  -- dtrev
#      Other   -- dset


double precision  ratio, mlo, mup, tmpl, tmpu

ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
#   set working parameters
n0 = nnull
n = nobs - nnull

info = 0
#   interchange the boundaries if necessary
if ( upp < low ) {
    mlo = low
    low = upp
    upp = mlo
}

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   check dimension
if ( n < 1 | n > ldq ) {
    info = -1
    return
}

###############################################################################
#   initialize golden section search for scrht
mlo = upp - ratio * (upp - low)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## in:     s,qraux,jpvt
#                 QR decomposition of S=FR
#          q          the reproducing kernels, of size (ldqr,ldqc,nq).
## out:    qwork  F^{T}(Q + mlo*M)F's.
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Calculate Qwork = Q + lambda M +++++++++++++++++++++++++++++++++++
## Decompose Qwork := F2^t Qwork F2
call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
              *nq, info, work, qwork, mlo)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   tridiagonalization
##          U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda) = T
##   on exit:
##          qwork(n0+1,n0+1) -
##                   diagonal:  diagonal elements of tridiag. transf.
##                   upper triangle:  off-diagonal of tridiag. transf.
##                   lower triangle:  overwritten by Householder factors.
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
if ( info != 0 )  return
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   z(lambda) : = U(lambda)^{T} z_{2}
## copy lower triangle of U T U^T into work
call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
## z(n0+2) := U^T F_2^T y(n0+2)
call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, mlo, tmpl,
*varht, info, twk, work)

###------------------------------------------------------------------
mup = low + ratio * (upp - low)
## Calculate Qwork = Q + lambda M +++++++++++++++++++++++++++++++++++
call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
              *nq, info, work, qwork, mlo)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   tridiagonalization
##          U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda) = T
##   on exit:
##          qwork(n0+1,n0+1) -
##                   diagonal:  diagonal elements of tridiag. transf.
##                   upper triangle:  off-diagonal of tridiag. transf.
##                   lower triangle:  overwritten by Householder factors.
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
if ( info != 0 )  return
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   z(lambda) : = U(lambda)^{T} z_{2}
## copy lower triangle of U T U^T into work
call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
## z(n0+2) := U^T F_2^T y(n0+2)
call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, mup, tmpu,
*varht, info, twk, work)

if ( info != 0 ) {
    info = -2
    return
}

#   golden section search for estimate of lambda
repeat {
    if ( mup - mlo < 1.d-7 )  break
    if ( tmpl < tmpu ) {
        upp = mup
        mup = mlo
        tmpu = tmpl
        mlo = upp - ratio * (upp - low)
        ## Calculate Qwork = Q + lambda M +++++++++++++++++++++++++++++++++++
        call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
                      *nq, info, work, qwork, mlo)
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ##   tridiagonalization
        ##          U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda) = T
        ##   on exit:
        ##          qwork(n0+1,n0+1) -
        ##                   diagonal:  diagonal elements of tridiag. transf.
        ##                   upper triangle:  off-diagonal of tridiag. transf.
        ##                   lower triangle:  overwritten by Householder factors.
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
        if ( info != 0 )  return
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ##   z(lambda) : = U(lambda)^{T} z_{2}
        ## copy lower triangle of U T U^T into work
        call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
        ## z(n0+2) := U^T F_2^T y(n0+2)
        call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
        call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, mlo, tmpl,
        *varht, info, twk, work)
    }
    else {
        low = mlo
        mlo = mup
        tmpl = tmpu
        mup = low + ratio * (upp - low)

        ## Calculate Qwork = Q + lambda M +++++++++++++++++++++++++++++++++++
        call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
                      *nq, info, work, qwork, mlo)
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ##   tridiagonalization
        ##          U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda) = T
        ##   on exit:
        ##          qwork(n0+1,n0+1) -
        ##                   diagonal:  diagonal elements of tridiag. transf.
        ##                   upper triangle:  off-diagonal of tridiag. transf.
        ##                   lower triangle:  overwritten by Householder factors.
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
        if ( info != 0 )  return
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ##   z(lambda) : = U(lambda)^{T} z_{2}
        ## copy lower triangle of U T U^T into work
        call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
        ## z(n0+2) := U^T F_2^T y(n0+2)
        call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
        call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, mup, tmpu,
        *varht, info, twk, work)

    }
}

#   compute the return value
nlaht = ( mup + mlo ) / 2.d0
## Calculate Qwork = Q + lambda M +++++++++++++++++++++++++++++++++++
call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
              *nq, info, work, qwork, mlo)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   tridiagonalization
##          U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda) = T
##   on exit:
##          qwork(n0+1,n0+1) -
##                   diagonal:  diagonal elements of tridiag. transf.
##                   upper triangle:  off-diagonal of tridiag. transf.
##                   lower triangle:  overwritten by Householder factors.
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
if ( info != 0 )  return
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##   z(lambda) : = U(lambda)^{T} z_{2}
## copy lower triangle of U T U^T into work
call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
## z(n0+2) := U^T F_2^T y(n0+2)
call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, nlaht, score,
*varht, info, twk, work)

return
end

#...............................................................................
