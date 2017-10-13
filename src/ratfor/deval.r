
#:::::::::::
#   deval
#:::::::::::

      subroutine  deval (vmu, s, lds, nint, qraux, nobs, nnull, tol, jpvt, M,
     *ldq, nq, q, z, y, low, upp, nlaht, score,
     *varht,info, work, twk, twk2, qwork)




      character       vmu
      integer           ldq, n, info, jpvt(*), lds,
     *nnull, nobs, nq, n, nint
      double precision  q(ldq,*), M(ldq,*), tol, z(*), y(*), low,
     *upp, nlaht, score, varht, twk(2,*), twk2(*), work(*),
     *qwork(ldq,*), qraux(*)

#  Purpose:  To evaluate GCV/GML function based on tridiagonal form and to
#      search minimum on an interval by equally spaced (in log10 scale) grid
#      search.

character*1       vmu
integer           ldq, n, nint, info
double precision  q(ldq,*), z(*), low, upp, nlaht, score(*), varht,_
                  twk(2,*), work(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      q          tidiagonal matrix in diagonal and super diagonal.
#      ldq        leading dimension of Q.
#      n          size of the matrix.
#      z          U^{T} F_{2}^{T} y.
#      nint       number of intervals (number of grids minus 1).
#      low        lower limit of log10(n*lambda).
#      upp        upper limit of log10(n*lambda).
#      varht      known variance if vmu=='u'.

#  On exit:
#      nlaht      the estimated log10(n*lambda).
#      score      the GCV/GML/URE score vector on grid points.
#      varht      the variance estimate at the estimated n*lambda.
#      info        0: normal termination.
#                 -1: dimension error.
#                 -2: tridiagonal form is not non-negative definite.
#                 -3: vmu or nint is out of scope.

#  Work arrays:
#      twk        array of length at least (2,n).
#      work       array of length at least (n).

#  Routines called directly:
#      Fortran -- dfloat
#      Blas    -- daxpy, dcopy
#      Rkpack  -- dtrev
#      Other   -- dset

#  Written:  Chong Gu, Statistics, Purdue, 12/29/91 latest version.

double precision  tmp, minscr, mlo, varhtwk
integer           j

info = 0

#   interchange boundaries if necessary
if ( upp < low ) {
    mlo = low
    low = upp
    upp = mlo
}

#   check job requests
if ( (vmu != 'v' & vmu != 'm' & vmu != 'u') | nint < 1 ) {
    info = -3
    return
}

#   check dimension
if ( 1 > n | n > ldq ) {
    info = -1
    return
}

#   evaluation
for (j=1;j<=nint+1;j=j+1) {
    tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)

    call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
                  *nq, info, work, qwork, mlo)
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

    ###-----------
        if ( score(j) <= minscr | j == 1 ) {
        minscr = score(j)
        nlaht = tmp
        varhtwk = varht
    }
}
varht = varhtwk

return
end
