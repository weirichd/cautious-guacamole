
#:::::::::::
#   dcore
#:::::::::::

subroutine  dcore (vmu, s, lds, M, ldmr, ldmc, nobs, nnull, qraux,_
                   jpvt, q, ldqr, ldqc, tol, y, z, 
                   job, limnla, nlaht, score, varht, info, twk,_
                   twk2, work, qwork)

#  Purpose:  To evaluate the GCV/GML score function at various trial values
#      of n*lambda using the tridiagonalization GCV/GML algorithm.  Perform
#      either golden section search or regular grid search for minimizing
#      n*lambda.

character*1       vmu
integer           ldqr, ldqc, ldmr, ldmc, lds, nobs, nnull, job, info,_
                  jpvt(*)
double precision  s(lds,*), q(ldqr,ldqc,*), M(ldmr,ldmc,*), tol, y(*), z(*), ,_
                  limnla(2), nlaht, score(*), varht, twk(2,*), twk2(*), work(*),_
                  qwork(ldmr,ldmc,*), qraux(*)


#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      s,qraux,jpvt
#                 QR decomposition of S=FR.
#      q          the reproducing kernels, of size (ldqr,ldqc,nq).
#      ldqr       leading dimension for rows of q.
#      ldqc       leading dimension for columns of q.
#      nq         number of Q's.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      tol        tolerance of truncation.
#      z          F^T y.
#      y          F^T y
#      job         0:  searching interval for nlaht chosen automatically.
#                 -1:  searching interval for nlaht provided by limnla.
#                 >0:  search regular grid points on [limnla(1),limnla(2)]:
#                        #(grids) = job + 1.
#      limnla     searching interval in log10 scale, see job.
#      varht      known variance if vmu=='u'.

#  On exit:
#      qwk      tridiagonal form of (Q + nlaht*M) in diagonal and superdiagonal of the
#                 corner, Householder factors in strict lower triangle of
#                 the corner.
#      y          F^{T} y.
#      z          diag(I, U^{T}) F^{T} y.
#      limnla     see limnla of entry.
#      nlaht      the estimated log10(n*lambda).
#      score      job <= 0 :  the GCV/GML/URE score at nlaht.
#                 job  > 0 :  the GCV/GML/URE score at the regular grid points.
#      varht      variance estimate.
#      info        0 :  normal termination.
#                 -1 :  dimension error.
#                 -2 :  F_{2}^{T}QF_{2} is not non-negative definite.
#                 -3 :  vmu is none of 'v', 'm', or 'u'.

#  Work arrays:
#      twk        of size at least (2,nobs-nnull).
#      work       of size at least (nobs-nnull).

#  Routines called directly:
#      Fortran -- dfloat, dlog10
#      Blas    -- dasum, dcopy
#      Linpack -- dqrsl
#      Rkpack  -- deval, dgold, dsytr

#  Written:  Chong Gu, Statistics, Purdue, latest version 3/24/92.

double precision  dum, low, upp, dasum, mchpr
integer           n0, n, j, nq, ldqr, ldqc

info = 0

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   check dimension
if ( nnull < 1 | nobs <= nnull | nobs > ldq ) {
    info = -1
    return
}

#   set searching range
if ( job == 0 ) {
    mchpr = 1.d0
    while ( 1.d0 + mchpr > 1.d0 )  mchpr = mchpr / 2.d0
    mchpr = mchpr * 2.d0
    limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
    limnla(1) = limnla(2) * mchpr
    limnla(2) = dlog10 (limnla(2))
    limnla(1) = dlog10 (limnla(1))
}

low = limnla(1)
upp = limnla(2)
if ( job <= 0 ) {
    #   compute score and estimate nlaht thru golden-section search
    call dgold (vmu, s, lds, qraux, nobs, nnull, tol, jpvt, M, ldmr, ldmc,
                1, q, ldqr, ldqc, n, z, y, low, upp, nlaht,_
                score(1), varht, info, twk, twk2, work, qwork)
    if ( vmu == 'v' )  score(1) = score(1) * dfloat (nobs) / dfloat (n)
    if ( vmu == 'm' )  score(1) = score(1) * dfloat (n) / dfloat (nobs)
    if ( vmu == 'u' )  score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht
}
else {
    #   regular grid evaluation
    call  deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, nlaht,_
                 score, varht, info, twk, work)
    dum = dfloat (nobs) / dfloat (n)
    for (j=1;j<=job+1;j=j+1) {
        if ( vmu == 'v' )  score(j) = score(j) * dum
        if ( vmu == 'm' )  score(j) = score(j) / dum
        if ( vmu == 'u' )  score(j) = score(j) / dum + 2.d0 * varht
    }
}

return
end

#...............................................................................

