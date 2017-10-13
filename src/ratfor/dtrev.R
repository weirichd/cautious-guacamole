
#:::::::::::
#   dtrev
#:::::::::::

subroutine  dtrev (vmu, t, ldt, M, ldq, n, z, score, varht, info,
                   work, twk)

#  Acronym:  Double-precision TRidiagonal EValuation.

#  Purpose:  To compute the GCV/GML function and the related variance
#      estimate from the tridiagonal matrix `t' and data vector `z'.

#  References:  1. Gu, Bates, Chen, and Wahba(1988), TR#823, Stat, UW-M.
#               2. Dongarra et al. (1979) LINPACK User's Guide. (Chap. 4)

character*1       vmu
integer           n, info, ldt, ldq
double precision  t(ldt,*), M(ldq,*), z(*), score, varht, work(*),
                  twk(*)

#  On entry:
#      vmu        'v':  GCV.
#                 'm':  GML.
#                 'u':  unbiased risk estimate.
#      t          the positive definite tridiagonal matrix  T,
#                 stored in packed form:
#                     t(1,2:n):  off-diagonal
#                     t(2,1:n):  diagonal.
#      ldt        leading dimension of t.
#      n          the dimension of the matrix.
#      z          the appropriately transformed data vector.
#      varht      known variance if vmu=='u'.

#  On exit:
#      score      the GCV/GML/URE score.
#      varht      \hat\sigma^{2}.
#      info         -3:  vmu is none of 'v', 'm', or 'u'.
#                 > -3:  as from LINPACK's `dpbfa'.

#  Work array:
#      work       of size at least (n).

#  Routines called directly:
#      Fortran -- dexp, dfloat, dlog
#      Blas    -- dasum, dcopy, ddot, dscal
#      Linpack -- dpbfa, dpbsl

double precision  nume, deno, tmp, alph, la, dasum, ddot
integer           j

info = 0

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

la = t(1,1)

#   standardize the matrix for numerical stability
alph = dfloat (n) / dasum (n, t(2,1), ldt)
call  dscal (n, alph, t(2,1), ldt)
call  dscal (n-1, alph, t(1,2), ldt)

#   decomposition
## of T(lambda)
## T(lambda) = U(lambda)^{T} [F^{T}( Q + lambda*M )F] U(lambda)
call  dpbfa (t, ldt, n, 1, info)
if ( info != 0 )  return

## set work := U^T F_2^T y
call  dcopy (n, z, 1, work, 1)

## solve T(lambda) x = z
## T^{-1} = U^T [F^{T}( Q + lambda*M )F]^{-1} U
## ## sets work :=  U^T[T(lambda)]^{-1} U U^T F2^T Y
##              := [F^{T}( Q + lambda*M )F]^{-1} F2^T Y
call  dpbsl (t, ldt, n, 1, work)

## compute work :=F2 [F2^T (Q + lambda M ) F2]^{-1}F2^T Y
call  dqrsl (s, lds, nobs, nnull, qraux, work, work, dum, work, dum, dum, 10000,_
             info)

## premultiply work by M ##--------------------------------------
## dgemm(TRANSA, TRANSB, nrow(op(A)), ncol(op(B)), ncol(op(A)),
##       ALPHA, A, LDA, B, LDB, BETA, C, LDC)
## C := alpha*op( A )*op( B ) + beta*C,

call dcopy (n, work, 1, twk, 1)
call dgemm ( "n", "n", ldm, 1, ldm, 1, M, ldm,
                   work, ldt, 0, twk, ldt)

#   GCV computation
if ( vmu == 'v' ) {
    tmp = 1.d0 / t(2,n) / t(2,n)
    deno = tmp
    for (j=n-1;j>0;j=j-1) {
        tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
        deno = deno + tmp
    }
    nume = ddot (n, work, 1, work, 1) / dfloat (n)
    deno = deno / dfloat (n)
    varht = alph * la * nume / deno
    score = nume / deno / deno
}

#   GML computation
if ( vmu == 'm' ) {
    deno = dlog (t(2,n))
    for (j=n-1;j>0;j=j-1)  deno = deno + dlog (t(2,j))
    nume = ddot (n, z, 1, work, 1) / dfloat (n)
    varht = alph * la * nume
    score = nume * dexp (2.d0 * deno / dfloat (n))
}

#   unbiased risk computation
if ( vmu == 'u' ) {
    nume = ddot (n, work, 1, work, 1) / dfloat (n)
    tmp = 1.d0 / t(2,n) / t(2,n)
    deno = tmp
    for (j=n-1;j>0;j=j-1) {
        tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
        deno = deno + tmp
    }
    deno = deno / dfloat (n)
    score = alph * alph * la * la * nume - 2.d0 * varht * alph * la * deno
}

return
end

#...............................................................................
