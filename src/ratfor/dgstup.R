      subroutine dgstup (s, M, ldm, lds, nobs, nnull, qraux, q, ldqr,
                   *ldqc, nq, info, work, qwk, nla)

      integer lds, nobs, nnull, ldqr, ldqc, nq, info, ldm
      double precision s(lds,*), qraux(*), q(ldqr,ldqc,*),
      *work(*), q(ldqr,ldqc,*), M(ldqr,ldqc,*)
      double precision dum, nla
      integer j
      info = 0


#   check dimension
if ( nobs < 1 | nobs > lds | nobs > ldqr | nobs > ldqc ) {
      info = -1
      return
}

for (j=1;j<=nq;j=j+1) {
      # set qwk = M
      dcopy(ldmr*ldmc,M,1,qwk,1)
      # set qwk = nla*qwk = nla*M
      dscal(ldmr*ldmc,nla,qwk,1) 
      # qwk = 1*q + qwk
      daxpy(ldq*ldq,1,q(1,1,j),1,qwk,1)
      dqrslm (s, lds, nobs, nnull, qraux, qwk(1,1,j), ldmr, 0,_
              info, work)
}
return
