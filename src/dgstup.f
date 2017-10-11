      subroutine dgstup (s, M, ldm, lds, nobs, nnull, qraux, q, ldqr,
     *ldqc, nq, info, work, qwk, nla)

      integer lds, nobs, nnull, ldqr, ldqc, nq, info, ldm
      double precision s(lds,*), qraux(*), q(ldqr,ldqc,*),
     *work(*), qwk(ldqr,ldqc,*), M(ldqr,ldqc,*)
      double precision dum, nla
      integer j
      info = 0

      if ( nobs .lt. 1 .or. nobs .gt. lds .or. nobs .gt. ldqr .or. nobs
     *.gt. ldqc )then
      info = -1
      return
      endif

      do j=1,nq,1
      call dcopy(ldqr*ldqc,M,1,qwk,1)
      call dscal(ldqr*ldqc,10.d0**nla,qwk,1)
      call daxpy(ldqr*ldqc,1,q(1,1,j),1,qwk,1)
      call dqrslm(s, lds, nobs, nnull, qraux, qwk(1,1,j), ldqr, 0,
     *info, work)
      end do
      return
      end
