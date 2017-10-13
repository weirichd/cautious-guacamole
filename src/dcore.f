      subroutine dcore (vmu, s, lds, M, nobs, nnull, qraux,
     *jpvt, q, ldq, tol, y, z, job, limnla, nlaht, score,
     *varht, info, work, twk, twk2, qwork)

      character vmu
      integer ldq, nobs, nnull, job, info
      double precision q(ldq,*), tol, z(*), limnla(2), nlaht, score(*),
     *varht, twk(2,*), work(*)
      double precision dum, low, upp, dasum, mchpr
      integer n0, n, j

      info = 0

      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif

      if( nnull .lt. 1 .or. nobs .le. nnull .or. nobs .gt. ldq )then
      info = -1
      return
      endif

      n0 = nnull
      n = nobs - nnull

      if( job .eq. 0 )then
      mchpr = 1.d0

23008 if( 1.d0 + mchpr .gt. 1.d0 )then
      mchpr = mchpr / 2.d0
      goto 23008
      endif
23009 continue
      mchpr = mchpr * 2.d0
      limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
      limnla(1) = limnla(2) * mchpr
      limnla(2) = dlog10 (limnla(2))
      limnla(1) = dlog10 (limnla(1))
      endif
      low = limnla(1)
      upp = limnla(2)
      if( job .le. 0 )then
      call dgold (vmu, s, lds, qraux, nobs, nnull, tol, jpvt, M,
     *ldq, 1, q, z, y, low, upp, nlaht, score(1), varht, info,
     *work, twk, twk2, qwork)
      if( vmu .eq. 'v' )then
      score(1) = score(1) * dfloat (nobs) / dfloat (n)
      endif
      if( vmu .eq. 'm' )then
      score(1) = score(1) * dfloat (n) / dfloat (nobs)
      endif
      if( vmu .eq. 'u' )then
      score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht
      endif
      else
      call deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, nla
     *ht, score, varht, info, work, twk, qwork)
      dum = dfloat (nobs) / dfloat (n)
      j=1
23018 if(.not.(j.le.job+1))goto 23020
      if( vmu .eq. 'v' )then
      score(j) = score(j) * dum
      endif
      if( vmu .eq. 'm' )then
      score(j) = score(j) / dum
      endif
      if( vmu .eq. 'u' )then
      score(j) = score(j) / dum + 2.d0 * varht
      endif
23019 j=j+1
      goto 23018
23020 continue
      endif
      return
      end
