C Output from Public domain Ratfor, version 1.0
      subroutine  dgold (vmu, s, lds, qraux, nobs, nnull, tol, jpvt, M,
     *ldq, nq, q, z, y, low, upp, nlaht, score,
     *varht,info, work, twk, twk2, qwork)
C Varible declares
      character         vmu
      integer           ldq, info, jpvt(*), lds,
     *nnull, nobs, n, nq
      double precision  q(ldq,*), M(ldm,*), tol, z(*), y(*),
     *low, upp, nlaht, score, varht, twk(2,*), twk2(*), work(*),
     *qwork(ldq,*), qraux(*)
      double precision ratio, mlo, mup, tmpl, tmpu

C Begin the code
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      n0 = nnull
      n = nobs - nnull
      info = 0
      if( upp .lt. low )then
      mlo = low
      low = upp
      upp = mlo
      endif
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif
      if( n .lt. 1 .or. n .gt. ldq )then
      info = -1
      return
      endif
      mlo = upp - ratio * (upp - low)

      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, 10.d0**mlo)
      call dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if ( info .ne. 0 ) then
      return
      endif
      call dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)

      call dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2, work,
     *y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
      if( info .ne. 0 )then
      return
      endif
      call dggold(vmu, M, qwork(n0+1,n0+1), ldq, n, z(n0+1),
     *low, upp, mlo, tmpl, varht, info, twk, work)


      mup = low + ratio * (upp - low)
      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, 10.d0**mup)
      call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)

      if ( info .ne. 0 )then
      return
      endif

      call dcopy (n-2, qqwork(n0+2,n0+1), ldq+1, work, 1)
      call dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2, work,
     *y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
      if( info .ne. 0 )then
      return
      endif
      call dggold (vmu, M, qwork(n0+1,n0+1), ldq, n, z(n0+1),
     *low, upp, mup, tmpu, varht, info, twk, work)

23010 continue
      if( mup - mlo .lt. 1.d-7 )then
      goto 23012
      endif
      if( tmpl .lt. tmpu )then
      upp = mup
      mup = mlo
      tmpu = tmpl
      mlo = upp - ratio * (upp - low)
      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, 10.d0**mlo)

      call dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if ( info .ne. 0 )then
      return
      endif
      call dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)

      call dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2, work,
     *y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
      if ( info .ne. 0 )then
      return
      endif
      call dggold(vmu, M, qwork(n0+1,n0+1), ldq, n, z(n0+1),
     *low, upp, mlo, tmpl, varht, info, twk, work)
      else
      low = mlo
      mlo = mup
      tmpl = tmpu
      mup = low + ratio * (upp - low)
      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, 10.d0**mup)

      call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if ( info .ne. 0 )then
      return
      endif
      call dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
      call dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2, work,
     *y(n0+2), dum, z(n0+2), dum, dum, dum, 01000, info)
      if( info .ne. 0 )then
      return
      endif
      call dggold (vmu, M, qwork(n0+1,n0+1), ldq, n, z(n0+1),
     *low, upp, mup, tmpu, varht, info, twk, work)

      endif
23011 goto 23010
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, 10.d0**nlaht)

      call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if ( info .ne. 0 )then
      return
      endif

      call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)

      call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2, work,
     *y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
      if( info .ne. 0 )then
      return
      endif
      
      call dggold(vmu, M, qwork(n0+1,n0+1), ldq, n, z(n0+1), low, upp,
     *nlaht, score, varht, info, twk, work)
      return
      end
