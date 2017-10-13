      subroutine deval (vmu, s, lds, nint, qraux, nobs, nnull, tol,
     *jpvt, M, ldq, nq, q, z, y, low, upp, nlaht, score, varht, info,
     *work, twk, twk2, qwork)

      ccharacter vmu
      integer ldq, n, info, jpvt(*), lds, nnull, nobs, nq, n, nint
      double precision  q(ldq,*), M(ldq,*), tol, z(*), y(*), low,
     *upp, nlaht, score, varht, twk(2,*), twk2(*), work(*),
     *qwork(ldq,*), qraux(*)
      double precision tmp, minscr, mlo, varhtwk
      integer j
      info = 0
      if( upp .lt. low )then
      mlo = low
      low = upp
      upp = mlo
      endif
      if( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u') .or. nint
     * .lt. 1 )then
      info = -3
      return
      endif
      if( 1 .gt. n .or. n .gt. ldq )then
      info = -1
      return
      endif
      j=1
23006 if(.not.(j.le.nint+1))goto 23008
      tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)
      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
                    *nq, info, work, qwork, tmp)
      call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if ( info != 0 ) then
      return
      endif
      call  dcopy (n-2, qwork(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, n-1, n-2,  work, y(n0+2),
     *dum, z(n0+2),dum, dum, dum, 01000, info)
      call dggold(vmu, M, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp,
     *mlo, score(j), varht, info, twk, work)
      if( score(j) .le. minscr .or. j .eq. 1 )then
      minscr = score(j)
      nlaht = tmp
      varhtwk = varht
      endif
23007 j=j+1
      goto 23006
23008 continue
      varht = varhtwk
      return
      end
