C Output from Public domain Ratfor, version 1.0
      subroutine dcoef (s, lds, nobs, nnull, M, qraux, jpvt, y, z, q, ldq,
     *nlaht, c, d, info, work, qwk, twk)

      integer lds, nobs, nnull, jpvt(*), ldq, info
      double precision s(lds,*), qraux(*), y(*), q(ldq,*), nlaht, c(*),
     *d(*), twk(*), M(ldq,*), work(*), qwk(ldq,*)
      double precision dum, ddot
      integer n, n0
      info = 0

      if( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds .or. nobs
     * .gt. ldq )then
      info = -1
      return
      endif

      n0 = nnull
      n = nobs - nnull

      call dgstup ( s, M, lds, nobs, nnull, qraux, q, ldq, nobs,
     *nq, info, work, qwork, nlaht)

      call  dsytr (qwork(n0+1,n0+1), ldq, n, tol, info, work)
      if( info .ne. 0 )then
      return
      endif

      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2, work,
     *y(n0+2),dum, z(n0+2),dum, dum, dum, 01000, info)

      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**nlaht

      call  dpbfa (twk, 2, n, 1, info)
      if( info .ne. 0 )then
      info = -2
      return
      endif

      call  dpbsl (twk, 2, n, 1, z(n0+1))
      call  dcopy (n-2, q(n0+2,n0+1), ldq+1, twk, 1)
      call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, twk, z(n0+2), z(n0+2),
     *dum, dum, dum, dum, 10000, info)

      call  dset (n0, 0.d0, c, 1)
      call  dcopy (n, z(n0+1), 1, c(n0+1), 1)
      call  dqrsl (s, lds, nobs, nnull, qraux, c, c, dum, dum, dum,
     *dum, 10000,info)

      j=1
23004 if(.not.(j.le.n0))goto 23006
      d(j) = z(j) - ddot (n, z(n0+1), 1, q(n0+1,j), 1)
23005 j=j+1
      goto 23004
23006 continue
      call dtrsl (s, lds, n0, d, 01, info)
      call dprmut (d, n0, jpvt, 1)
      return
      end
