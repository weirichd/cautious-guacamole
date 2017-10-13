      subroutine  dsidr (vmu, s, lds, nobs, nnull, y, z, M, q,
     *ldq, tol, job, limnla, nlaht, score, varht, c, d, qraux, jpvt,
     *wk, twk, twk2, qwk, info)

      character vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info
      double precision  s(lds,*), y(*), z(*), M(ldq,*), q(ldq,*),
     *tol, limnla(2), nlaht, score(*), varht, c(*), d(*), qraux(*),
     *wk(*), qwk(ldq,*), twk(*), twk2(*)

      info = 0
      if( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds .or. nobs
     * .gt. ldq )then
      info = -1
      return
      endif

      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
      info = -3
      return
      endif

      call  dstup (s, lds, nobs, nnull, qraux, jpvt, y, info, wk)
      if( info .ne. 0 )then
      return
      endif

      call  dcopy(nobs, y, 1, z, 1)
      call  dcore(vmu, s, lds, M, nobs, nnull, qraux, jpvt, q, ldq,
     *tol, y, z, job, limnla, nlaht, score, varht, info, wk, twk,
     *twk2, qwk)
      if( info .ne. 0 )then
      return
      endif
    
      call dcoef (s, lds, nobs, nnull, M, qraux, jpvt, y, z, q, ldq,
     *nlaht, c, d, info, work, qwk, twk)

      return
      end
