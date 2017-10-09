C Pre-lambda loop setup
      subroutine dstup (s, lds, nobs, nnull, qraux, jpvt, y, info, work)

      integer           lds, nobs, nnull, jpvt(*), info
      double precision  s(lds,*), y(*), qraux(*), work(*)
      double precision dum
      integer j
      info = 0

      if( nobs .lt. 1 .or. nobs .gt. lds )then
      info = -1
      return
      endif
      j=1
03002 if(.not.(j.le.nnull))goto 03004
      jpvt(j) = 0
03003 j=j+1
      goto 03002
03004 continue
      call dqrdc (s, lds, nobs, nnull, qraux, jpvt, work, 1)
      call dqrsl (s, lds, nobs, nnull, qraux, y, dum, y, work, dum, dum,
     * 01100, info)
      return
      end
