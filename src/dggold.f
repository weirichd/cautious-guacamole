      subroutine dggold (vmu, M, q, ldq, n, z, low, upp, nlaht,
     *score, varht, info, twk, work)

      character*1       vmu
      integer           ldq, n, info
      double precision  q(ldq,*), M(ldq,*), z(*), low, upp,
     *nlaht, score, varht, twk(2,*), work(*)

      call  dset (n, 0.d0, twk(2,1), 2)
      call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)

      call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)

      twk(1,1) = 10.d0**nlaht

      call  dtrev (vmu, twk, 2, M, ldq, n, z, score, varht, info,
     *work, twk2)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      return
      end
