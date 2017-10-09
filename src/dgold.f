C Output from Public domain Ratfor, version 1.0
      subroutine  dgold (vmu, s, lds, qraux, nobs, nnull, tol, jpvt, M,
     *ldmr, ldmc, nq, q, ldqr, ldqc, z, y, low, upp, nlaht, score,
     *varht,info, twk, twk2, work, qwork)
C Varible declares
      character         vmu
      integer           ldmr, ldmc, ldqr, ldqc, info, jpvt(*), lds,
     *nnull, nobs, ncol, nq
      double precision  q(ldqr,ldqc,*), M(ldmr,ldmc,*), tol, z(*), y(*),
     *low, upp, nlaht, score, varht, twk(2,*), twk2(*), work(*),
     *qwork(ldmr,ldmc,*), qraux(*)
      double precision ratio, mlo, mup, tmpl, tmpu
C Begin the code
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      n0 = nnull
      ncol = nobs - nnull
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
      if( ncol .lt. 1 .or. ncol .gt. ldq )then
      info = -1
      return
      endif
      mlo = upp - ratio * (upp - low)
      do j=1,nq,1
      call dcopy(ldmr*ldmc,M,1,qwork,1)
      call dscal(ldmr*ldmc,mlo,qwork,1)
      call daxpy(ldq*ldq,1,q(1,1,j),1,qwork,1)
      call dqrslm (s, lds, nobs, nnull, qraux, qwork(1,1,j), ldmr, 0,
     *info, work)
      end do

      call  dsytr (qwork(n0+1,n0+1), ldmr, ncol, tol, info, work)
      if ( info .ne. 0 )  return
      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2,  work, y(n0+2),
     *dum, z(n0+2),dum, dum, dum, 01000, info)
      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call  dtrev (vmu, twk, 2, M, nobs, ncol, z, mlo, varht, info,
     *work, twk2)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      mup = low + ratio * (upp - low)
      do j=1,nq,1
      call dcopy(ldmr*ldmc,M,1,qwork,1)
      call dscal(ldmr*ldmc,mup,qwork,1)
      call daxpy(ldq*ldq,1,q(1,1,j),1,qwork,1)
      call dqrslm (s, lds, nobs, nnull, qraux, qwork(1,1,j), ldmr, 0,
     *info, work)
      end do

      call  dsytr (qwork(n0+1,n0+1), ldmr, ncol, tol, info, work)
      if ( info .ne. 0 )  return
      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2,  work, y(n0+2),
     *dum, z(n0+2),dum, dum, dum, 01000, info)
      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call  dtrev (vmu, twk, 2, M, nobs, ncol, z, mup, varht, info,
     *work, twk2)
      if( info .ne. 0 )then
      info = -2
      return
      endif
23010 continue
      if( mup - mlo .lt. 1.d-7 )then
      goto 23012
      endif
      if( tmpl .lt. tmpu )then
      upp = mup
      mup = mlo
      tmpu = tmpl
      mlo = upp - ratio * (upp - low)
      do j=1,nq,1
      call dcopy(ldmr*ldmc,M,1,qwork,1)
      call dscal(ldmr*ldmc,mlo,qwork,1)
      call daxpy(ldq*ldq,1,q(1,1,j),1,qwork,1)
      call dqrslm (s, lds, nobs, nnull, qraux, qwork(1,1,j), ldmr, 0,
     *info, work)
      end do

      call  dsytr (qwork(n0+1,n0+1), ldmr, ncol, tol, info, work)
      if ( info .ne. 0 )  return
      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2,  work,
     *y(n0+2), dum, z(n0+2),dum, dum, dum, 01000, info)
      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, ncol, z, tmpl, varht, info, work)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      else
      low = mlo
      mlo = mup
      tmpl = tmpu
      mup = low + ratio * (upp - low)
      do j=1,nq,1
      call dcopy(ldmr*ldmc,M,1,qwork,1)
      call dscal(ldmr*ldmc,mup,qwork,1)
      call daxpy(ldq*ldq,1,q(1,1,j),1,qwork,1)
      call dqrslm (s, lds, nobs, nnull, qraux, qwork(1,1,j), ldmr, 0,
     *info, work)
      end do

      call  dsytr (qwork(n0+1,n0+1), ldmr, ncol, tol, info, work)
      if ( info .ne. 0 )  return
      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2,  work, y(n0+2),
     *dum, z(n0+2),dum, dum, dum, 01000, info)
      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call  dtrev (vmu, twk, 2, ncol, z, tmpl, varht, info, work, twk2)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      endif
23011 goto 23010
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      do j=1,nq,1
      call dcopy(ldmr*ldmc,M,1,qwork,1)
      call dscal(ldmr*ldmc,nlaht,qwork,1)
      call daxpy(ldq*ldq,1,q(1,1,j),1,qwork,1)
      call dqrslm (s, lds, nobs, nnull, qraux, qwork(1,1,j), ldmr, 0,
     *info, work)
      end do

      call  dsytr (qwork(n0+1,n0+1), ldmr, ncol, tol, info, work)
      if ( info .ne. 0 )  return
      call  dcopy (ncol-2, q(n0+2,n0+1), ldq+1, work, 1)
      call  dqrsl (qwork(n0+2,n0+1), ldq, ncol-1, ncol-2,  work, y(n0+2),
     *dum, z(n0+2),dum, dum, dum, 01000, info)
      call  dset (ncol, 0.d0, twk(2,1), 2)
      call  daxpy (ncol, 1.d0, q, ldq+1, twk(2,1), 2)
      call  dcopy (ncol-1, qwork(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call  dtrev (vmu, twk, 2, M, nobs, ncol, z, nlaht, varht, info,
     *work, twk2)
      if( info .ne. 0 )then
      info = -2
      return
      endif
      return
      end
