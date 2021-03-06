      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
      M = MOD(N,7)
      IF (M.NE.0) THEN
      DO I = 1,M
      DY(I) = DX(I)
      END DO
      IF (N.LT.7) RETURN
      END IF
      MP1 = M + 1
      DO I = MP1,N,7
      DY(I) = DX(I)
      DY(I+1) = DX(I+1)
      DY(I+2) = DX(I+2)
      DY(I+3) = DX(I+3)
      DY(I+4) = DX(I+4)
      DY(I+5) = DX(I+5)
      DY(I+6) = DX(I+6)
      END DO
      ELSE
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
      DY(IY) = DX(IX)
      IX = IX + INCX
      IY = IY + INCY
      END DO
      END IF
      RETURN
      END
