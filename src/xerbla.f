      SUBROUTINE XERBLA( SRNAME, INFO )
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
      INTRINSIC          LEN_TRIM

      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO

      STOP
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2,
     *' had ','an illegal value' )
      END
