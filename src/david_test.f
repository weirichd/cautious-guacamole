      program test
      
      integer i
      double precision lambda, Q(10), M(10), Qwork(10)
      
      lambda = 6d0

      do i = 1, 10
           M(i) = i * 5d-1
           Q(i) = 10
      enddo

17    format(10f7.2)

      write(*,*) "lambda ="
      write(*,17) lambda

      write(*,*) "M = ["
      write(*,17) (M(i), i = 1, 10)
      write(*,*) "]"

      write(*,*) "Q = ["
      write(*,17) (Q(i), i = 1, 10)
      write(*,*) "]"

      write(*,*)
      write(*,*) "Calling dcopy..."
      call dcopy(10, M, 1, Qwork, 1)

      write(*,*) "Qwork = ["
      write(*,17) (Qwork(i), i = 1, 10)
      write(*,*) "]"

      write(*,*)
      write(*,*) "Calling dscal..."
      call dscal(10, lambda, Qwork, 1)

      write(*,*) "Qwork = ["
      write(*,17) (Qwork(i), i = 1, 10)
      write(*,*) "]"

C This doesn't work! Here we pass the integer 1 to DA in daxpy
      write(*,*)
      write(*,*) "Calling daxpy the first time... this does NOT work..."
      call daxpy(10, 1, Q, 1, Qwork, 1)

      write(*,*) "Qwork = ["
      write(*,17) (Qwork(i), i = 1, 10)
      write(*,*) "]"

C You need to pass a double literal to DA instead. When we use 1.0d
C instead it works just as required
      write(*,*)
      write(*,*) "Calling daxpy the second time... this does work..."
      call daxpy(10, 1.0d0, Q, 1, Qwork, 1)

      write(*,*) "Qwork = ["
      write(*,17) (Qwork(i), i = 1, 10)
      write(*,*) "]"

      end program
