      subroutine dpbsl(abd,lda,n,m,b)
      implicit none
      integer lda,n,m
      double precision abd(lda,n),b(n)
      double precision ddot,t
      integer k,kb,la,lb,lm
      do 10 k = 1,n
      lm = min0(k-1,m)
      la = m + 1 - lm
      lb = k - lm
      t = ddot(lm,abd(la,k),1,b(lb),1)
      b(k) = (b(k) - t)/abd(m+1,k)
   10 continue
      do 20 kb = 1,n
      k = n + 1 - kb
      lm = min0(k-1,m)
      la = m + 1 - lm
      lb = k - lm
      b(k) = b(k)/abd(m+1,k)
      t = -b(k)
      call daxpy(lm,t,abd(la,k),1,b(lb),1)
   20 continue
      return
      end
