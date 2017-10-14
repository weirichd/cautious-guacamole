      subroutine dpbfa(abd,lda,n,m,info)
      implicit none
      integer lda,n,m,info
      double precision abd(lda,n)
      double precision ddot,t
      double precision s
      integer ik,j,jk,k,mu
      do 30 j = 1, n
      info = j
      s = 0.0e0
      ik = m + 1
      jk = max0(j-m,1)
      mu = max0(m+2-j,1)
      if (m .lt. mu) go to 20
      do 10 k = mu, m

      t = abd(k,j) - ddot(k-mu,abd(ik,jk),1,abd(mu,j),1)
      t = t/abd(m+1,jk)
      abd(k,j) = t
      s = s + t*t
      ik = ik - 1
      jk = jk + 1
   10 continue
   20 continue
      s = abd(m+1,j) - s
      if (s .le. 0.0e0) go to 40
      abd(m+1,j) = sqrt(s)
   30 continue
      info = 0
   40 continue
      return
      end
