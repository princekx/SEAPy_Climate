      subroutine disco(t,n,k2,b,mnxmy)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  subroutine disco calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  the first dimension specification of the matrix b must be the same as
c  in the calling program; h must have a dimension specification at
c  least 2*k+2.
c
      include 'param_interp.ptr'
c
      dimension t(n),b(mnxmy,k2),h(k2max2)
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j+k
          prod = 1.
          do 20 i=j,jk
            prod = prod*h(i)
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
