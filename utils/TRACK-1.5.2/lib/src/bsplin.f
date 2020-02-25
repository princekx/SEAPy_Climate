      subroutine bsplin(t,n,k,x,l,h,hh)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  subroutine bsplin evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  the dimension specifications of the following arrays must be
c  at least h(k+1),hh(k).
      dimension t(n),h(k+1),hh(k)
      h(1) = 1.
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
