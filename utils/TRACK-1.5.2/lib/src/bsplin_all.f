      subroutine bsplin_all(t,n,k,x,l,h)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  subroutine bsplin_all evaluates all the non-zero b-splines of
c  degree k > 0 at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  the dimension specifications of the following arrays must be
c  at least h(k+1),hh(k).
      dimension t(n),h(k+1,k+1)
      h(1,1) = 1.
      do 20 j=1,k
        h(1,j+1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = h(i,j)/(t(li)-t(lj))
          h(i,j+1) = h(i,j+1)+f*(t(li)-x)
          h(i+1,j+1) = f*(x-t(lj))
  20  continue
      return
      end
