      subroutine cossin(piv,wi,ww,cos,sin)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  subroutine cossin calculates the parameters of a givens
c  transformation without square roots.
      store = piv*wi
      dd = ww+store*piv
      cos = ww/dd
      sin = store/dd
      ww = dd
      wi = cos*wi
      return
      end
