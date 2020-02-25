      subroutine rotate(piv,cos,sin,a,b)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  subroutine rotate applies a givens rotation to a and b.
      store = b
      b = cos*store+sin*a
      a = a-piv*store
      return
      end
