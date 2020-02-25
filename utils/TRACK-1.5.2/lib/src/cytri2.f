c
      subroutine cytri2(a,n,b,c)
      implicit double precision (A-H, O-Z)
c subroutine cytri2 solves a linear n x n system
c         a * c = b
c where matrix a is a cyclic tridiagonal matrix, decomposed
c using subroutine cytri1.
      dimension a(6,n),b(n),c(n)
      c(1) = b(1)*a(4,1)
      sum = c(1)*a(5,1)
      n1 = n-1
      do 10 i=2,n1
         c(i) = (b(i)-a(1,i)*c(i-1))*a(4,i)
         sum = sum+c(i)*a(5,i)
  10  continue
      cc = (b(n)-sum)*a(4,n)
      c(n) = cc
      c(n1) = c(n1)-cc*a(6,n1)
      j = n1
      do 20 i=3,n
         j1 = j-1
         c(j1) = c(j1)-c(j)*a(3,j1)*a(4,j1)-cc*a(6,j1)
         j = j1
  20  continue
      return
      end
