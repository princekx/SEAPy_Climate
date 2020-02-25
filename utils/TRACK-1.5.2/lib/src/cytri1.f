c
      subroutine cytri1(a,n)
      implicit double precision (A-H, O-Z)
c (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
c elements stored as follows
c
c    ! a(2,1) a(3,1)                                    a(1,1)  !
c    ! a(1,2) a(2,2) a(3,2)                                     !
c    !        a(1,3) a(2,3) a(3,3)                              !
c    !               ...............                            !
c    !                               a(1,n-1) a(2,n-1) a(3,n-1) !
c    ! a(3,n)                                  a(1,n)   a(2,n)  !
c
      dimension a(6,n)
      n2 = n-2
      beta = 1./a(2,1)
      gamma = a(3,n)
      teta = a(1,1)*beta
      a(4,1) = beta
      a(5,1) = gamma
      a(6,1) = teta
      sum = gamma*teta
      do 10 i=2,n2
         v = a(3,i-1)*beta
         aa = a(1,i)
         beta = 1./(a(2,i)-aa*v)
         gamma = -gamma*v
         teta = -teta*aa*beta
         a(4,i) = beta
         a(5,i) = gamma
         a(6,i) = teta
         sum = sum+gamma*teta
  10  continue
      n1 = n-1
      v = a(3,n2)*beta
      aa = a(1,n1)
      beta = 1./(a(2,n1)-aa*v)
      gamma = a(1,n)-gamma*v
      teta = (a(3,n1)-teta*aa)*beta
      a(4,n1) = beta
      a(5,n1) = gamma
      a(6,n1) = teta
      a(4,n) = 1./(a(2,n)-(sum+gamma*teta))
      return
      end
