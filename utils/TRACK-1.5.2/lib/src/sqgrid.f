      subroutine sqgrid(mx,my,z,mxy,nx,ny,c,ncof,sq,fpx,fpy,mnxmy,
     <                  spx,spy,bx,by,nrx,nry)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
c  given the set of function values z(i,j) defined on the rectangular
c  grid (x(i),y(j)),i=1,2,...,mx;j=1,2,...,my and given also a spline
c  approximation s(x,y) of degree kx in x and ky in y with a set of knots
c  tx(i),i=1,2,...,nx in the x-direction and ty(j),j=1,2,...,ny in the
c  y-direction, sqgrid calculates the quantities
c    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
c    sq = sumi=1,mx(sumj=1,my(res(i,j)))
c    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
c                  tx(r+kx) <= x(i) <= tx(r+kx+1)
c    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
c                  ty(r+ky) <= y(j) <= ty(r+ky+1)
c  for this purpose sqgrid makes use of the information stored in the
c  common block values, i.e. the value of the non-zero b-splines at the
c  different grid points.
c
c  calling sequence:
c     call sqgrid(mx,my,z,kx,ky,nx,ny,c,sq,fpx,fpy)
c
c  input parameters:
c    mx,my,z,kx,ky,nx,ny,c: see grismo.
c
c  output parameters:
c    sq  : real value, giving the total sum of squared residuals.
c    fpx : array,length nx, giving the sum of squared residuals for the
c          different knot intervals in the x-direction.
c    fpy : array,length ny, giving the sum of squared residuals for the
c          different knot intervals in the y-direction.
c
c
      include 'param_interp.ptr'
c
      dimension z(mxy),c(ncof),fpx(nx),fpy(ny),spx(mx,kx1),nrx(mx),
     <          spy(my,ky1),nry(my),bx(mnxmy,kx2),by(mnxmy,ky2)
c  initialization
      sq = 0.
      do 10 i=1,nx
        fpx(i) = 0.
  10  continue
      do 20 i=1,ny
        fpy(i) = 0.
  20  continue
      nk1y = ny-ky1
      iz = 0
      nroldx = 0
c  main loop for the different grid points.
      do 70 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 60 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
c  evaluate s(x,y) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (x,y), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
          do 40 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
            do 30 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
  30        continue
            k1 = k1+nk1y
  40      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          sq = sq+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*0.5
          if(numy.eq.nroldy) go to 50
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
  50      nroldy = numy
          if(numx.eq.nroldx) go to 60
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
  60    continue
        nroldx = numx
  70  continue
      return
      end
