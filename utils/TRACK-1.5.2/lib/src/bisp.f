      subroutine bisp(bp,d1,tx,nx,ty,ny,c,ncof,x,y,idif)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
c
c  function bisp evaluates a two-dimensional spline function s(x,y)
c  which is given in its b-spline representation.
c
c  calling sequence:
c     value = bisp(bp,d1,tx,nx,ty,ny,c,ncof,x,y,idif)
c
c  input parameters:
c    bp   : value of spline at the point (x,y)
c    d1   : value of derivatives at the point (x,y)
c    tx   : array,length nx, which contains the position of the knots in
c           the x-direction,i.e. the position of the interior knots
c           tx(kx+2),...,tx(nx-kx-1), as well as the position of the
c           knots tx(1),..,tx(kx+1) and tx(nx-kx),..,tx(nx) which are
c           needed for the b-spline representation of s(x,y).
c    nx   : integer value, giving the total number of knots in the x-
c           direction.
c    ty   : array,lengt ny, which contains the position of the knots in
c           the y-direction,i.e. the position of the interior knots
c           ty(ky+2),..,ty(ny-ky-1), as well as the position of the
c           knots ty(1),..,ty(ky+1) and ty(ny-ky),..,ty(ny) which are
c           needed for the b-spline representation of s(x,y).
c    ny   : integer value, containing the total number of knots in the
c           y-direction.
c    kx,ky: integer values, giving the degree of the spline in each
c           direction.
c    c    : array,length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients of s(x,y).
c    x,y  : real values,which contain the coordinates of the point where
c           s(x,y) must be evaluated.
c    idif : idif = 1; calculate 1st derivatives, idif = 2; calulate
c           dxx and dyy as well.
c
c  output parameter:
c    value: real variable, which contains the value of s(x,y).
c
c  restrictions:
c      tx(kx+1) <= x <= tx(nx-kx)  ;  ty(ky+1) <= y <= ty(ny-ky)
c    if these conditions are not fulfilled, value is put to zero.
c
c  other subroutines required: bsplin.
c
c  the following auxiliary arrays must have dimension specifications at
c  least  hx(kx+1), hy(ky+1).
c
      include 'param_interp.ptr'
c
      dimension tx(nx),ty(ny),c(ncof),hx(kx1,kx1),hy(ky1,ky1),
     < cc(kmax2),ccx(kx1,ky1), ccy(kx1,ky1), d1(4)
c
      nk1x = nx-kx1
      nk1y = ny-ky1
      bp = 0.
c  test whether the absciss of the given point lies inside the limits
c  tx(kx+1) and tx(nx-kx).
      if(x.lt.tx(kx1) .or. x.gt.tx(nk1x+1))then
         if((tx(kx1) - x) .le. tolkn) then
            x = tx(kx1)
         else if ((x - tx(nk1x+1)) .le. tolkn) then
            x = tx(nk1x+1)
         else
            write(nout, 900)
            go to 400
         endif
      endif
c  test whether the ordinate of the given point lies inside the limits
c  ty(ky+1) and ty(ny-ky).
      if(y.lt.ty(ky1) .or. y.gt.ty(nk1y+1))then
         if((ty(ky1) - y) .le. tolkn) then
            y = tx(ky1)
         else if ((y - tx(nk1y+1)) .le. tolkn) then
            y = tx(nk1y+1)
         else
            write(nout, 910)
            go to 400
         endif
      endif
      lx = kx1
      l1 = lx+1
      ly = ky1
      l2 = ly+1
c  search for knot interval tx(lx) <= x <  tx(lx+1)
  10  if(x.lt.tx(l1) .or. lx.eq.nk1x) go to 20
      lx = l1
      l1 = lx+1
      go to 10
c  search for knot interval ty(ly) <= y <  ty(ly+1).
  20  if(y.lt.ty(l2) .or. ly.eq.nk1y) go to 30
      ly = l2
      l2 = ly+1
      go to 20
c  evaluate the (kx+1) non-zero b-splines at x.
  30  call bsplin_all(tx,nx,kx,x,lx,hx)
c  evaluate the (ky+1) non-zero b-splines at y.
      call bsplin_all(ty,ny,ky,y,ly,hy)
c  find the value of s(x,y) by making the sum of the cross products of
c  the non-zero b-splines in x- and y-direction, multiplied with the
c  appropriate coefficients.
      i1 = (lx-kx1)*nk1y+ly-ky1
      do 50 i=1,kx1
        hxi = hx(i,kx1)
        j1 = i1
        do 40 j=1,ky1
          j1 = j1+1
          bp = bp+hxi*hy(j,ky1)*c(j1)
          ccx(i,j) = c(j1)
          ccy(i,j) = c(j1)
  40    continue
        i1 = i1+nk1y
  50  continue
c
c calculate all first and second order (dxx, dyy) derivatives if required
c
      if(idif .ge. 1) then
c calculate coefficient differance table in x- direction
         k = kx
         do 60 j=1, ky1
           do 65 i=1, kx1
              cc(i) = ccx(i,j)
  65       continue
           l = lx - kx
           do 70 i=2, kx1
              l = l + 1
              fac = tx(l+k) - tx(l)
              if(fac .le. 0.) go to 70
              ccx(i,j) = (cc(i) - cc(i-1))*k/fac
  70       continue
  60     continue
c compute first x- derivative
         d1(1) = 0.
         do 80 i=1,k
            hxi = hx(i,k)
            ii = i + 1
            do 90 j=1,ky1
              d1(1) = d1(1) + hxi*hy(j,ky1)*ccx(ii,j)
  90       continue
  80     continue
c calculate second derivative in the x- direction
         if(idif .eq. 2) then
            k = kx - 1
            do 100 j=1, ky1
              do 110 i=1, kx1
                 cc(i) = ccx(i,j)
  110         continue
              l = lx - kx
              do 120 i=2, kx1
                 l = l + 1
                 fac = tx(l+k) - tx(l)
                 if(fac .le. 0.) go to 120
                 ccx(i,j) = (cc(i) - cc(i-1))*k/fac
c                 print*,i,j,ccx(i,j)
  120         continue
  100        continue
c compute second x- derivative
            d1(3) = 0.
            do 130 i=1,k
               hxi = hx(i,k)
               ii = i + 2
               do 140 j=1,ky1
                 d1(3) = d1(3) + hxi*hy(j,ky1)*ccx(ii,j)
c                 print*,hxi, hy(j,ky1), ccx(ii,j)
  140          continue
  130       continue
         endif
c calculate coefficient differance table in y- direction
         k = ky
         do 200 j=1, kx1
           do 210 i=1, ky1
              cc(i) = ccy(j,i)
 210       continue
           l = ly - ky
           do 220 i=2, ky1
              l = l + 1
              fac = ty(l+k) - ty(l)
              if(fac .le. 0.) go to 220
              ccy(j,i) = (cc(i) - cc(i-1))*k/fac
 220       continue
 200     continue
c compute first y- derivative
         d1(2) = 0.
         do 230 i=1,k
            hyi = hy(i,k)
            ii = i + 1
            do 240 j=1,kx1
               d1(2) = d1(2) + hyi*hx(j,kx1)*ccy(j,ii)
 240        continue
 230     continue 
c calculate second derivative in the y- direction
         if(idif .eq. 2) then
            k = ky - 1
            do 300 j=1, kx1
              do 310 i=1, ky1
                 cc(i) = ccy(j,i)
 310          continue
              l = ly - ky
              do 320 i=2, ky1
                 l = l + 1
                 fac = ty(l+k) - ty(l)
                 if(fac .le. 0.) go to 320
                 ccy(j,i) = (cc(i) - cc(i-1))*k/fac
 320          continue
 300        continue
c compute second y- derivative
            d1(4) = 0.
            do 330 i=1,k
               hyi = hy(i,k)
               ii = i + 2
               do 340 j=1,kx1
                 d1(4) = d1(4) + hyi*hx(j,kx1)*ccy(j,ii)
 340           continue
 330        continue
         endif
      endif
 400  return
 900  format("***Warning***, ordinate outside of knot limits for" /
     +       "               B-spline evaluation in X direction." )
 910  format("***Warning***, ordinate outside of knot limits for" /
     +       "               B-spline evaluation in Y direction.")
      end
