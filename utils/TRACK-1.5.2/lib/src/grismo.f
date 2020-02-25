      subroutine grismo(x,mx,y,my,z,mxy,tx,nx,ty,ny,p,c,ncf,iflag,
     <                  h, right, mnxmy, q, mynx, ax, ay, d, mnxy,
     <                  spx, spy, bx, by, nrx, nry)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
c  given the set of function values z(i,j) defined on the rectangular
c  grid (x(i),y(j)),i=1,2,...,mx;j=1,2,...,my, grismo determines a
c  smooth spline approximation of degree kx in x and ky in y with a set
c  of knots tx(i),i=1,2,...,nx in the x-direction and ty(j),j=1,2,...,ny
c  in the y-direction. this spline sp(x,y) depends on the value of the
c  parameter p (p>0) and has the following properties:
c    - if p tends to infinity, sp(x,y) becomes the least-squares spline
c      of given degree with given knots.
c    - if p tends to zero, sp(x,y) becomes the least-squares polynomial
c      of degree kx in x and ky in y.
c    - the function  f(p)=sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)
c      is continuous and strictly decreasing for p>0.
c  grismo offers the possibility of computing several smoothing splines
c  according to different p-values, without having to repeat all
c  computations. in that case the user must provide with a common block
c  values. the smoothing spline is given in its b-spline representation
c  and can be evaluated by means of function bisp.
c
c  calling sequence:
c     call grismo(x,mx,y,my,z,mxy,tx,nx,ty,ny,p,c,ncf,iflag)
c
c  input parameters:
c    x    : array, length mx, containing the different absciss values
c           of the rectangular grid.
c    mx   : integer, giving the number of abscissae x(i).
c    y    : array, length my, containing the different ordinate values
c           of the rectangular grid.
c    my   : integer, giving the number of ordinates y(j).
c    z    : array, minimum length mx*my, containing the data values z(i,j)
c           in the order z(1,1),z(1,2),..,z(1,my),z(2,1),..,z(2,my),..
c           z(mx,1),...,z(mx,my).
c    kx,ky: integer values, giving the degree of approximation in each
c           direction.
c    tx   : array,length nx, which contains the position of the knots in
c           the x-direction,i.e. the position of the interior knots
c           tx(kx+2),...,tx(nx-kx-1), as well as the position of the
c           knots tx(1),..,tx(kx+1) and tx(nx-kx),..,tx(nx) which are
c           needed for the b-spline representation of sp(x,y).
c    nx   : integer value, giving the total number of knots in the x-
c           direction.
c    ty   : array,lengt ny, which contains the position of the knots in
c           the y-direction,i.e. the position of the interior knots
c           ty(ky+2),..,ty(ny-ky-1), as well as the position of the
c           knots ty(1),..,ty(ky+1) and ty(ny-ky),..,ty(ny) which are
c           needed for the b-spline representation of sp(x,y).
c    ny   : integer value, containing the total number of knots in the
c           y-direction.
c    p    : real value, giving the positive smoothing parameter p.
c           p<=0  is assumed to stand for p=infinity. so in that case
c           the least-squares spline approximation is determined.
c    iflag: integer flag.
c      iflag < 0 : first call of the routine with a new set of data.
c      iflag = 0 : second call of the routine in case p<=0 at the first
c                  call.
c      iflag > 0 : second call of the routine in case p>0 at the first
c                  call. also all next calls.
c    if iflag >= 0 the user must provide with a common block values.
c
c  output parameter:
c    c    : array,length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients of sp(x,y).
c
c  restrictions:
c    1) mx > kx > 0 ; my > ky > 0.
c    2) tx(kx+1) <= x(r) < x(r+1) <= tx(nx-kx) , r=1,2,..,mx-1.
c       ty(ky+1) <= y(r) < y(r+1) <= ty(ny-ky) , r=1,2,..,my-1.
c    3) restrictions on the position of the knots; this can be verified
c       by means of subroutine check,i.e. after the calls
c             call check(x,mx,tx,nx,kx,ier)     and
c             call check(y,my,ty,ny,ky,ier)
c       the output parameter ier must be equal to zero.
c
c  other subroutines required:
c    disco, cossin, rotate and bsplin.
c  the following auxiliary arrays must have dimension specifications
c  at least  h(max(kx,ky)+2), right(max(nx,my)),q(my*nx),ax(nx,kx+1),
c  ay(ny,ky+1),d(max(nx,ny)).
c
      include 'param_interp.ptr'
c
      dimension x(mx),y(my),z(mxy),tx(nx),ty(ny),c(ncf),h(kmax2)
     < ,right(mnxmy),q(mynx),ax(nx,kx1),ay(ny,ky1),d(mnxy)
     < ,spx(mx,kx1),nrx(mx),spy(my,ky1),nry(my),bx(mnxmy,kx2)
     < ,by(mnxmy,ky2)
c  common/values/spx(mx,kx+1),nrx(mx),spy(my,ky+1),nry(my),
c                bx(max(nx,ny),kx+2),by(max(nx,ny),ky+2)
c    spx  : array,which contains the value of the non-zero b-splines
c           at each absciss x(i).
c    nrx  : integer array, which indicates which b-splines are non-zero
c           at each absciss x(i).
c    spy  : array, which contains the value of the non-zero b-splines
c           at each ordinate y(j).
c    nry  : integer array, which indicates which b-splines are non-zero
c           at each ordinate y(j).
c    bx   : array, which contains the discontinuity jumps of the
c           derivatives (order kx) of the b-splines at the interior knots
c           tx(kx+2),...,tx(nx-kx-1).
c    by   : array,which contains the discontinuity jumps of the
c           derivatives (order ky) of the b-splines at the interior knots
c           ty(ky+2),...,ty(ny-ky-1).
c
      external bsplin_c
c
c  the b-spline coefficients of the smoothing spline are calculated as
c  the least-squares solution of the over-determined linear system of
c  equations  (ay) c (ax)' = q       where
c
c               |     (spx)      |            |     (spy)      |
c        (ax) = | -------------- |     (ay) = | -------------- |
c               | sqrt(1/p) (bx) |            | sqrt(1/p) (by) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
c                b-spline coefficients.
c       z      : the my x mx matrix which contains the function values.
c       spx,spy: the mx x (nx-kx-1), resp. my x (ny-ky-1) observation
c                matrices according to the least-squares problems in
c                the x-, resp. y-direction.
c       bx,by  : the (nx-2*kx-1) x (nx-kx-1),resp. (ny-2*ky-1) x (ny-ky-1)
c                matrices which contain the discontinuity jumps of the
c                derivatives of the b-splines in the x-,resp. y-direction
      nk1x = nx-kx1
      nk1y = ny-ky1
      if(p.gt.0.) pinv = 1./p
c  it depends on the value of iflag and p whether the matrices (spx),
c  (spy), (bx) and (by) still must be determined.
      if(iflag) 10,100,120
c  calculate the non-zero elements of the matrix (spx) which is the
c  observation matrix according to the least-squares spline approximation
c  problem in the x-direction.
  10  l = kx1
      l1 = kx2
      number = 0
      do 50 it=1,mx
        arg = x(it)
  20    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 30
        l = l1
        l1 = l+1
        number = number+1
        go to 20
  30    call bsplin_c(tx,nx,kx,arg,l,h)
        do 40 i=1,kx1
          spx(it,i) = h(i)
  40    continue
        nrx(it) = number
  50  continue
c  calculate the non-zero elements of the matrix (spy) which is the
c  observation matrix according to the least-squares spline approximation
c  problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call bsplin_c(ty,ny,ky,arg,l,h)
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      if(p.le.0.) go to 120
c  calculate the non-zero elements of the matrix (bx).
 100  if(nx.eq.2*kx1) go to 110
      call disco(tx,nx,kx2,bx,mnxmy)
c  calculate the non-zero elements of the matrix (by).
 110  if(ny.eq.2*ky1) go to 120
      call disco(ty,ny,ky2,by,mnxmy)
c  reduce the matrix (ax) to a unit upper triangular form (rx) using
c  givens rotations without square roots. apply the same transformations
c  to the rows of matrix q to obtain the my x (nx-kx-1) matrix g.
c  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
c  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
        d(i) = 0.
        do 140 j=1,kx1
          ax(i,j) = 0.
 140  continue
      l = 0
      nrold = 0
c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      do 270 it=1,mx
        number = nrx(it)
 150    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 260
        ibandx = kx2
c  fetch a new row of matrix (bx).
        n1 = nrold+1
        do 160 j=1,kx2
          h(j) = bx(n1,j)
 160    continue
c  find the appropriate column of q.
        do 170 j=1,my
          right(j) = 0.
 170    continue
        wi = pinv
        irot = nrold
        go to 210
c  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
c  find the appropriate column of q.
        do 200 j=1,my
          l = l+1
          right(j) = z(l)
 200    continue
        wi = 1.
        irot = number
c  rotate the new row of matrix (ax) into triangle.
 210    do 240 i=1,ibandx
          if(wi.eq.0.) go to 250
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 240
c  calculate the parameters of the givens transformation.
          call cossin(piv,wi,d(irot),cos,sin)
c  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
            call rotate(piv,cos,sin,right(j),q(iq))
 220      continue
c  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 0
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
            call rotate(piv,cos,sin,h(j),ax(irot,i2))
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
 260    nrold = nrold+1
        go to 150
 270  continue
c  reduce the matrix (ay) to a unit upper triangular form (ry) using
c  givens rotations without square roots. apply the same transformations
c  to the columns of matrix g to obtain the (ny-ky-1) x (nx-kx-1) matrix
c  h. store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
c  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue
      do 290 i=1,nk1y
        d(i) = 0.
        do 290 j=1,ky1
          ay(i,j) = 0.
 290  continue
      nrold = 0
c  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      do 420 it=1,my
        number = nry(it)
 300    if(nrold.eq.number) go to 330
        if(p.le.0.) go to 410
        ibandy = ky2
c  fetch a new row of matrix (by).
        n1 = nrold+1
        do 310 j=1,ky2
          h(j) = by(n1,j)
 310    continue
c  find the appropiate row of g.
        do 320 j=1,nk1x
          right(j) = 0.
 320    continue
        wi = pinv
        irot = nrold
        go to 360
c  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
c  find the appropiate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        wi = 1.
        irot = number
c  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          if(wi.eq.0.) go to 400
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call cossin(piv,wi,d(irot),cos,sin)
c  apply that transformation to the colums of matrix g.
          ic = irot
          do 370 j=1,nk1x
            call rotate(piv,cos,sin,right(j),c(ic))
            ic = ic+nk1y
 370      continue
c  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 0
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
            call rotate(piv,cos,sin,h(j),ay(irot,i2))
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = nrold+1
        go to 300
 420  continue
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (ry) c (rx)' = h.
c  first step: solve the system  (ry) (c1) = h.
      k = -1
      do 450 i=1,nk1x
        k = k+nk1y
        k1 = k
        j = nk1y-1
        do 440 i1=2,nk1y
          store = c(k1)
          l2 = min0(ibandy-1,i1-1)
          l1 = k1
          do 430 l3=1,l2
            l1 = l1+1
            store = store-c(l1)*ay(j,l3)
 430      continue
          c(k1) = store
          j = j-1
          k1 = k1-1
 440    continue
 450  continue
c  second step: solve the system  c (rx)' = (c1).
      k = ncof-2*nk1y
      do 480 j=1,nk1y
        k = k+1
        k1 = k
        i = nk1x-1
        do 470 i1=2,nk1x
          store = c(k1)
          l2 = min0(ibandx-1,i1-1)
          l1 = k1
          do 460 l3=1,l2
            l1 = l1+nk1y
            store = store-c(l1)*ax(i,l3)
 460      continue
          c(k1) = store
          i = i-1
          k1 = k1-nk1y
 470    continue
 480  continue
      return
      end
