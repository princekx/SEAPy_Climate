      subroutine smoopy(x,mx,y,my,z,mxy,xb,xe,yb,ye,kxx,kyy,s,
     < nx,tx,ny,ty,nmax,c,ncof,nxest,nyest,nrdatx,nrdaty,fpintx,
     < fpinty,nplusx,nplusy,lastdi,fp0,fpold,reducx,reducy,fp,
     < iopt,ier)
c
c      implicit none
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
c  given the set of function values z(i,j) defined on the rectangular
c  grid (x(i),y(j)),i=1,2,...,mx;j=1,2,...,my, smoopy determines a
c  smooth spline approximation of degree kx in the x-direction and ky
c  in the y-direction for the interval xb<=x<=xe;yb<=y<=ye. the number
c  of knots in each direction and their position tx(i),i=1,2,...,nx;
c  ty(j),j=1,2,...,ny, is chosen automatically by the routine. the
c  smoothness of s(x,y) is achieved by minimalizing the discontinuity
c  jumps of the derivatives of the spline at the knots. the amount of
c  smoothness of sp(x,y) is determined by the condition that
c    fp=sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)) <= s,
c  with s a given non-negative constant.
c  the fit is given in its b-spline representation and can be evaluated
c  by means of function bisp.
c
c  author :
c	Paul Dierckx
c	Department of Computer Science, K.U.Leuven,
c	Celestijnenlaan 200A, B-3030 Leuven, Belgium.
c
c  references :
c	SIAM J.Numerical Analysis 19(1982),1286-1304.
c	Report TW53, Dept. Computer Science, K.U.Leuven (1980).
c
c  calling sequence:
c     call smoopy(x,mx,y,my,z,xb,xe,yb,ye,kx,ky,s,nx,tx,ny,ty,
c     c,fp,iopt,ier)
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
c    xb,xe: real values, containing the boundaries of the approximation
c           domain in the x-direction.
c    yb,ye: real values, containing the boundaries of the approximation
c           domain in the y-direction.
c    kx,ky: integer values, giving the degree of approximation in each
c           direction.
c    s    : real value, containing the smoothing factor.
c    iopt : integer flag, which takes the value 0 or 1.
c      iopt = 0: the routine will restart all computations.
c      iopt = 1: the routine will start with the knots found at the last
c                call of the routine. if iopt=1 the output parameters
c                tx,nx,ty and ny must be input parameters as well. if
c                iopt=1 the user must provide with a common block
c     common/optopt/nrdatx(nmax),nrdaty(nmax),fp0,fpold,reducx,reducy,
c     nplusx,nplusy,lastdi
c
c  output parameters:
c    tx   : array,length nmax(see data init.stat.), which contains the
c          position of the knots in the x-direction,i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1), as well as the
c          position of the knots tx(1)=tx(2)=..=tx(kx+1)=xb and tx(nx-kx)
c          =..=tx(nx)=xe which are needed for the b-spline representation
c    nx   : integer value, giving the total number of knots in the x-
c           direction.
c    ty   : array, length nmax, which contains the position of the knots
c           in the y-direction.
c    ny   : integer value, containing the total number of knots in the
c           y-direction.
c    c    : array, length ncof(see data init.stat.),which contains the
c           b-spline coefficients.
c    fp   : real value,which contains the weighted sum of squared residus
c    ier  : error message.
c      ier = 0: normal return
c      ier =-1: normal return; s(x,y) is an interpolating spline.
c      ier =-2: normal return; s(x,y) is the least-squares polynomial.
c      ier > 0: abnormal termination.
c        ier = 1: the required storage space exceeds the available space,
c                 specified by the parameters nxest and nyest( data init.)
c                 probably causes: s,nxest or nyest too small.
c        ier = 2: a theoretically impossible behaviour of the function
c                 f(p) was found during the iteration process.
c                 probably causes: tol too small.
c        ier = 3: the maximum allowable number of iterations to find the
c                 root of f(p)=s has been reached.
c                 probably causes: maxit or tol too small.
c        ier =10: some of the input data are invalid(see restrictions)
c
c  restrictions:
c    1) mx > kx > 0 ; my > ky > 0.
c    2) xb <= x(r) < x(r+1) <= xe ,r=1,2,...,mx-1.
c       yb <= y(r) < y(r+1) <= ye ,r=1,2,...,my-1.
c    3) s >= 0.
c    4) nxest >= 2*kx+2 ; nyest >= 2*ky+2.
c
c  other subroutines required:
c    grismo,sqgrid,ration,disco,cossin,rotate,bsplin and nknot.
c
      include 'param_interp.ptr'
c
      dimension x(mx),y(my),z(mxy),tx(nmax),ty(nmax),c(ncof),
     < nrdatx(nmax),nrdaty(nmax),fpintx(nmax),fpinty(nmax)
c
      external grismo_c

c
c    nrdatx : integer array, length nmax, which gives the number of
c             absciss values x(i) inside each knot interval tx(j)<x<tx(j+1)
c    nrdaty : integer array,length nmax, which gives the number of or-
c             dinate values y(i) inside each knot interval ty(j)<y<ty(j+1)
c    fp0    : real value, which contains the sum of squares of residual
c             right hand sides for the least-squares polynomial of
c             degree kx in x and ky in y.
c    fpold  : real value, which contains the sum of squares of residual
c             right hand sides for the least-squares spline which
c             corresponds to the last found set of knots but one.
c    reducx : real value, which gives the reduction in the sum of squares
c             of residual right hand sides according to the last addition
c             of knots in the x-direction.
c    reducy : real value, which gives the reduction in the sum of squares
c             of residual right hand sides according to the last addition
c             of knots in the y-direction.
c    nplusx : integer value, which contains the number of knots we
c             added in the x-direction the last time.
c    nplusy : integer value, which contains the number of knots we
c             added in the y-direction the last time.
c    lastdi : integer value, which denotes whether the last added knots
c             were located in the x- or in the y-direction.
c
c  data initialization statement to specify
c    tol  : the requested relative accuracy for the root of f(p)=s.
c    maxit: the maximum allowable number of iterations to find the root.
c    nxest: over-estimates for the numbers nx and ny. these parameters
c    nyest  must be set by the user to indicate the storage space
c           available to the routine. the dimension specifications in
c             smoopy:z(mx*my),tx,ty,nrdatx,nrdaty,fpintx,fpinty(nmax),
c                    c(ncof),
c             grismo:z(mx*my),c(ncof),h(kmax+2),right(max(nx,my)),
c                    q(my*nx),ax(nx,kx+1),ay(ny,ky+1),d(nmax),
c             sqgrid:z(mx*my),c(ncof)
c             disco :b(nmax,...),h(2*kmax+2),
c             bsplin:h(kmax+1),hh(kmax),
c             nknot :t,fpint,nrdata(nmax),
c             common/values/spx(mx,kx+1),nrx(mx),spy(my,ky+1),nry(my),
c               bx(nmax,kx+2),by(nmax,ky+2)
c           depend (implicitly) on mx,my,kx,ky,nx and ny,i.e.
c             nmax=max(nx,ny),kmax=max(kx,ky),ncof=(nx-kx-1)*(ny-ky-1)
c           since nx and ny are unknown at the time the user sets up the
c           dimension information an over-estimate of these arrays will
c           generally be made. the following remarks will help the user
c             (1) 2*kx+2<=nx<=mx+kx , 2*ky+2<=ny<=my+ky
c             (2) the smaller the value of s, the greater nx and ny
c                 will be.
c             (3) normally nx=mx/2 and ny=my/2 are over-estimates.
c  before starting computations a data check is made. if the input data
c  are invalid,controle is immediately repassed to the driver program
c
c check the compatability of the spline approx. for C and fortran
c
      kmax = max0(kx, ky)
      kmax2=kmax + 2
      k2max2=2*kmax + 2

      if(kxx .ne. kx .or. kyy .ne. ky) then
         ier = 20
         return
      endif

      ier = 0
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(x,y). c
c  if the sum f(p=inf)<=s we accept the choice of knots. otherwise we  c
c  increase their number by one. the initial set of knots depends on   c
c  the value of s and iopt                                             c
c    if s=0 we have spline interpolation; in that case the number of   c
c     knots in the x-direction equals nx=nmaxx=mx+kx+1 and in the y-   c
c     direction ny=nmaxy=my+ky+1.                                      c
c    if s>0 and                                                        c
c      iopt=0 we first compute the least-squares polynomial of degree  c
c       kx in x and ky in y: nx=nminx=2*kx+2 ; ny=nminy=2*ky+2.        c
c      iopt=1 we start with the set of knots found at the last call    c
c       of the routine, except for the case that s > fp0; then we      c
c       compute the least-squares polynomial directly.                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  find nmaxx and nmaxy which denote the number of knots in x- and y-
c  direction in case of spline interpolation.
      nmaxx = mx+kx1
      nmaxy = my+ky1
c  find nxe and nye which denote the maximum number of knots
c  allowed in each direction
      nxe = min0(nmaxx,nxest)
      nye = min0(nmaxy,nyest)
      if(s.gt.tols) goto 100
c  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
c  test whether the required storage space exceeds the available one.
      if(ny.gt.nyest .or. nx.gt.nxest) go to 420
c  find the position of the interior knots in case of interpolation.
c  the knots in the x-direction.
      mk1 = mx-kx1
      if(mk1.eq.0) go to 60
      k3 = kx/2
      i = kx1+1
      j = k3+2
      if(k3*2.eq.kx) go to 40
      do 30 l=1,mk1
        tx(i) = x(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
  40  do 50 l=1,mk1
        tx(i) = (x(j)-x(j-1))*0.5
        i = i+1
        j = j+1
  50  continue
c  the knots in the y-direction.
  60  mk1 = my-ky1
      if(mk1.eq.0) go to 120
      k3 = ky/2
      i = ky1+1
      j = k3+2
      if(k3*2.eq.ky) go to 80
      do 70 l=1,mk1
        ty(i) = y(j)
        i = i+1
        j = j+1
  70  continue
      go to 120
  80  do 90 l=1,mk1
        ty(i) = (y(j)-y(j-1))*0.5
        i = i+1
        j = j+1
  90  continue
      go to 120
c  if s > 0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s >= fp0, we start computing the least-squares
c  polynomial of degree kx in x and ky in y (which is a spline without
c  interior knots). if iopt=1 and fp0 > s we start computing the least-
c  squares spline according to the set of knots found at the last call
c  of the routine.
 100  if(iopt.le.0) go to 110
      if(fp0.gt.s) go to 120
 110  nx = nminx
      ny = nminy
      nrdatx(1) = mx-2
      nrdaty(1) = my-2
      lastdi = 0
 120  mm = mx+my
      iflag = -1
      p = -1.
c  main loop for the different sets of knots. mm=mx+my is a save upper
c  bound for the number of trials.
      do 270 iter=1,mm
        if(nx.eq.nminx .and. ny.eq.nminy) ier = -2
c  find nrintx (nrinty) which is the number of knot intervals in the
c  x-direction (y-direction).
        nrintx = nx-nminx+1
        nrinty = ny-nminy+1
c  find ncof, the number of b-spline coefficients for the current set
c  of knots.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        i = nx
        do 130 j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
 130    continue
        i = ny
        do 140 j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
 140    continue
c  find the least-squares spline sinf(x,y).
        call grismo_c(z,tx,nx,ty,ny,p,c,fp,fpintx,
     <                 fpinty,iflag)
        if(ier.eq.-2) fp0 = fp
c  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430
c  test whether the least-squares spline is an acceptable solution.
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  increase the number of knots.
c  if nx=nxe and ny=nye we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nx.eq.nxe .and. ny.eq.nye) go to 420
        ier = 0
c  find the maximum sum of squared residuals for the knot intervals
c  in the x-direction.
        fpmaxx = 0.
        do 150 i=1,nrintx
          if(fpintx(i).gt.fpmaxx) fpmaxx = fpintx(i)
 150    continue
c  find the maximum sum of squared residuals for the knot intervals
c  in the y-direction.
        fpmaxy = 0.
        do 160 i=1,nrinty
          if(fpinty(i).gt.fpmaxy) fpmaxy = fpinty(i)
 160    continue
c  adjust the parameter reducx or reducy according to the direction
c  in which the last added knots were located.
        if(lastdi) 170,190,180
 170    reducx = fpold-fp
        go to 190
 180    reducy = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 190    fpold = fp
c  test whether we are going to add knots in the x- or in the y-direction
        if((fpmaxx.lt.fpmaxy .and. ny.lt.nye) .or. nx.eq.nxe)
     <    go to 230
c  find nplusx, the number of knots we are going to add in the x-direction
        if(nx.gt.nminx) go to 200
        nplusx = 1
        go to 210
 200    npl1 = nplusx*2
        if(reducx.gt.acc) npl1 = dble(nplusx)*fpms/reducx
        nplusx = min0(nplusx*2,max0(npl1,nplusx/2,1))
 210    lastdi = -1
        do 220 l=1,nplusx
c  add a new knot in the x-direction
          call nknot(x,mx,tx,nx,nmax,fpintx,nrdatx,nrintx)
c  test whether we cannot further increase the number of knots in the
c  x-direction.
          if(nx.eq.nxe) go to 270
 220    continue
        go to 270
c  find nplusy, the number of knots we are going to add in the y-direction
 230    if(ny.gt.nminy) go to 240
        nplusy = 1
        go to 250
 240    npl1 = nplusy*2
        if(reducy.gt.acc) npl1 = dble(nplusy)*fpms/reducy
        nplusy = min0(nplusy*2,max0(npl1,nplusy/2,1))
 250    lastdi = 1
        do 260 l=1,nplusy
c  add a new knot in the y-direction.
          call nknot(y,my,ty,ny,nmax,fpinty,nrdaty,nrinty)
c  test whether we cannot further increase the number of knots in the
c  y-direction.
          if(ny.eq.nye) go to 270
 260    continue
c  restart the computations with the new set of knots.
 270  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.-2) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
c  this smoothing spline varies with the parameter p in such a way thatc
c    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. iteratively we then have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -1.
      f3 = fpms
      p = -f1/f3
      iflag = 0
      icheck = 0
c  iteration process to find the root of f(p)=s.
      do 320 iter = 1,maxit
c  find the smoothing spline sp(x,y).
        call grismo_c(z,tx,nx,ty,ny,p,c,fp,fpintx,
     <                 fpinty,iflag)
        iflag = 1
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(icheck.ne.0) go to 310
        if((f2-f3).gt.acc) go to 305
c  our initial choice of p is too large.
        p = p*0.1e-02
        p3 = p2
        f3 = f2
        go to 320
 305    if((f1-f2).gt.acc) go to 310
c  our initial choice of p is too small
        p = p*0.1e+04
        p1 = p2
        f1 = f2
        go to 320
c  test whether the iteration process proceeds as theoretically
c  expected.
 310    if(f2.ge.f1 .or. f2.le.f3) go to 410
        icheck = 1
c  find the new value of p.
        p = ration(p1,f1,p2,f2,p3,f3)
 320  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
