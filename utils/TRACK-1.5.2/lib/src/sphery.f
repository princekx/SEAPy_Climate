c
      subroutine sphery(u,mu,v,mv,r,muv,vb,ve,mumin,r0,r1,s,nu,
     1 tu,nv,tv,nmax,c,ncof,nuest,nvest,nrdatu,nrdatv,fpintu,fpintv,
     2 nplusu, nplusv,lastdi,lastu0,lastu1,fp0,fpold,reducu,reducv,
     3 step,dr,fp,iopt,ider,ier)
c      implicit none
      implicit double precision (A-H, O-Z)
      implicit integer (i-n)

c  given the function values r(i,j) on the latitude-longitude grid
c  (u(i),v(j)), i=1,...,mu ; j=1,...,mv , sphery determines a smooth
c  bicubic spline approximation on the rectangular domain 0<=u<=pi,
c  0<=v<=2*pi. this approximation s(u,v) will satisfy the properties
c
c    (1) s(0,v) = s(0,0) = dr(1)
c
c        d s(0,v)           d s(0,0)           d s(0,pi/2)
c    (2) -------- = cos(v)* -------- + sin(v)* -----------
c        d u                d u                d u
c
c                 = cos(v)*dr(2)+sin(v)*dr(3)
c                                                      0 <= v <= 2*pi
c    (3) s(pi,v) = s(pi,0) = dr(4)
c
c        d s(pi,v)           d s(pi,0)           d s(pi,pi/2)
c    (4) -------- = cos(v)*  --------- + sin(v)* ------------
c        d u                 d u                 d u
c
c                 = cos(v)*dr(5)+sin(v)*dr(6)
c
c  and will be periodic in the variable v, i.e.
c
c         j          j
c        d s(u,0)   d s(u,2*pi)
c    (5) -------- = -----------   0 <=u<= pi , j=0,1,2
c           j          j
c        d v        d v
c
c  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
c  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
c  smoothness of s(u,v) is achieved by minimalizing the discontinuity
c  jumps of the derivatives of the spline at the knots. the amount of
c  smoothness of s(u,v) is determined by the condition that
c  fp=sumi=1,mu(sumj=1,mv((r(i,j)-s(u(i),v(j)))**2))+(r0-s(0,0))**2
c  + (r1-s(pi,0))**2 <= s, with s a given non-negative constant.
c  the fit s(u,v) is given in its b-spline representation and can be
c  evaluated by means of function bisp.
c
c  author :
c	Paul Dierckx
c	Department of Computer Science, K.U.Leuven,
c	Celestijnenlaan 200A, B-3030 Leuven, Belgium.
c
c  references :
c	Report TW73, Dept. Computer Science, K.U.Leuven (1985).
c
c  calling sequence:
c     call sphery(u,mu,v,mv,r,r0,r1,iopt,ider,s,nu,tu,nv,tv,c,fp,ier)
c
c  input parameters:
c    u    : array, length mu, containing the different u-values
c           of the radius-angle grid.
c    mu   : integer, giving the number of latitude values u(i).
c           if mu exceeds 50, see data initialisation statement.
c    v    : array, length mv, containing the different v-values
c           of the radius-angle grid. v(1) must be 0.
c    mv   : integer, giving the number of longitude values v(j).
c           if mv exceeds 50, see data initialisation statement.
c    r    : array, minimum length mu*mv, containing the values r(i,j)
c           in the order r(1,1),r(1,2),..,r(1,mv),r(2,1),..,r(2,mv),...
c           r(mu,1),...,r(mu,mv).
c    r0   : real value, containing the data value at the pole u=0.
c           if (ider(1) = 0 or 1 ).
c    r1   : real value, containing the data value at the pole u=pi.
c           if (ider(3) = 0 or 1 ).
c    iopt : integer array, specifying different options.
c     iopt(3): integer flag which specifies the requested order of
c              continuity at the pole u=pi,i.e.
c               if iopt(2)=0 condition (3) must be fulfilled,
c               if iopt(2)=1 conditions (3)+(4) must be fulfilled.
c     iopt(2): integer flag which specifies the requested order of
c              continuity at the pole u=0.,i.e.
c               if iopt(2)=0 condition (1) must be fulfilled,
c               if iopt(2)=1 conditions (1)+(2) must be fulfilled.
c     iopt(1): integer flag which specifies the mode of computation,i.e.
c               if iopt(1)=0 the routine will restart all computations,
c               if iopt(1)=1 the routine will start with the knots found
c                at the last call. in that case the output parameters
c                tu,nu,tv,nv must be input parameters as well and
c                the user must provide the common block sphsph.
c    ider : integer array, specifying the information at the poles.
c     ider(1): integer flag which specifies whether (ider(1)=0 or 1) or
c              not (ider(1)=-1) there is a data value r0 at the pole u=0
c              if ider(1)=1, r0 will be considered to be the right func-
c              tion value, and it will be fitted exactly ( dr(1)=r0 ).
c              if ider(1)=0, r0 will be considered as a data value just
c              like the other values r(i,j).
c     ider(2): integer flag which specifies whether (ider(2)=1) or not
c              (ider(2)=0) the approximation has vanishing derivatives
c              dr(2) and dr(3) at the pole u=0.
c     ider(3): integer flag which specifies whether (ider(3)=0 or 1) or
c              not (ider(3)=-1)there is a data value r1 at the pole u=pi
c              if ider(3)=1, r1 will be considered to be the right func-
c              tion value, and it will be fitted exactly ( dr(4)=r1 ).
c              if ider(3)=0, r1 will be considered as a data value just
c              like the other values r(i,j).
c     ider(4): integer flag which specifies whether (ider(4)=1) or not
c              (ider(4)=0) the approximation has vanishing derivatives
c              dr(5) and dr(6) at the pole u=pi.
c    s    : real value, containing the smoothing factor.
c
c  output parameters:
c    tu   : array,length nmax(see data init.stat.), which contains the
c           position of the knots in the u-direction,i.e. the position
c           of the interior knots tu(5),...,tu(nu-4), as well as the
c           position of the boundary knots tu(1),...tu(4) and tu(nu-3),
c           ...tu(nu) which are needed for the b-spline representation.
c    nu   : integer value, giving the total number of knots in the u-
c           direction.
c    tv   : array, length nmax, which contains the position of the knots
c           in the v-direction.
c    nv   : integer value, containing the total number of knots in the
c           v-direction.
c    c    : array, length ncof(see data init.stat.),which contains the
c           b-spline coefficients.
c    fp   : real value,containing the sum of squared residuals.
c    ier  : error message.
c      ier = 0: normal return
c      ier =-1: normal return; s(u,v) is an interpolating spline.
c      ier =-2: normal return; s(u,v) is the least-squares polynomial.
c      ier > 0: abnormal termination.
c        ier = 1: the required storage space exceeds the available space
c                 specified by nuest and nvest (see data initialization)
c                 probably causes: s,nuest or nvest too small.
c        ier = 2: a theoretically impossible behaviour of the function
c                 f(p) was found during the iteration process.
c                 probably causes: tol too small.
c        ier = 3: the maximum allowable number of iterations to find the
c                 root of f(p)=s has been reached.
c                 probably causes: maxit or tol too small.
c        ier =10: some of the input data are invalid (see restrictions).
c
c  restrictions:
c    1) 0 <= iopt(2) <= 1 ; 0 <= iopt(3) <= 1.
c    2) -1 <= ider(1) <= 1 ; 0 <= ider(2) <= 1.
c       -1 <= ider(3) <= 1 ; 0 <= ider(4) <= 1.
c    3) mu >= mumin ; mv >= 4.
c          mumin = max(1,mu0+mu1)
c          mu0 = 1-ider(2) if ider(1) >= 0
c              = 2-ider(2) if ider(1) <  0
c          mu1 = 1-ider(4) if ider(3) >= 0
c              = 2-ider(4) if ider(3) <  0
c    4) 0 < u(i) < u(i+1) < pi ,i=1,2,...,mu-1.
c    5) v(j) < v(j+1) < 2*pi , j=1,2,...,mv-1.
c       v(1) = 0
c    6) s >= 0.
c    7) nuest >= 8 ; nvest >= 11.
c
c  other subroutines required:
c    optsph,grisph,symsys,sqgrid,ration,cossin,rotate,bsplin,knot,
c    cytri1 and cytri2
c
      include 'param_interp.ptr'
c
      dimension u(mu),v(mv),r(muv),tu(nmax),tv(nmax),c(ncof),
     1 nrdatu(nmax),nrdatv(nmax),fpintu(nmax),fpintv(nmax),iopt(3),
     2 ider(4),idd(4),dr(6),drr(6),step(2)

       external sqgrid_c

c   common/sphsph/nrdatu(nmax),nrdatv(nmax),fp0,fpold,reducu,reducv,
c     nplusu,nplusv,lastdi,step(2),lastu0,lastu1,dr(6)
c    nrdatu : integer array, length nmax, which gives the number of
c             values u(i) inside each knot interval tu(j)<u<tu(j+1).
c    nrdatv : integer array,length nmax, which gives the number of
c             values v(i) inside each knot interval tv(j)<v<tv(j+1).
c    fp0    : real value, which contains the sum of squares of residual
c             right hand sides for the least-squares polynomial which
c             satisfies the boundary constraints.
c    fpold  : real value, which contains the sum of squares of residual
c             right hand sides for the least-squares spline which
c             corresponds to the last found set of knots but one.
c    reducu : real value, which gives the reduction in the sum of
c             squares of residual right hand sides according to the
c             last addition of knots in the u-direction.
c    reducv : real value, which gives the reduction in the sum of
c             squares of residual right hand sides according to the
c             last addition of knots in the v-direction.
c    nplusu : integer value, which contains the number of knots we
c             added in the u-direction the last time.
c    nplusv : integer value, which contains the number of knots we
c             added in the v-direction the last time.
c    lastdi : integer value, which denotes whether the last added knots
c             were located in the u- or in the v-direction.
c    step(1): real value, indicating the range (r0-step(1),r0+step(1))
c             of possible values for dr(1).
c    step(2): real value, indicating the range (r1-step(2),r1+step(2))
c             of possible values for dr(4).
c    lastu0 : integer value, indicating the number of data points used
c             to calculate the parameter step(1).
c    lastu1 : integer value, indicating the number of data points used
c             to calculate the parameter step(2).
c    dr     : real array, containing the derivative values dr(i) of
c             the least-squares spline with the current set of knots.
c  data initialization statement to specify
c    tol  : the requested relative accuracy for the root of f(p)=s.
c    maxit: the maximum allowable number of iterations to find the root.
c    nuest: over-estimates for the numbers nu and nv. these parameters
c    nvest  must be set by the user to indicate the storage space
c           available to the routine. the dimension specifications in
c             sphery:r(mu*mv),tu,tv,nrdatu,nrdatv,fpintu,fpintv(nmax),
c                    c(ncof),
c             optsph:r(mu*mv),c(ncof),
c             grisph:r(mu*mv),c(ncof),right(max(nu,mvv)),d(nu),dd(nv),
c                    q(nu*mvv),au(4,nu),av1(6,nv),av2(4,nv),c0(nv),
c                    a0(2,mv),b0(2,mv),c1(nv),a1(2,mv),b1(2,mv)
c             sqgrid:r(mu*mv),c(ncof)
c             knot :t,fpint,nrdata(nmax),
c             common/values/spu(mu,4),nru(mu),spv(mv,4),nrv(mv),
c               bu(nu,5),bv(nv,5),cosi(2,nv)
c           depend (implicitly) on mu,mv,nu and nv,i.e.
c             nmax=max(nu,nv),ncof=(nu-4)*(nv-4),mvv = mv+nv
c           since nu and nv are unknown at the time the user sets up the
c           dimension information an over-estimate of these arrays will
c           generally be made. the following remarks will help the user
c             (1) 8<=nu<=mu+5+iopt(2)+iopt(3) , 8<=nv<=mv+7
c             (2) the smaller the value of s, the greater nu and nv
c                 will be.
c             (3) normally nu=mu/2 and nv=mv/2 are over-estimates.

c  before starting computations a data check is made. if the input data
c  are invalid,controle is immediately repassed to the driver program
c

      kmax = max0(kx, ky)
      kmax2=kmax + 2
      k2max2=2*kmax + 2

      if(kx .ne. 3 .or. ky .ne. 3) then
         ier = 100
         go to 440
      endif
c
c  if not given,we compute an estimate for r0
      m = muv
      if(ider(1).ge.0) go to 30
      r0 = 0.
      do 25 i=1,mv
         r0 = r0+r(i)
  25  continue
      r0 = r0/dble(mv)
c  if not given,we compute an estimate for r1
  30  if(ider(3).ge.0) go to 40
      r1 = 0.
      j = m
      do 35 i=1,mv
         r1 = r1+r(j)
         j = j-1
  35  continue
      r1 = r1/dble(mv)
c  determine the range of r-values
  40  rmin = r0
      rmax = r1
      do 45 i=1,m
         if(r(i).lt.rmin) rmin = r(i)
         if(r(i).gt.rmax) rmax = r(i)
  45  continue
c  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+6+iopt(2)+iopt(3)
      nvmax = mv+7
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(s .lt. 0.)then
        goto 440 
      else if(s .lt. 0.001) then
        go to 50
      else
        go to 100
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(u,v). c
c  if the sum f(p=inf)<=s we accept the choice of knots. otherwise we  c
c  increase their number . the initial set of knots depends on the     c
c  value of s and iopt                                                 c
c    if s=0 we have spline interpolation; in that case the number of   c
c     knots in the u-direction equals nu=numax=mu+6+iopt(2)+iopt(3)    c
c     and in the v-direction nv=nvmax=mv+7.                            c
c    if s>0 and                                                        c
c      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
c       spline without interior knots : nu=8 ; nv=8.                   c
c      iopt(1)=1 we start with the set of knots found at the last call c
c       of the routine, except for the case that s > fp0; then we      c
c       compute the least-squares polynomial directly.                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  if s = 0, s(u,v) is an interpolating spline.
  50  ier = -1
      nu = numax
      nv = nvmax
c  test whether the required storage space exceeds the available one.
      if(nu.gt.nuest .or. nv.gt.nvest) go to 420
c  find the position of the knots in the v-direction.
      do 55 l=1,mv
        tv(l+3) = v(l)
  55  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 60 i=1,3
         tv(i) = v(l1)-pi2
         tv(l2) = v(i+1)+pi2
         l1 = l1+1
         l2 = l2+1
  60  continue
c  if not all the derivative values dr(i) are given, we will first
c  estimate these values by computing a least-squares spline
      dr(1) = r0
      dr(2) = 0.
      dr(3) = 0.
      dr(4) = r1
      dr(5) = 0.
      dr(6) = 0.
      idd(1) = ider(1)
      if(idd(1).eq.0) idd(1) = 1
      idd(2) = ider(2)
      idd(3) = ider(3)
      if(idd(3).eq.0) idd(3) = 1
      idd(4) = ider(4)
      if(ider(1).lt.0 .or. ider(3).lt.0) go to 65
      if(iopt(2).ne.0 .and. ider(2).eq.0) go to 65
      if(iopt(3).eq.0 .or. ider(4).ne.0) go to 85
c we set up the knots in the u-direction for computing the least-squares
c spline.
  65  i1 = 3
      i2 = mu-2
      nu = 4
      do 70 i=1,mu
         if(i1.gt.i2) go to 75
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  70  continue
  75  do 80 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = pi
  80  continue
      step(1) = rmax-rmin
      step(2) = step(1)
c we compute the least-squares spline for estimating the derivatives.
      ppp=-1.
      call optsph(u,mu,v,mv,r,muv,r0,r1,tu,nu,tv,nv,ppp,iopt,idd,-1,
     1  step,dr,c,ncof,fp)
c if all the derivatives at the poles are known, we compute the
c interpolating spline.
c we set up the knots in the u-direction, needed for interpolation.
  85  nn = numax-8
      if(nn.eq.0) go to 97
      ju = 2-iopt(2)
      do 90 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  90  continue
      nu = numax
      l = nu
      do 93 i=1,4
         tu(i) = 0.
         tu(l) = pi
         l = l-1
  93  continue
c we compute the interpolating spline.
      ppp = -1.
  97  call optsph(u,mu,v,mv,r,muv,r0,r1,tu,nu,tv,nv,ppp,iopt,idd,-1,
     1 step,dr,c,ncof,fp)
      lastu0 = mu
      lastu1 = mu
      go to 440
c  if s > 0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s >= fp0,we start computing the least-squares
c  polynomial. if iopt=1 and fp0 > s we start computing the least-
c  squares spline according to the set of knots found at the last call
c  of the routine.
 100  ier = 0
c  find nue and nve which denote the maximum number of knots
c  allowed in each direction
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)
      if(iopt(1).eq.0) go to 105
      step(1) = -step(1)
      step(2) = -step(2)
      if(fp0.gt.s) go to 110
 105  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      idd(3) = ider(3)
      idd(4) = 1
      dr(1) = r0
      dr(2) = 0.
      dr(3) = 0.
      dr(4) = r1
      dr(5) = 0.
      dr(6) = 0.
      step(1) = rmax-rmin
      step(2) = step(1)
      nu = 8
      nv = 8
      nrdatu(1) = mu-2+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      lastu0 = mu
      lastu1 = mu
c  main loop for the different sets of knots. mm=mu+mv is a save upper
c  bound for the number of trials.
 110  mm = mu+mv
      iflag = -1
      p = -1.
      do 270 iter=1,mm
c  find nrintu (nrintv) which is the number of knot intervals in the
c  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
        i = nu
        do 115 j=1,4
          tu(j) = 0.
          tu(i) = pi
          i = i-1
 115    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 120 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-pi2
          tv(l3) = tv(l1)+pi2
 120    continue
        if(nrintu.eq.1) go to 140
c  find an estimate of the range of possible values for the optimal
c  derivatives at the pole u=0.
        ktu = nrdatu(1)+2-iopt(2)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu0) go to 130
         rmin = r0
         rmax = r1
         l = mv*ktu
         do 125 i=1,l
            if(r(i).lt.rmin) rmin = r(i)
            if(r(i).gt.rmax) rmax = r(i)
 125     continue
         step(1) = rmax-rmin
         lastu0 = ktu
c  find an estimate of the range of possible values for the optimal
c  derivatives at the pole u=pi.
 130    ktu = nrdatu(nrintu)+2-iopt(3)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu1) go to 140
         rmin = r0
         rmax = r1
         l = mv*ktu
         j = m
         do 135 i=1,l
            if(r(j).lt.rmin) rmin = r(j)
            if(r(j).gt.rmax) rmax = r(j)
            j = j-1
 135     continue
         step(2) = rmax-rmin
         lastu1 = ktu
c  find the least-squares spline sinf(u,v).
 140    call optsph(u,mu,v,mv,r,muv,r0,r1,tu,nu,tv,nv,p,iopt,idd,iflag,
     2  step,dr,c,ncof,fp)
        if(step(1).lt.0.) step(1) = -step(1)
        if(step(2).lt.0.) step(2) = -step(2)
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  increase the number of knots.
c  if nu=nue and nv=nve we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
c  calculate for each knot interval tu(j+3)<=u<=tu(j+4) (tv(j+3)<=
c  v<=tv(j+4)) the sum of squared residuals fpintu(j),j=1,2,...,nu-7
c  (fpintv(j),j=1,2,...,nv-7) for the data points having their u-
c  value (v-value) belonging to that interval.
        call sqgrid_c(r,nu,nv,c,sq,fpintu,fpintv)
        if(ider(1).eq.0) fpintu(1) = fpintu(1)+(r0-dr(1))**2
        if(ider(3).eq.0) fpintu(nrintu) = fpintu(nrintu)+(r1-dr(4))**2
        ier = 0
c  adjust the parameter reducu or reducv according to the direction
c  in which the last added knots were located.
        if(lastdi .lt. 0)then
         go to 160
        else if(lastdi .eq. 0) then
         go to 150
        else
         go to 170
        endif
 150     nplv = 3
         idd(2) = ider(2)
         idd(4) = ider(4)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.8) go to 180
        npl1 = nplusu*2
        if(reducu.gt.acc) npl1 = dble(nplusu)*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv.eq.8) go to 190
        npl1 = nplusv*2
        if(reducv.gt.acc) npl1 = dble(nplusv)*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
c  test whether we are going to add knots in the u- or v-direction.
 190    if(nplu-nplv .lt. 0) then
          go to 210 
        else if(nplu-nplv .gt. 0) then
          go to 230
        endif
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
         if(lastdi.ge.0) lastdi = 0
         lastdi = lastdi-1
         if(lastdi.gt.(-5)) go to 215
         lastdi = -1
         nplv = nplv/2
         go to 230
c  addition in the u-direction.
 215    nplusu = nplu
        istart = 0
        if(iopt(2).eq.0) istart = 1
        do 220 l=1,nplusu
c  add a new knot in the u-direction
          call nknot(u,mu,tu,nu,nmax,fpintu,nrdatu,nrintu,istart)
c  test whether we cannot further increase the number of knots in the
c  u-direction.
          if(nu.eq.nue) go to 270
 220    continue
        go to 270
 230    if(nv.eq.nve) go to 210
         if(lastdi.le.0) lastdi = 0
         lastdi = lastdi+1
         if(lastdi.lt.5) go to 235
         lastdi = 1
         nplu = nplu/2
         go to 210
c  addition in the v-direction.
 235    nplusv = nplv
        do 240 l=1,nplusv
c  add a new knot in the v-direction.
          call nknot(v,mv,tv,nv,nmax,fpintv,nrdatv,nrintv,1)
c  test whether we cannot further increase the number of knots in the
c  v-direction.
          if(nv.eq.nve) go to 270
 240    continue
c  restart the computations with the new set of knots.
 270  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
c  test whether the least-squares spline is a solution of our
c  approximation problem.
      call sqgrid_c(r,nu,nv,c,fp,fpintu,fpintv)
      if(ider(1).eq.0) fp = fp+(r0-dr(1))**2
      if(ider(3).eq.0) fp = fp+(r1-dr(4))**2
      fpms = fp-s
      if(abs(fpms).lt.acc) go to 440
      if(fpms.gt.0.) go to 410
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(u,v)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
c  this smoothing spline depends on the parameter p in such a way that c
c    f(p) = sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. then iteratively we have to determine the     c
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
      do 305 i=1,6
         drr(i) = dr(i)
 305  continue
      iflag = 0
      icheck = 0
c  iteration process to find the root of f(p)=s.
      do 320 iter = 1,maxit
c  find the smoothing spline sp(u,v).
        call optsph(u,mu,v,mv,r,muv,r0,r1,tu,nu,tv,nv,p,iopt,idd,iflag,
     2  step,drr,c,ncof,sq)
        iflag = 1
c  calculate the corresponding value f(p).
        call sqgrid_c(r,nu,nv,c,fp,fpintu,fpintv)
c  test whether the approximation sp(u,v) is an acceptable solution.
        if(ider(1).eq.0) fp = fp+(r0-drr(1))**2
        if(ider(3).eq.0) fp = fp+(r1-drr(4))**2
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(icheck.ne.0) go to 315
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p = p*0.1e-02
        p3 = p2
        f3 = f2
        go to 320
 310    if((f1-f2).gt.acc) go to 315
c  our initial choice of p is too small
        p = p*0.1e+04
        p1 = p2
        f1 = f2
        go to 320
c  test whether the iteration process proceeds as theoretically
c  expected.
 315    if(f2.ge.f1 .or. f2.le.f3) go to 410
        icheck = 1
c  find the new value of p.
        p = ration(p1,f1,p2,f2,p3,f3)
 320  continue
c  error codes.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
 440  return
      end
