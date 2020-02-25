c
      subroutine grisph(u,mu,v,mv,r,muv,dr,iop0,iop1,tu,nu,tv,nv,
     2 p,iflag,iback,c,ncf,sq,right,mnumv,q,mvnu,au,av1,av2,a0,a1,
     3 b0,b1,c0,c1,d,dd,nmax,spu,spv,bu,bv,cosi,nru,nrv)
c      implicit none
      implicit integer (I-N)
      implicit double precision (A-H, O-Z)
c  given the set of function values r(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, grisph determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) has the same properties as the one considered in sub-
c  routine optsph, the only difference being that all values dr(i) must
c  be provided here.
c
c  calling sequence:
c     call grisph(u,mu,v,mv,r,dr,iop0,iop1,tu,nu,tv,nv,p,iflag,iback,c,
c                 sq)
c
c  input parameters:
c    u,mu,v,mv,r : see sphery
c    dr          : derivative values dr(i) (see optsph)
c    iop0,iop1   : the parameters iopt(2) and iopt(3) of sphery.
c    tu,nu,tv,nv : see sphery
c    p,iflag     : see optsph
c    iback       : integer flag.
c      iback = 0 : the routine will return the b-spline coefficients.
c      iback > 0 : the routine will only return the sum of squared
c                  residuals , not the b-spline coefficients.
c
c  output parameters:
c    c    : array,length (nu-4)*(nv-4), which contains the
c           b-spline coefficients of sp(u,v). (if iback=0)
c    sq   : real value, giving the sum of squared residuals .
c
c  restrictions: see optsph
c
c  other subroutines required:
c    cytri1,cytri2,cossin, rotate and bsplin.
c
c  the following auxiliary arrays must have dimension specifications
c  at least  right(max(nu,mvv)),q(mvv*nu),au(4,nu),av1(6,nv),av2(4,nv),
c  d(nu),dd(nv),cc(nv),a0,a1,b0,b1(2,mv), with mvv=mv+nv.
      dimension u(mu),v(mv),r(muv),tu(nu),tv(nv),c(ncf),h(8) 
     2 ,d(nu),right(mnumv),q(mvnu),au(4,nu),av1(6,nv),av2(4,nv),h1(5)
     3 ,h2(4),c0(nv),c1(nv),a0(2,mv),a1(2,mv),b0(2,mv),b1(2,mv)
     4 ,dd(nv),dr(6), spu(mu,4),nru(mu),spv(mv,4),nrv(mv),bu(nmax,5),
     5  bv(nmax,5), cosi(2,nmax)

      external bsplin_c


c  common/values/spu(mu,4),nru(mu),spv(mv,4),nrv(mv),bu(nu,5),bv(nv,5),
c   cosi(2,nv)
c    spu  : array,which contains the value of the non-zero b-splines
c           at each u-value u(i).
c    nru  : integer array, which indicates which b-splines are non-zero
c           at each u-value u(i).
c    spv  : array, which contains the value of the non-zero b-splines
c           at each v-value v(j).
c    nrv  : integer array, which indicates which b-splines are non-zero
c           at each v-value v(j).
c    bu   : array, which contains the discontinuity jumps of the
c           derivatives (order 3) of the b-splines at the interior knots
c           tu(5),...,tu(nu-4).
c    bv   : array,which contains the discontinuity jumps of the
c           derivatives (order 3) of the b-splines at the interior knots
c           tv(5),...,tv(nv-4).
c    cosi : array, which contains the b-spline coefficients of the
c           periodic cubic spline interpolants for cos(v) and sin(v).

c  let
c               !     (spu)      !            !     (spv)      !
c        (au) = ! -------------- !     (av) = ! -------------- !
c               ! sqrt(1/p) (bu) !            ! sqrt(1/p) (bv) !
c
c                                ! r  ' 0 !
c                            q = ! ------ !
c                                ! 0  ' 0 !
c
c  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
c                coefficients.
c       r      : the mu x mv matrix which contains the function values.
c       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
c                according to the least-squares problems in the u-,resp.
c                v-direction.
c       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
c                containing the discontinuity jumps of the derivatives
c                of the b-splines in the u-,resp.v-variable at the knots
c  the b-spline coefficients of the smoothing spline are then calculated
c  as the least-squares solution of the following over-determined linear
c  system of equations
c
c  (1)  (av) c (au)' = q
c
c  subject to the constraints
c
c  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
c
c  (3)  if iop0 = 0  c(1,j) = dr(1)
c          iop0 = 1  c(1,j) = dr(1)
c                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))*
c                            tu(5)/3. = c0(j) , j=1,2,...nv-4
c
c  (4)  if iop1 = 0  c(nu-4,j) = dr(4)
c          iop1 = 1  c(nu-4,j) = dr(4)
c                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j))
c                                *(tu(nu-4)-tu(nu-3))/3. = c1(j)
c
c  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-2
      if(p.gt.0.) pinv = 1./sqrt(p)
c  it depends on the value of iflag and p whether the matrices (spu),
c  (spv), (cosi), (bu) and (bv) still must be determined.
      if(iflag) 10,90,150
c  calculate the non-zero elements of the matrix (spu) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the u-direction.
  10  l = 4
      l1 = 5
      number = 0
      do 30 it=1,mu
        arg = u(it)
  15    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 15
  20    call bsplin_c(tu,nu,3,arg,l,h)
        do 25 i=1,4
          spu(it,i) = h(i)
  25    continue
        nru(it) = number
  30  continue
c  calculate the non-zero elements of the matrix (spv) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the v-direction.
      l = 4
      l1 = 5
      number = 0
      do 50 it=1,mv
        arg = v(it)
  35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40
        l = l1
        l1 = l+1
        number = number+1
        go to 35
  40    call bsplin_c(tv,nv,3,arg,l,h)
        do 45 i=1,4
          spv(it,i) = h(i)
  45    continue
        nrv(it) = number
  50  continue
      if(iop0.eq.0 .and. iop1.eq.0) go to 85
c  calculate the coefficients of the interpolating splines for cos(v)
c  and sin(v).
      do 55 i=1,nv4
         cosi(1,i) = 0.
         cosi(2,i) = 0.
  55  continue
      if(nv7.lt.4) go to 85
      do 65 i=1,nv7
         l = i+3
         arg = tv(l)
         call bsplin_c(tv,nv,3,arg,l,h)
         do 60 j=1,3
            av1(j,i) = h(j)
  60     continue
         cosi(1,i) = cos(arg)
         cosi(2,i) = sin(arg)
  65  continue
      call cytri1(av1,nv7)
      do 80 j=1,2
         do 70 i=1,nv7
            right(i) = cosi(j,i)
  70     continue
         call cytri2(av1,nv7,right,right)
         do 75 i=1,nv7
            cosi(j,i+1) = right(i)
  75     continue
         cosi(j,1) = cosi(j,nv7+1)
         cosi(j,nv7+2) = cosi(j,2)
         cosi(j,nv4) = cosi(j,3)
  80  continue
  85  if(p.le.0.) go to  150
c  calculate the non-zero elements of the matrix (bu).
  90  if(nu8.eq.0) go to 120
      fac = dble(nu7)/(tu(nu-3)-tu(4))
      do 110 i=1,nu8
         l = i+4
         k0 = i
         k1 = l
         do 95 j=1,4
           h(j) = tu(l)-tu(k0)
           k0 = k0+1
           k1 = k1+1
           h(j+4) = tu(l)-tu(k1)
  95     continue
         lp = i
         do 105 j=1,5
           prod = h(j)
           k0 = j
           do 100 k=1,3
             k0 = k0+1
             prod = prod*h(k0)*fac
 100       continue
           bu(i,j) = (tu(lp+4)-tu(lp))/prod
           lp = lp+1
 105     continue
 110  continue
c  calculate the non-zero elements of the matrix (bv).
 120  if(nv8.eq.0) go to 150
      fac = dble(nv7)/(tv(nv-3)-tv(4))
      do 140 i=1,nv8
         l = i+4
         k0 = i
         k1 = l
         do 125 j=1,4
           h(j) = tv(l)-tv(k0)
           k0 = k0+1
           k1 = k1+1
           h(j+4) = tv(l)-tv(k1)
 125     continue
         lp = i
         do 135 j=1,5
           prod = h(j)
           k0 = j
           do 130 k=1,3
             k0 = k0+1
             prod = prod*h(k0)*fac
 130       continue
           bv(i,j) = (tv(lp+4)-tv(lp))/prod
           lp = lp+1
 135     continue
 140  continue
c  substituting (2),(3) and (4) into (1), we obtain the overdetermined
c  system
c         (5)  (avv) (cc) (auu)' = (qq)
c  from which the nuu*nv7 remaining coefficients
c         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
c  the elements of (cc), are then determined in the least-squares sense.
c  simultaneously, we compute the resulting sum of squared residuals sq.
 150  dr01 = dr(1)
      dr11 = dr(4)
      do 155 i=1,mv
         a0(1,i) = dr01
         a1(1,i) = dr11
 155  continue
      if(nv8.eq.0 .or. p.le.0.) go to 165
      do 160 i=1,nv8
         b0(1,i) = 0.
         b1(1,i) = 0.
 160  continue
 165  mvv = mv
      if(iop0.eq.0) go to 195
      fac = (tu(5)-tu(4))/3.
      dr02 = dr(2)*fac
      dr03 = dr(3)*fac
      do 170 i=1,nv4
         c0(i) = dr01+dr02*cosi(1,i)+dr03*cosi(2,i)
 170  continue
      do 180 i=1,mv
         number = nrv(i)
         fac = 0.
         do 175 j=1,4
            number = number+1
            fac = fac+c0(number)*spv(i,j)
 175     continue
         a0(2,i) = fac
 180  continue
      if(nv8.eq.0 .or. p.le.0.) go to 195
      do 190 i=1,nv8
         number = i
         fac = 0.
         do 185 j=1,5
            fac = fac+c0(number)*bv(i,j)
            number = number+1
 185     continue
         b0(2,i) = fac*pinv
 190  continue
      mvv = mv+nv8
 195  if(iop1.eq.0) go to 225
      fac = (tu(nu4)-tu(nu4+1))/3.
      dr12 = dr(5)*fac
      dr13 = dr(6)*fac
      do 200 i=1,nv4
         c1(i) = dr11+dr12*cosi(1,i)+dr13*cosi(2,i)
 200  continue
      do 210 i=1,mv
         number = nrv(i)
         fac = 0.
         do 205 j=1,4
            number = number+1
            fac = fac+c1(number)*spv(i,j)
 205     continue
         a1(2,i) = fac
 210  continue
      if(nv8.eq.0 .or. p.le.0.) go to 225
      do 220 i=1,nv8
         number = i
         fac = 0.
         do 215 j=1,5
            fac = fac+c1(number)*bv(i,j)
            number = number+1
 215     continue
         b1(2,i) = fac*pinv
 220  continue
      mvv = mv+nv8
c  we first determine the matrices (auu) and (qq). then we reduce the
c  matrix (auu) to an unit upper triangular form (ru) using givens
c  rotations without square roots. we apply the same transformations to
c  the rows of matrix qq to obtain the mv x nuu matrix g.
c  we store matrix (ru) into au and g into q.
 225  l = mvv*nuu
c  initialization.
      sq = 0.
      if(l.eq.0) go to 245
      do 230 i=1,l
        q(i) = 0.
 230  continue
      do 240 i=1,nuu
        d(i) = 0.
        do 240 j=1,4
          au(j,i) = 0.
 240  continue
      l = 0
 245  nrold = 0
      n1 = nrold+1
      do 420 it=1,mu
        number = nru(it)
c  find the appropriate column of q.
 250    do 260 j=1,mvv
           right(j) = 0.
 260    continue
        if(nrold.eq.number) go to 280
        if(p.le.0.) go to 410
c  fetch a new row of matrix (bu).
        do 270 j=1,5
          h(j) = bu(n1,j)*pinv
 270    continue
        i0 = 1
        i1 = 5
        go to 310
c  fetch a new row of matrix (spu).
 280    do 290 j=1,4
          h(j) = spu(it,j)
 290    continue
c  find the appropriate column of q.
        do 300 j=1,mv
          l = l+1
          right(j) = r(l)
 300    continue
        i0 = 1
        i1 = 4
 310    j0 = n1
        j1 = nu7-number
c  take into account that we eliminate the constraints (3)
 315     if(j0-1.gt.iop0) go to 335
         fac0 = h(i0)
         do 320 j=1,mv
            right(j) = right(j)-fac0*a0(j0,j)
 320     continue
         if(mv.eq.mvv) go to 330
         j = mv
         do 325 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac0*b0(j0,jj)
 325     continue
 330     j0 = j0+1
         i0 = i0+1
         go to 315
c  take into account that we eliminate the constraints (4)
 335     if(j1-1.gt.iop1) go to 360
         fac1 = h(i1)
         do 340 j=1,mv
            right(j) = right(j)-fac1*a1(j1,j)
 340     continue
         if(mv.eq.mvv) go to 350
         j = mv
         do 345 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac1*b1(j1,jj)
 345     continue
 350     j1 = j1+1
         i1 = i1-1
         go to 335
 360     irot = nrold-iop0-1
         if(irot.lt.0) irot = 0
c  rotate the new row of matrix (auu) into triangle.
        wi = 1.
        if(i0.gt.i1) go to 390
        do 385 i=i0,i1
          if(wi.eq.0.) go to 400
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 385
c  calculate the parameters of the givens transformation.
          call cossin(piv,wi,d(irot),co,si)
c  apply that transformation to the rows of matrix (qq).
          iq = (irot-1)*mvv
          do 370 j=1,mvv
            iq = iq+1
            call rotate(piv,co,si,right(j),q(iq))
 370      continue
c  apply that transformation to the columns of (auu).
          if(i.eq.i1) go to 385
          i2 = 0
          i3 = i+1
          do 380 j=i3,i1
            i2 = i2+1
            call rotate(piv,co,si,h(j),au(i2,irot))
 380      continue
 385    continue
c  we update the sum of squared residuals.
 390    do 395 j=1,mvv
          sq = sq+wi*right(j)**2
 395    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = n1
        n1 = n1+1
        go to 250
 420  continue
      if(nuu.eq.0) go to 890
c  we determine the matrix (avv) and then we reduce her to an unit
c  upper triangular form (rv) using givens rotations without square
c  roots. we apply the same transformations to the columns of matrix
c  g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h.
c  we store matrix (rv) into av1 and av2, h into c.
c  the nv7 x nv7 triangular unit upper matrix (rv) has the form
c              ! av1 '     !
c       (rv) = !     ' av2 !
c              !  0  '     !
c  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper
c  triangular matrix of bandwidth 5.
      ncof = nuu*nv7
c  initialization.
      do 430 i=1,ncof
        c(i) = 0.
 430  continue
      do 440 i=1,nv4
        dd(i) = 0.
        do 440 j=1,4
          av1(j,i) = 0.
          av2(j,i) = 0.
 440  continue
      jper = 0
      nrold = 0
      do 770 it=1,mv
        number = nrv(it)
 450    if(nrold.eq.number) go to 480
        if(p.le.0.) go to 760
c  fetch a new row of matrix (bv).
        n1 = nrold+1
        do 460 j=1,5
          h(j) = bv(n1,j)*pinv
 460    continue
c  find the appropiate row of g.
        do 465 j=1,nuu
          right(j) = 0.
 465    continue
        if(mv.eq.mvv) go to 510
        l = mv+n1
        do 470 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 470    continue
        go to 510
c  fetch a new row of matrix (spv)
 480    h(5) = 0.
        do 490 j=1,4
          h(j) = spv(it,j)
 490    continue
c  find the appropiate row of g.
        l = it
        do 500 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 500    continue
 510    wi = 1.
c  test whether there are non-zero values in the new row of (avv)
c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
         if(nrold.lt.nv11) go to 710
         if(jper.ne.0) go to 550
c  initialize the matrix (av2).
         jk = nv11
         do 540 i=1,4
            ik = jk
            do 520 j=1,4
               if(ik.le.0) go to 530
               av2(i,ik) = av1(j,ik)
               ik = ik-1
 520        continue
 530        jk = jk+1
 540     continue
         jper = 1
c  if one of the non-zero elements of the new row corresponds to one of
c  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
c  (2) for setting up this row of (avv). the row is stored in h1( the
c  part with respect to av1) and h2 (the part with respect to av2).
 550     do 560 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 560     continue
         h1(5) = 0.
         j = nrold-nv11
         do 600 i=1,5
            j = j+1
            l0 = j
 570        l1 = l0-4
            if(l1.le.0) go to 590
            if(l1.le.nv11) go to 580
            l0 = l1-nv11
            go to 570
 580        h1(l1) = h(i)
            go to 600
 590        h2(l0) = h2(l0) + h(i)
 600     continue
c  rotate the new row of (avv) into triangle.
         if(nv11.le.0) go to 670
c  rotations with the rows 1,2,...,nv11 of (avv).
         do 660 j=1,nv11
            piv = h1(1)
            i2 = min0(nv11-j,4)
            if(piv.eq.0.) go to 640
c  calculate the parameters of the givens transformation.
            call cossin(piv,wi,dd(j),co,si)
c  apply that transformation to the columns of matrix g.
            ic = j
            do 610 i=1,nuu
               call rotate(piv,co,si,right(i),c(ic))
               ic = ic+nv7
 610        continue
c  apply that transformation to the rows of (avv) with respect to av2.
            do 620 i=1,4
               call rotate(piv,co,si,h2(i),av2(i,j))
 620        continue
c  apply that transformation to the rows of (avv) with respect to av1.
            if(i2.eq.0) go to 670
            do 630 i=1,i2
               call rotate(piv,co,si,h1(i+1),av1(i,j))
 630        continue
 640        do 650 i=1,i2
               h1(i) = h1(i+1)
 650        continue
            h1(i2+1) = 0.
 660     continue
c  rotations with the rows nv11+1,...,nv7 of avv.
 670     do 700 j=1,4
            if(wi.eq.0.) go to 750
            ij = nv11+j
            if(ij.le.0) go to 700
            piv = h2(j)
            if(piv.eq.0.) go to 700
c  calculate the parameters of the givens transformation.
            call cossin(piv,wi,dd(ij),co,si)
c  apply that transformation to the columns of matrix g.
            ic = ij
            do 680 i=1,nuu
               call rotate(piv,co,si,right(i),c(ic))
               ic = ic+nv7
 680        continue
            if(j.eq.4) go to 700
c  apply that transformation to the rows of (avv) with respect to av2.
            j1 = j+1
            do 690 i=j1,4
               call rotate(piv,co,si,h2(i),av2(i,ij))
 690        continue
 700     continue
c  we update the sum of squared residuals.
         do 705 i=1,nuu
           sq = sq+wi*d(i)*right(i)**2
 705     continue
         go to 750
c  rotation into triangle of the new row of (avv), in case the elements
c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
 710     irot =nrold
         do 740 i=1,5
            if(wi.eq.0.) go to 750
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 740
c  calculate the parameters of the givens transformation.
            call cossin(piv,wi,dd(irot),co,si)
c  apply that transformation to the columns of matrix g.
            ic = irot
            do 720 j=1,nuu
               call rotate(piv,co,si,right(j),c(ic))
               ic = ic+nv7
 720        continue
c  apply that transformation to the rows of (avv).
            if(i.eq.5) go to 740
            i2 = 0
            i3 = i+1
            do 730 j=i3,5
               i2 = i2+1
               call rotate(piv,co,si,h(j),av1(i2,irot))
 730        continue
 740     continue
c  we update the sum of squared residuals.
         do 745 i=1,nuu
           sq = sq+wi*d(i)*right(i)**2
 745     continue
 750     if(nrold.eq.number) go to 770
 760     nrold = nrold+1
         go to 450
 770  continue
c  test whether the b-spline coefficients must be determined.
      if(iback.ne.0) return
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (rv) (cc) (ru)' = h.
c  first step: solve the system  (rv) (c1) = h.
      if(nv7.eq.1) go to 850
      ic = -1
      do 840 k=1,nuu
         l = nv7-1
         ic = ic+nv7
         iq = ic
         do 790 i=2,4
            store = c(iq)
            j = 6-i
            l0 = iq
            do 780 l1 = j,4
               l0 = l0+1
               store = store-c(l0)*av2(l1,l)
 780        continue
            c(iq) = store
            iq = iq-1
            l = l-1
            if(l.eq.0) go to 840
 790     continue
         do 810 i=1,nv11
            store = c(iq)
            l0 = ic+1
            l1 = 4
            do 800 j=1,4
               store = store-c(l0)*av2(l1,l)
               l0 = l0-1
               l1 = l1-1
 800        continue
            c(iq) = store
            iq = iq-1
            l = l-1
 810     continue
         if(nv11.eq.1) go to 840
         iq = ic-4
         l = nv11-1
         do 830 j=2,nv11
            store = c(iq)
            i1 = min0(4,j-1)
            l0 = iq
            do 820 l1=1,i1
               l0 = l0+1
               store = store-c(l0)*av1(l1,l)
 820        continue
            c(iq) = store
            iq = iq-1
            l = l-1
 830     continue
 840  continue
 850  if(nuu.le.1) go to 890
c  second step: solve the system  (cc) (ru)' = (c1).
      k = ncof-2*nv7
      do 880 j=1,nv7
        k = k+1
        k1 = k
        i = nuu-1
        do 870 i1=2,nuu
          store = c(k1)
          l2 = min0(4,i1-1)
          l1 = k1
          do 860 l3=1,l2
            l1 = l1+nv7
            store = store-c(l1)*au(l3,i)
 860      continue
          c(k1) = store
          i = i-1
          k1 = k1-nv7
 870    continue
 880  continue
c  calculate from the conditions (2)-(3)-(4), the remaining b-spline
c  coefficients.
 890  ncof = nu4*nv4
      j = ncof
      do 900 l=1,nv4
         q(l) = dr01
         q(j) = dr11
         j = j-1
 900  continue
      i = nv4
      j = 0
      if(iop0.eq.0) go to 920
      do 910 l=1,nv4
         i = i+1
         q(i) = c0(l)
 910  continue
 920  if(nuu.eq.0) go to 960
      do 950 l=1,nuu
         ii = i
         do 930 k=1,nv7
            i = i+1
            j = j+1
            q(i) = c(j)
 930     continue
         do 940 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 940     continue
 950  continue
 960  if(iop1.eq.0) go to 980
      do 970 l=1,nv4
         i = i+1
         q(i) = c1(l)
 970  continue
 980  do 990 i=1,ncof
         c(i) = q(i)
 990  continue
      return
      end

