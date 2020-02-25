c
      subroutine optsph(u,mu,v,mv,r,muv,r0,r1,tu,nu,tv,nv,p,iopt,ider,
     1 iflag,step,dr,c,ncof,sq)
c
c      implicit none
      implicit double precision (A-H, O-Z)
      implicit integer (i-n)
c
c  given the set of function values r(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, optsph determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) will be periodic in the variable v and will satisfy
c  the following constraints
c
c     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3)
c
c     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3)
c
c  and (if iopt(2) = 1)
c
c     d s(tu(1),v)
c     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  and (if iopt(3) = 1)
c
c     d s(tu(nu),v)
c     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  where the parameters dr(i) correspond to the derivative values at the
c  poles as defined in subroutine sphery.
c
c  the b-spline coefficients of sp(u,v) are determined as the least-
c  squares solution  of an overdetermined linear system which depends
c  on the value of p and on the values dr(i),i=1,...,6. the correspond-
c  ing sum of squared residuals sq is a simple quadratic function in
c  the variables dr(i). these may or may not be provided. the values
c  dr(i) which are not given will be determined so as to minimize the
c  resulting sum of squared residuals sq. in that case the user must
c  provide some initial guess dr(i) and some estimate (dr(i)-step,
c  dr(i)+step) of the range of possible values for these latter.
c
c  sp(u,v) also depends on the parameter p (p>0) in such a way that
c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
c      with given knots, satisfying the constraints.
c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
c      satisfying the constraints.
c    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)
c      is continuous and strictly decreasing for p>0.
c
c  calling sequence:
c     call optsph(u,mu,v,mv,r,r0,r1,tu,nu,tv,nv,p,iopt,ider,iflag,
c                 step,dr,c,sq)
c
c  input parameters:
c    u,mu,v,mv,r,r0,r1 : see sphery.
c    tu,nu,tv,nv       : see sphery.
c    p    : real value, giving the positive smoothing parameter p.
c           p<=0  is assumed to stand for p=infinity. so in that case
c           the least-squares spline approximation is determined.
c    iopt              : see sphery ; iopt(1) is not used here.
c    ider : integer array, specifying which derivatives dr(i) are given
c           and which are to be determined, i.e.
c            ider(1) = 1 : dr(1) is given.
c                    = 0 : dr(1) must be determined; data value r0 given
c                    =-1 : dr(1) must be determined; no value r0 given
c            ider(2) = 1 : dr(2) and dr(3) are given
c                    = 0 : dr(2) and dr(3) must be determined.
c            ider(3) = 1 : dr(4) is given.
c                    = 0 : dr(4) must be determined; data value r1 given
c                    =-1 : dr(4) must be determined; no value r1 given
c            ider(4) = 1 : dr(5) and dr(6) are given
c                    = 0 : dr(5) and dr(6) must be determined.
c    iflag: integer flag.
c      iflag < 0 : first call of the routine with a new set of data.
c      iflag = 0 : second call of the routine in case p<=0 at the first
c                  call.
c      iflag > 0 : second call of the routine in case p>0 at the first
c                  call. also all next calls.
c      if iflag >= 0 the user must provide with a common block values.
c                    (see grisph)
c    step : real array, indicating whether all derivative values dr(i),
c           i=1,2,3 are given (step(1)<=0) or not (step(1)>0) . if
c           step(1)>0 it must give an estimate of the range of possible
c           values for the unknown dr(i),i=1,2,3.; similar comments
c           for step(2) with respect to the values dr(i),i=4,5,6.
c    dr   : real array, containing either the given value for the
c           derivative dr(i) or an estimate of its optimal value.
c
c  output parameters:
c    dr   : real array, containing the optimal values of the derivatives
c           dr(i) ; the given values dr(i) remain unchanged.
c    c    : array,length (nu-4)*(nv-4), which contains the
c           b-spline coefficients of sp(u,v).
c    sq   : real value, giving the sum of squared residuals f(p).
c
c  restrictions:
c    1),2),3),4)and 5) : see sphery
c    6) restrictions on the position of the knots.
c         tu(i) < tu(i+1) ,i=4,...,nu-4.
c         tu(i) = tu(4) = 0 ; tu(nu-i) = tu(nu) = pi , i=1,2,3.
c         tv(j) < tv(j+1) ,j=4,...,nv-4.
c         tv(nv+j-7)-tv(j) = tv(nv-3)-tv(4), j=1,...,7.
c       there must be at least one subset of distinct u(i)-values,i.e.
c       uu(j), and at least one subset of distinct v(i)-values,i.e.
c       vv(j), such that
c         tu(j+ib) < uu(j) < tu(j+4+ib), j=1,2,...,nu-4-ib-ie
c         tv(j+jb) < vv(j) < tv(j+4+jb), j=1,2,...,nv-7
c       where jb= either 0,1,2 or 3, ib=iopt(2)+1, ie=iopt(3)+1
c
c  other subroutines required:
c    grisph,symsys,cytri1,cytri2,cossin,rotate and bsplin.
      dimension u(mu),v(mv),r(muv),dr(6),tu(nu),tv(nv),c(ncof),
     1 step(2),delta(6),drr(6),sum(6),nr(6),a(6,6),g(6),iopt(3),
     2 ider(4)

       external grisph_c
c
c  we calculate the smoothing spline sp(u,v) according to the input
c  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0 = ider(1)
      id1 = ider(3)
      call grisph_c(r,dr,iop0,iop1,tu,nu,tv,nv,p,iflag,0,c,sq)
      sq0 = 0.
      sq1 = 0.
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
c in case all derivative values dr(i) are given (step<=0) or in case
c we have spline interpolation, we accept this spline as a solution.
      if(sq.le.0.) return
      if(step(1).le.0. .and. step(2).le.0.) return
      do 10 i=1,6
        drr(i) = dr(i)
  10  continue
c number denotes the number of derivative values dr(i) that still must
c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 20
      number = 1
      nr(1) = 1
      delta(1) = step(1)
  20  if(iop0.eq.0) go to 30
      if(ider(2).ne.0) go to 30
      step2 = step(1)*3./(tu(5)-tu(4))
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  30  if(id1.gt.0) go to 40
      number = number+1
      nr(number) = 4
      delta(number) = step(2)
  40  if(iop1.eq.0) go to 50
      if(ider(4).ne.0) go to 50
      step2 = step(2)*3./(tu(nu)-tu(nu-4))
      nr(number+1) = 5
      nr(number+2) = 6
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  50  if(number.eq.0) return
c the sum of squared residulas sq is a quadratic polynomial in the
c parameters g(j). we determine the unknown coefficients of this
c polymomial by calculating (number+1)*(number+2)/2 different splines
c according to specific values for g(j).
      do 60 i=1,number
         l = nr(i)
         step1 = delta(i)
         drr(l) = dr(l)+step1
         call grisph_c(r,drr,iop0,iop1,tu,nu,tv,nv,p,1,1,c,sum(i))
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1
         call grisph_c(r,drr,iop0,iop1,tu,nu,tv,nv,p,1,1,c,sqq)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sqq = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-2.*sq)/step1**2
         if(a(i,i).le.0.) go to 110
         g(i) = (sqq-sum(i))/(step1+step1)
  60  continue
      if(number.eq.1) go to 90
      do 80 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         drr(l1) = dr(l1)+step1
         i1 = i-1
         do 70 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            drr(l2) = dr(l2)+step2
            call grisph_c(r,drr,iop0,iop1,tu,nu,tv,nv,p,1,1,
     2      c,sqq)
            if(id0.eq.0) sq0 = (r0-drr(1))**2
            if(id1.eq.0) sq1 = (r1-drr(4))**2
            sqq = sqq+sq0+sq1
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            drr(l2) = dr(l2)
  70     continue
         drr(l1) = dr(l1)
  80  continue
c the optimal values g(j) are found as the solution of the system
c d (sq) / d (g(j)) = 0 , j=1,...,number.
  90  call symsys(a,number,g)
      do 100 i=1,number
         l = nr(i)
         dr(l) = dr(l)+g(i)
 100  continue
c we determine the spline sp(u,v) according to the optimal values g(j).
 110  call grisph_c(r,dr,iop0,iop1,tu,nu,tv,nv,p,iflag,
     1 0,c,sq)
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end
