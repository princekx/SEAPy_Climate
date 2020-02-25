
      PROGRAM velvorc

C     NB DX=1.875,dy=1.25 degrees
C     FOR THIS PROGRAM TO WORK IT REQUIRES DRS LIBARIES.
C     THIS WILL COMPILE ON THE HADLEY CENTRE HP'S USING
C     f77 +U77 velvorc2.f ~hadls/umdiagn/lib/libdrs.a
c     ANYWHERE ELSE YOU WILL NEDD TO GET HOLD OF THESE LIBARIES AND
C     CHANGE THE PROGRAM ACCORDINGLY
      INCLUDE '/users/nutis/kih/local/include/drsdef.H'
C     A = RADIUS OF EARTH, DPHI,DLAN AREA MERIDIONAL & LONGITUDINAL GRID
C     SPACINGS
      REAL a,DPHI,PI,DLAN
C     NX,NY = NUMBER OF GRID POINTS
      INTEGER NFILE,NX,NY,TL,t,I3
      PARAMETER(NFILE=3,NX=96,NY=72)
      INTEGER INFILE(NFILE)
      PARAMETER(a=6.371E6,DPHI=((2.5*3.141592)/180),PI=3.141592)
      PARAMETER(DLAN=((PI*3.75)/180))
      REAL UARR(NX,NY),VARR(NX,NY),MEANU(NX,NY),DVDLN(NX,NY)
      REAL MNPHI(NY),DUDPHI(NX,NY),VORT(NX,NY),VPHI(NY),MEANV(NX,NY)
      REAL R2
      CHARACTER*23 DRSFILE(NFILE)
c     EDIT TL TO THE NUMBER OF TIME STEPS IN YOUR DATA
      TL=736
c     INPUT FILES FOR U AND V OUTPUT FILE NAME/ROUTE
      DRSFILE(1)='UM_amip2_U_mjjaso95.dic'
      DRSFILE(2)='UM_amip2_V_mjjaso95.dic'
      DRSFILE(3)='UM_amip2_VOR_mjjaso95'
c      DRSFILE(1)='/home/hc0200/haddn/ut678_12_01'
c      DRSFILE(2)='/home/hc0200/haddn/vt678_12_01'
c      DRSFILE(3)='/home/hc0200/haddn/vorcfile'
c     I3 IS THE TIME OF THE FIRST TIME STEP LESS 0.25 EG FOR THIS YEAR DATA
C     STARTS AT 1371.25
      I3 = 143886
C     READ IN U DATA into UARR
      DO t=1,TL
      INFILE(1)=10
      R2=(t-1) * 6.0 + I3 
      ierr=aslun(INFILE(1),DRSFILE(1),INFILE(1)+1,' ',1)
      ierr=cluvdb()
      ierr = setname(' ','ua',' ',' ',' ')
      ierr = setdim(1,'longitude',' ',NX,1.875,358.125)
      ierr = setdim(2,'latitude',' ',NY,88.75,-88.75)
      ierr = setdim(3,'time',' ',1,R2,R2)
      ierr = getdat(INFILE(1),UARR,NX*NY*4)
      ierr = cllun(INFILE(1))
C     READ IN V DATA into VARR
      INFILE(2)=20
      ierr=aslun(INFILE(2),DRSFILE(2),INFILE(2)+1,' ',1)
      ierr=cluvdb ()
      ierr = setname(' ','va',' ',' ',' ')
      ierr = setdim(1,'longitude',' ',NX,1.875,358.125)
      ierr = setdim(2,'latitude',' ',NY,88.75,-88.75)
      ierr = setdim(3,'time',' ',1,R2,R2)
      ierr = getdat(INFILE(2),VARR,NX*NY*4)
      ierr = cllun (INFILE(2))
C     TAKE AN AVERAGE OF PAIRS of x adjacent u values
      DO J=1,NY
      DO I=1,NX-1
        MEANU(I,J)=(UARR(I,J)+UARR(I+1,J))/2
      ENDDO
      ENDDO
      DO J=1,NY
        MEANU(NX,J)=(UARR(NX,J)+UARR(1,J))/2
      ENDDO
C      TAKE AN AVERAGE OF PAIRS of y adjacent v values  
      DO I=1,NX
      DO J=1,NY-1
        MEANV(I,J)=(VARR(I,J+1)+VARR(I,J))/2
      ENDDO
      ENDDO
      DO I=1,NX
        MEANV(I,NY)=(VARR(I,NY)+VARR(I,1))/2
      ENDDO
C     Now find dv/dlanda
      DO J=1,NY
      DO I=1,NX-1
        DVDLN(I,J)=(MEANV(I+1,J)-MEANV(I,J))/DLAN
      ENDDO
      ENDDO
      DO J=1,NY
        DVDLN(NX,J)=(MEANV(1,J)-MEANV(NX,J))/DLAN
      ENDDO
C     NOW d(ucos(phi)/dphi
C     THIS BIT IS VERY GRID SPECIFIC. ITS WORKING OUT THE LATTITUDE OF
C      THE MIDDLE POINT OF TWO ADJACENT U VALUES
      DO J=1,NY
        MNPHI(J)=(89.375-((J-1)*1.25))*((2*PI)/360)
      ENDDO
      DO I=1,NX
      DO J=1,NY-1
C       ******put in rads and chk pi!!!!!!!!*****
        DUDPHI(I,J)=((MEANU(I,J)*COS(MNPHI(J)))
     :-(MEANU(I,J+1)*COS(MNPHI(J+1))))/DPHI
     
      ENDDO
      ENDDO
      DO I=1,NX
        DUDPHI(I,NY)=((MEANU(I,NY)*COS(MNPHI(NY)))
     : -(MEANU(I,1)*COS(MNPHI(1))))/DPHI
      ENDDO
C     NOW find vorticities
      DO J=1,NY
      VPHI(J)=1.559888-(DPHI/2)-((J-1)*DPHI)
      ENDDO
      DO I=1,NX
      DO J=1,NY-1
      VORT(I,J)=(DVDLN(I,J)/(a*COS(VPHI(J))))
     :-(DUDPHI(I,J)/(a*COS(VPHI(J))))
      ENDDO
      ENDDO
C     Now save data in DRS file
      print*,DVDLN(1,1),DUDPHI(1,1),a*COS(VPHI(1)),(MNPHI(38))
       PRINT*,DPHI,MEANU(1,1),VORT(1,1)
      INFILE(3)=40
      IF(t.EQ.1)THEN
      ierr=aslun(INFILE(3),DRSFILE(3),INFILE(3)+1,' ',IDRS_CREATE)
      ELSE
      ierr=aslun(INFILE(3),DRSFILE(3),INFILE(3)+1,' ',IDRS_EXTEND)
      ENDIF
      ierr=cluvdb ()
      ierr = setname(' ','VORT','vorticity','s-1',' ')
      ierr = setdim(1,'longitude','degrees',NX,3.75,360.000)
      ierr = setdim(2,'lattitude','degrees',NY-1,87.5,-87.5)
      ierr = setdim(3,'pressure','mb',1,850,850)
      ierr = setdim(4,'time',' ',1,R2,R2)
      ierr = putdat(INFILE(3),VORT)
      ierr = cllun(INFILE(3))
       ENDDO
      stop
      END


















