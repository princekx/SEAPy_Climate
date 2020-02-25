
/* fortran preprocessor file for stats. plot using uniras */

#define     ABS_SWITCH     0

#ifndef     HAVEUNIRAS

      SUBROUTINE FILTERPL()

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                     C
C  Dummy routine for FILTERPL if no UNIRAS available.                 C
C                                                                     C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      INCLUDE 'paramc.ptr'

#ifdef    NAMELISTS

      READ(NAMCHN, NML=PARAMC)
      REWIND(20)

#endif

      WRITE(NUT, 900)

      RETURN
  900 FORMAT('/***ERROR***, no plotting possible without UNIRAS/')

      END

#else

      SUBROUTINE FILTERPL(XGRID, YGRID, Z, IX, IY, XP, YP, VF, DFIL, K,
     +                    CMI ,CMXG, CMYG, DCM, X1U, X2U, Y1U, Y2U, FR,
     +                    TH ,XFEAT, YFEAT, FENUM, X, Y, ZZ, IIAN, XAN,
     +                    YAN, XOR, YOR, ZOR, XBB, YBB, VBB, IBB, IBN,
     +                    IPBB)
      IMPLICIT NONE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                       C
C Program to plot the filtered data following pre-processing stage.     C
C May need modifying for ireggular data.                                C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INTEGER IX, IY, K, DFIL, DCM, CMI(DCM), IXD, IYD, ICDIM
      INTEGER ICM, I, J, X1U, X2U, Y1U, Y2U, FR, FENUM
      INTEGER IIAN, XAN(IIAN), YAN(IIAN)
      INTEGER XP(DFIL), YP(DFIL)
      INTEGER IPBB, IBN, IBB(IBN)
      REAL VF(DFIL), TH, PERIOD, CTOL, SCY
      REAL CMXG(DCM), CMYG(DCM), XFEAT(FENUM), YFEAT(FENUM)
      REAL XOR(FENUM), YOR(FENUM), ZOR(2*FENUM)
      REAL XBB(IBN), YBB(IBN), VBB(IBN)
      REAL PII

      PARAMETER(PII=57.295828)
      PARAMETER(ICDIM=100)

C
C------------------------------------------------------------------------
C
C Include the parameter file.
C
      INCLUDE 'paramf.ptr'
      INCLUDE 'paramc.ptr'
C
      REAL XGRID(IX), YGRID(IY), Z(IX, IY), FNUM
      REAL XMIN, XMAX, YMIN, YMAX, XSI, YSI, XM, 
     +     XOFF, YOFF, XSIZE, YSIZE, HEIGHT, ZL(2)
      REAL X(IX*IY), Y(IX*IY), ZZ(IX*IY), ZDUM(1)
      REAL YEXAG, MAXVEC, XSC, YSC
      INTEGER IC, TX, TY
      INTEGER SDARRY(ICDIM), ICNW
      LOGICAL FSTAT

      INTEGER IPID, IPD
      LOGICAL IPIDA

      INTEGER GETPID, SYSTEM
      CHARACTER NEWFIL*70, UDIR*70, CPID*8
      CHARACTER ANS*3

      DATA IPIDA / .FALSE. /

#ifdef    NAMELISTS

      READ(NAMCHN, NML=PARAMF)
      READ(NAMCHN, NML=PARAMC)
      REWIND(NAMCHN)

#endif

C
C---------------------------------------------------------------------
C
C Get process ID for unique output directory for graphics
C
      IF( .NOT. IPIDA) THEN
         WRITE(NUT, 792)
         READ(NIN, *) ANS
         IPIDA = .TRUE.
         IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
            IPID = GETPID()
            WRITE(NUT, 795) IPID
            WRITE(CPID,'(I8.8)') IPID
            UDIR = 'mkdir ' // PATH // '/' // CPID
            NEWFIL=PATH // '/' // CPID // '/' // FILNMF
            ipd=system(UDIR)

         ELSE
            NEWFIL = PATH // '/' // FILNMF
         ENDIF

      ENDIF

      print*,NEWFIL

      FSTAT= .TRUE. 
      ICM = 1
      FNUM=FR
      ZL(1)=ZINIT
      ZL(2)=STEP
      PERIOD = XGRID(IX) - XGRID(1)
      CTOL = PERIOD / 4.0

      YEXAG = 1.0
      MAXVEC = 0.001



C
C-----------------------------------------------------------------------
C
C Initialize plotting arrays
C

      IF(NU .EQ. 1) THEN

        DO 1 J=1, IY
         DO 2 I=1,IX
          IC = IX*(J-1) + I
          ZZ(IC) = 0.0
          X(IC) = XGRID(I)
          Y(IC) = YGRID(J) 
    2    CONTINUE
    1   CONTINUE 
        DO 3 I=1, DFIL
#if  ABS_SWITCH == 0
          ZZ(IX*(YP(I)+Y1U-2) + XP(I)+X1U-1)=SCALE*ABS(VF(I)) 
#else
          ZZ(IX*(YP(I)+Y1U-2) + XP(I)+X1U-1)=SCALE*VF(I)
#endif
    3   CONTINUE
      ELSE
        DO 4 J=1, IY
         DO 5 I=1, IX
          Z(I,J)=-1.0
    5    CONTINUE
    4   CONTINUE

        DO 6 J=Y1U, Y2U
         DO 7 I=X1U, X2U
           Z(I,J)=0.0
    7    CONTINUE
    6   CONTINUE
        DO 8 I=1,DFIL
#if  ABS_SWITCH == 0
          Z(XP(I)+X1U-1, YP(I)+Y1U-1) = SCALE*ABS(VF(I))
#else
          Z(XP(I)+X1U-1, YP(I)+Y1U-1) = SCALE*VF(I)
#endif
    8   CONTINUE
      ENDIF

C
C
C------------------------------------------------------------------------
C
C Set max. and min. values of the grid.
C
      XMIN=XGRID(1)
      XMAX=XGRID(IX)
      YMIN=YGRID(1)
      YMAX=YGRID(IY)   
           
C
C=========================================================================
C
C Plot data using UNIRAS.
C
C
C------------------------------------------------------------------------
C
C Check if file already exists for plotting.
C
      INQUIRE(FILE=NEWFIL, EXIST=FSTAT)
C
      CALL GROUTE(DEVICE)
      CALL GFILE(10, NOUT, NEWFIL)
      CALL GOPEN
      IF(FSTAT)  CALL GSEGWK(0)
      CALL GSEGCR(K)

      CALL GRPSIZ(XSI,YSI)
      XM=MIN(XSI,YSI)
      XSIZE=0.8*XM
      YSIZE=XSIZE
      XOFF=0.2*(XSI-XSIZE)
      YOFF=0.7*(YSI-YSIZE)

      IF(SCALEY .LE. 0.) THEN
         SCY = SCALEX * (YMAX - YMIN) / (XMAX - XMIN)
      ELSE
         SCY = SCALEY
      ENDIF

      CALL GLIMIT(XMIN, XMAX, YMIN, YMAX, 0.0, 0.0)
      CALL GVPORT(XOFF, YOFF, SCALEX*XSIZE, SCY*YSIZE)
      CALL RRECT(XOFF, YOFF, SCALEX*XSIZE,SCY*YSIZE,0 ,2)
      
      CALL GRIDEF(XGRID, IX, YGRID, IY, ZDUM, 1)

      IF(NU .EQ. 1) THEN
        CALL GSMTH(-2)
        CALL GINTP1(X, Y, ZZ, IX*IY, Z, IX, IY)
        DO 20 I=1, IY
         TY = (I-Y1U)*(Y2U-I)
         DO 30 J=1,IX
           TX = (J-X1U)*(X2U-J)
           IF((TY .LT. 0 .OR. TX .LT. 0) .AND. Z(J, I) .LE. 0.0 )
     +        Z(J, I) = -1.0
   30    CONTINUE
   20   CONTINUE
      ENDIF

      HEIGHT = 0.02*MIN(XSIZE,YSIZE)
      CALL GSCALE
      
      CWI(1) = CONLW
      CLSTY(1) = CONLS

C
C------------------------------------------------------------------------
C
C Plot processed data.
C
      CALL RUNDEF(VALUE,COL)
      CALL RCLASS(ZL,ICNUM,-5)
      
      IF(ITAB .GT. 6) CALL COLDEF(ITAB)
      
      CALL RSHADE(ITAB,0)

      IF(.NOT. IPBB .EQ. 2) THEN
         IF(ICONT .EQ. 0) THEN
            CALL GCONA(CHS*HEIGHT, 1, -2.0*XSIZE, 3)
            CALL GCONAB(9999, IFR)

C            CALL RFMODE('MIX')

            IF(INCT .GT. 0) THEN
               CALL RCLASS(ZNCT, INCT, 0)
               CALL RSHADE(ITAB, 0)
               CALL GCNR2S(Z, IX, IY)
            ENDIF

            CALL RCLASS(ZL,ICNUM,-5)
            CALL RSHADE(ITAB,0)

            CALL GCNR2V(Z, IX, IY)    
            CALL RQSHAD(SDARRY, ICNW)
            CALL GCONCO(SDARRY, ICNW)
            CALL GCNR2V(Z, IX, IY) 

         ELSE IF(ICONT .EQ. 1) THEN
	    CALL GCONSM(CSMTH)
	    CALL GCONDA(CLSTY, -1)
            CALL GCONWI(CWI, 1)
            CALL GCNR2S(Z, IX, IY)
         ENDIF
C
         IF(DOVLAY .EQ. 1) THEN
            CALL GWICOL(0.006*XSIZE, 2)
C         CALL GDOT(XP, YP, DFIL)
            DO 200 J=1, IY
              DO 300 I=1, IX
                IF(Z(I,J).GT.0.0) THEN
                  CALL GDOT(XGRID(I), YGRID(J), 1)
                ENDIF
  300         CONTINUE
  200       CONTINUE

         ELSE

            CALL GWICOL(0.006*XSIZE, 2)
            DO 400 J=1, IY
             DO 500 I=1, IX
               IYD = (J-Y1U)*(Y2U-J)
               IXD = (I-X1U)*(X2U-I)
               IF((IYD .LT. 0 .OR. IXD .LT. 0) .AND. Z(I,J).GT.0.0) 
     +         THEN
                 CALL GDOT(XGRID(I), YGRID(J), 1)
               ENDIF
  500        CONTINUE
  400       CONTINUE

         ENDIF

      ENDIF


      IF(IIAN .GT. 0) THEN

        CALL GWICOL(0.006*XSIZE, 32)

        DO 510 I=1, IIAN
          CALL GDOT(XGRID(XAN(I)), YGRID(YAN(I)), 1)

  510   CONTINUE

      ENDIF

      IF(NCSC .EQ. 0) THEN
         CALL GCOSCL(XMAX+0.1*(XMAX-XMIN), YMIN)
      ELSE IF(NCSC .GT. 0) THEN
         CALL GCOSCS(XMAX+0.01*(XMAX-XMIN), YMIN)
      ENDIF

C
C------------------------------------------------------------------------
C
C Plot country map.
C
      CALL GWICOL(0.001*XSIZE, 33)
   40 CONTINUE
      IF(ICM+1 .LE. DCM) THEN
       CALL GVECT(CMXG(ICM), CMYG(ICM), 0)
   50  CONTINUE   
       IF(ICM+1 .LE. DCM) THEN
         IF(ABS(CMI(ICM)) .EQ. ABS(CMI(ICM+1)) .AND.
     +      ABS(CMXG(ICM)-CMXG(ICM+1)) .LT. CTOL) THEN
            IF(CMI(ICM) .GT. 0 .AND. CMI(ICM+1) .GT. 0)
     +         CALL GVECT(CMXG(ICM+1), CMYG(ICM+1), 1)
            ICM = ICM+1
            GOTO 50
         ELSE
            ICM = ICM+1
            GOTO 40
         ENDIF
       ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C Plot feature points.
C

      IF(IIAN .GT. 0) THEN
         CALL GWICOL(0.008*XSIZE, 4)
      ELSEIF(IPBB .EQ. 2) THEN
         CALL GWICOL(0.008*XSIZE, 4)
      ELSE
         CALL GWICOL(0.006*XSIZE, 32)
      ENDIF
      CALL GDOT(XFEAT, YFEAT, FENUM)


C
C-----------------------------------------------------------------------
C
C Plot orientation vectors.
C

      CALL RUNDEF(RDEF, IDEF)

      IF(IIAN .GT. 0) THEN

         DO 60 I=1, FENUM
            ZOR(2*I-1) = SQRT(XOR(I) * XOR(I) + YOR(I) * YOR(I) 
     +                        * YEXAG * YEXAG)

            MAXVEC = MAX(MAXVEC, ZOR(2*I-1))

            ZOR(2*I) = ATAN2(YOR(I) * YEXAG, XOR(I)) * PII

            IF(ZOR(2*I-1) .LT. SUPVEC) THEN
              ZOR(2*I-1) = RDEF
            ENDIF


   60    CONTINUE

         CALL GSYMDA(1)
         CALL GSYMSY(2, 2)


         XSC = MAXVEC / VECSCL

         YSC = XSC * YSIZE / (XSIZE * YEXAG)

         CALL GSYMSC(1, 1, XSC)

         CALL GSYMDI(XFEAT, YFEAT, ZOR, 2, FENUM)

      ENDIF

C
C-----------------------------------------------------------------------
C
C PLOT BOUNDARY
C

      IF(IBN .GT. 0) THEN
         ICM = 1
         CALL GWICOL(0.002*XSIZE, 32)
         IF(IPBB .EQ. 2) CALL GWICOL(0.002*XSIZE, 4)
   70    CONTINUE
         IF(ICM+1 .LE. IBN) THEN
           CALL GVECT(XBB(ICM), YBB(ICM), 0)
   80      CONTINUE 
           IF(IBB(ICM) .EQ. IBB(ICM+1) .AND.
     +        ABS(XBB(ICM)-XBB(ICM+1)) .LT. CTOL) THEN
              CALL GVECT(XBB(ICM+1), YBB(ICM+1), 1)
              ICM = ICM+1
              GOTO 80
           ELSE
              ICM = ICM+1
              GOTO 70
           ENDIF
         ENDIF
      ENDIF

C
C-----------------------------------------------------------------------
C
C PLOT AND LABEL AXES
C
      CALL RAXDIS(6, 1, 0)
      IF(NAL .EQ. 0) CALL RAXDIS(4, 0, 0)
      IF(NOAX .EQ. 0) THEN
         CALL RAXTEA(6, 2, 2, RDEF, IDEF, IDEF, RDEF)

         CALL RAXTEX(6, -1, XLAB, RDEF, RDEF, 0.6*HEIGHT)
         CALL RAXIS(1, YMIN, HEIGHT, 1)
         CALL RAXTEX(6, -1, YLAB, RDEF, RDEF, 0.6*HEIGHT)
         CALL RAXIS(2, XMIN, HEIGHT, 1)
      ENDIF


C
C------------------------------------------------------------------------
C
C Plot labeling 
C
      CALL RTXCOL(2,2)
      CALL RTXHEI(2.0)
C
C For data
C
      CALL RTX(-1, 'FRAME ', XMIN, YMAX + 0.1*(YMAX - YMIN))
      CALL RTXNC(FNUM, 0)
      CALL RTXC(-1, '  THRESHOLD = ')
      CALL RTXNBO(5, 'E')
      CALL RTXNC(TH, 4) 
C
C For Title
C
      CALL RTXHEI(3.0)
      CALL RTX(-2, TITLE1, XMIN, YMAX + 0.2*(YMAX - YMIN))

      CALL GEMPTY
      CALL GSEGCL(K)
      CALL GCLOSE
      RETURN
  792 FORMAT('Do you want a unique output directory for the graphics?')
  795 FORMAT('Process Id is ', I8.8);
      END


#endif


   
