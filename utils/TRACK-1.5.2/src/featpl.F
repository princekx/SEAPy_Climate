
/* fortran preprocessor file for stats. plot using uniras */

#ifndef     HAVEUNIRAS

      SUBROUTINE FEATPL()

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                     C
C  Dummy routine for FEATPL if no UNIRAS available.                   C
C                                                                     C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      INCLUDE 'paramc.ptr'


#ifdef    NAMELISTS

      READ(NAMCHN, NML=PARAMC)
      REWIND(NAMCHN)

#endif

      WRITE(NUT, 900)

      RETURN
  900 FORMAT('/***ERROR***, no plotting possible without UNIRAS/')

      END

#else

      SUBROUTINE FEATPL(XGRID, YGRID, IX, IY, XP, YP, ZP, TID, DFIL, CMI
     +                 ,CMXG, CMYG, DCM, K, W1, W2, DMAX, PHIMAX, IZM
     +                 ,XMN, YMN, XMX, YMX, DELB, IARR, NTR, IPLVEC,
     +                 XOR, YOR, ZOR) 
      IMPLICIT NONE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                       C
C Program to plot the feature points.                                   C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INTEGER IX, IY, DFIL, DCM, CMI(DCM), IZM, DELB, NTR
      INTEGER ICM, IC, I, TID(DFIL), K, IARR, ICOL
      INTEGER IPLVEC

      REAL W1, W2, DMAX, PHIMAX, PERIOD, CTOL

      REAL XP(DFIL), YP(DFIL), ZP(DFIL)
      REAL XOR(DFIL), YOR(DFIL), ZOR(2*DFIL)
      REAL CMXG(DCM), CMYG(DCM)
      REAL XMN, YMN, XMX, YMX, PTX, PTY, XD, YD, RAT, XS, SCL
      REAL YEXAG, MAXVEC
      REAL XSC, YSC
      CHARACTER TITLE*100, FIL*70
      CHARACTER ANS*3
C
C------------------------------------------------------------------------
C
C Include the parameter file.
C
      INCLUDE 'paramf.ptr'
      INCLUDE 'paramc.ptr'
C
      REAL XGRID(IX), YGRID(IY)
      REAL XMIN, XMAX, YMIN, YMAX, XSI, YSI, XM, 
     +     XOFF, YOFF, XSIZE, YSIZE, HEIGHT, DM, ZL(2), SCY
      REAL PII, PI2
      INTEGER IANS, IATYP1, IATYP2
      LOGICAL FSTAT, SWTCH1, SWTCH2, SWTCH
      CHARACTER TITSUB*100

      PARAMETER(PII=57.295828, PI2=1.570796)

      INTEGER IPID

      INTEGER GETPID, SYSTEM
      CHARACTER NEWFIL*70, UDIR*70, CPID*8

      INTEGER     NCOL, N
      REAL        COLOR(3, 201), CLR(3), HATCH(5)
      CHARACTER CTFIL*70


#ifdef    NAMELISTS

      READ(NAMCHN, NML=PARAMF)
      READ(NAMCHN, NML=PARAMC)
      REWIND(NAMCHN)

#endif

C
C--------------------------------------------------------------------
C
C  Read color table information if required.
C
      IF(IACT .EQ. 1) THEN
         COLOR(1,1) = 0
         COLOR(2,1) = 0
         COLOR(3,1) = 0
         CTFIL=PATH // '../data/rgb_ct.dat'  
         OPEN(30, FILE=CTFIL, STATUS='OLD')
         READ(30, *) NCOL
         DO 2 I=1,NCOL
            READ(30, *) (COLOR(N,I+1),N=1,3)
    2    CONTINUE
         CLOSE(30)
      ENDIF
C
C---------------------------------------------------------------------
C
C Get process ID for unique output directory for graphics
C
      WRITE(NUT, 792)
      READ(NIN, *) ANS
      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN

         IPID = GETPID()
         WRITE(NUT, 795) IPID
         WRITE(CPID,'(I8.8)') IPID
         UDIR = 'mkdir ' // PATH // '/' // CPID
         ipid=system(UDIR)
         IF(K .LT. 0) THEN
            NEWFIL=PATH // '/' // CPID // '/' // SPLICE
            SPLICE = NEWFIL
         ELSE
            NEWFIL=PATH // '/' // CPID // '/' // FILFET
            FILFET = NEWFIL
         ENDIF
      ELSE
         IF(K .LT. 0) THEN
            SPLICE = PATH // '/' // SPLICE
         ELSE
            FILFET = PATH // '/' // FILFET
         ENDIF
      ENDIF

      FSTAT= .TRUE. 
      SCL=SCALE
      PERIOD = XGRID(IX) - XGRID(1)
      CTOL = PERIOD / 4.0

      IF(IARR .EQ. 0) THEN
         IATYP1 = ARRHTP
         IATYP2 = 0
      ELSE
         IATYP1 = 0
         IATYP2 = ARRHTP
      ENDIF
      
      IF(NVEC .EQ. 0) THEN
         IATYP1 = 0
         IATYP2 = 0
      ENDIF

      WRITE(NUT, 850)
      READ(NIN, *)ANS
      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
        WRITE(NUT, 860)
        READ(NIN, *) SCL 
      ENDIF
      ICOL = ICNUM 
      ZL(1) = ZINIT
      ZL(2) = STEP
      WRITE(NUT, 790)ICNUM, ZINIT, STEP
      WRITE(NUT, 800)
      READ(NIN, *)ANS
      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
        WRITE(NUT, 810)
        READ(NIN, *) ICOL, ZL(1), ZL(2)
      ENDIF
    
C
C------------------------------------------------------------------------
C
C Check if output file's already exists for plotting.
C

      IF(K .LT. 0) THEN
        INQUIRE(FILE=SPLICE, EXIST=FSTAT)
        FIL = SPLICE
        TITLE = TITLE3
        K = 1
        CALL GROUTE(DUMMY)
        CALL GFILE(10, NOUT, FIL)
        CALL GOPEN
    5   CONTINUE
        CALL GSEGEX(K, IANS)
        IF(IANS .GT. 0) THEN
          K = K + 1
          GOTO 5
        ENDIF
        CALL GCLOSE
      ELSE 
        INQUIRE(FILE=FILFET, EXIST=FSTAT)
        FIL = FILFET
        TITLE = TITLE2
        CALL GROUTE(DUMMY)
        CALL GFILE(10, NOUT, FIL)
        CALL GOPEN
    6   CONTINUE
        CALL GSEGEX(K, IANS)
        IF(IANS .GT. 0) THEN
          K = K + 1
          GOTO 6
        ENDIF
        CALL GCLOSE
      ENDIF

      DM = SQRT(DMAX)
C
C
C------------------------------------------------------------------------
C
C Set max. and min. values of the grid.
C

      IF(IZM .EQ. 0) THEN
        XMIN=XGRID(1)
        XMAX=XGRID(IX)
        YMIN=YGRID(1)
        YMAX=YGRID(IY)
      ELSE IF(IZM .EQ. 1) THEN
        XD = XMX-XMN
        YD = YMX-YMN
        RAT = REAL(DELB)/REAL(IX)
        XMIN=XMN-RAT*XD
        XMAX=XMX+RAT*XD
        RAT = REAL(DELB)/REAL(IY)
        YMIN=YMN-RAT*YD
        YMAX=YMX+RAT*YD
      ENDIF
C
C------------------------------------------------------------------------
C
C Put data into the plotting array
C
      DO 10 I=1, DFIL
        ZP(I) = SCL*ZP(I)
   10 CONTINUE
C
C=========================================================================
C
C Plot data using UNIRAS.
C
      CALL GROUTE(DEVICE)
      CALL GFILE(10, NOUT, FIL)
      CALL GOPEN
      IF(FSTAT)  CALL GSEGWK(0)
      CALL GSEGCR(K)
      CALL GRPSIZ(XSI,YSI)
      XM=MIN(XSI,YSI)
      XSIZE=0.85*XM
      YSIZE=XSIZE
      IF(SCALEY .LE. 0.) THEN
         SCY = SCALEX * (YMAX - YMIN) / (XMAX - XMIN)
      ELSE
         SCY = SCALEY
      ENDIF

      XOFF=SOFFX*(XSI-XSIZE)
      YOFF=SOFFY*(YSI-YSIZE)
      CALL GLIMIT(XMIN, XMAX, YMIN, YMAX, 0.0, 0.0)
      CALL GVPORT(XOFF, YOFF, SCALEX*XSIZE, SCY*YSIZE)

      IF(IZM .EQ. 1) CALL GWBOX(SCALEX*XSIZE, SCY*YSIZE, 0.0)
      CALL GSCALE
      HEIGHT = 0.02*AMIN1(XSIZE,YSIZE)
      CALL RRECT(XMIN, YMIN, XMAX, YMAX, 0, 0)
      CALL RRECT(XMN, YMN, XMX, YMX, 1, 0)
      CALL RUNDEF(RDEF, IDEF)
C
C-----------------------------------------------------------------------
C
C Define color table
C

      IF(IACT .EQ. 1) THEN 
         CALL RCTCRE(2,NCOL,3,2,255.)
         CALL RCTA(2)
         DO 30 I=1,NCOL+1 
            CLR(1) = COLOR(1, I)
            CLR(2) = COLOR(2, I)
            CLR(3) = COLOR(3, I)
            CALL RCOLOR(0,I,CLR,HATCH)
   30    CONTINUE
      ENDIF
C
C------------------------------------------------------------------------
C
C Plot country map.
C
      ICM = 1
c      XS = 0.001*XSIZE
      XS = 0.0025*XSIZE
   40 CONTINUE
      IF(ICM+1 .LE. DCM) THEN 
       PTX = (CMXG(ICM) - XMN)*(XMX - CMXG(ICM))
       PTY = (CMYG(ICM) - YMN)*(YMX - CMYG(ICM))
       IF((PTX .GE. 0.0) .AND. (PTY .GE. 0.0)) THEN
         SWTCH1 = .TRUE.
       ELSE
         SWTCH1 = .FALSE.
       ENDIF
   50  CONTINUE   
       IC = ICM + 1

       PTX = (CMXG(IC) - XMN)*(XMX - CMXG(IC))
       PTY = (CMYG(IC) - YMN)*(YMX - CMYG(IC))
       IF((PTX .GE. 0.0) .AND. (PTY .GE. 0.0)) THEN
         SWTCH2 = .TRUE.
       ELSE
         SWTCH2 = .FALSE.
       ENDIF
       IF(SWTCH1 .OR. SWTCH2 ) THEN
         SWTCH = .TRUE.
       ELSE 
         SWTCH = .FALSE.
       ENDIF
       IF(CMI(ICM) .EQ. CMI(IC) .AND. 
     +    ABS(CMXG(ICM)-CMXG(IC)) .LT. CTOL) THEN
          IF(IZM .EQ. 1) THEN
             IF(SWTCH)THEN
               CALL GWICOL(XS, 0)
               CALL GVECT(CMXG(ICM), CMYG(ICM), 0)
               CALL GVECT(CMXG(IC), CMYG(IC), 1)
             ENDIF
          ELSE
             IF(SWTCH) THEN
               CALL GWICOL(XS, 0)
               CALL GVECT(CMXG(ICM), CMYG(ICM), 0)
               CALL GVECT(CMXG(IC), CMYG(IC), 1)
             ELSE
               CALL GWICOL(XS, 1)
               CALL GVECT(CMXG(ICM), CMYG(ICM), 0)
               CALL GVECT(CMXG(IC), CMYG(IC), 1)
             ENDIF
          ENDIF
          ICM = ICM+1
          SWTCH1 = SWTCH2
          GOTO 50
       ELSE
          ICM = ICM+1
          GOTO 40
        ENDIF
      ENDIF

C
C-----------------------------------------------------------------------
C
C Plot feature points.
C

      IF(NTR .EQ. 1 ) GOTO 500

      ICM = 1
      SWTCH = .TRUE.
  200 CONTINUE


      IF(ICM .LE. DFIL) THEN
       IF(IZ) THEN
          PTX = (XP(ICM) - XMN)*(XMX - XP(ICM))
          PTY = (YP(ICM) - YMN)*(YMX - YP(ICM))

          IF((PTX .GE. 0.) .AND. (PTY .GE. 0.)) THEN
             SWTCH1 = .TRUE.
          ELSE
             SWTCH1 = .FALSE.
             ZP(ICM) = RDEF
          ENDIF
       ENDIF
  300  CONTINUE
       IC = ICM + 1
       IF(IC .LE. DFIL) THEN 
         IF(TID(ICM) .EQ. TID(IC)) THEN

            IF(IZ) THEN
              PTX = (XP(IC) - XMN)*(XMX - XP(IC))
              PTY = (YP(IC) - YMN)*(YMX - YP(IC))
              IF((PTX .GE. 0.0) .AND. (PTY .GE. 0.0)) THEN
                 SWTCH2 = .TRUE.
              ELSE
                 SWTCH2 = .FALSE.
                 ZP(IC) = RDEF
              ENDIF
              IF(SWTCH1 .OR. SWTCH2 ) THEN
                 SWTCH = .TRUE.
              ELSE 
                 SWTCH = .FALSE.
              ENDIF
            ENDIF

            IF(SWTCH) THEN
C
C Check for wrap-round
C 
              IF(ABS(XP(ICM)-XP(IC)) .GT. CTOL) THEN
C                 IF(PERIOD-XP(ICM) .LT. TOLB) THEN 
C                   CALL GARROW(XP(IC), YP(IC), 3, XP(ICM)-PERIOD, 
C     +                    YP(ICM), 0, 27, 27, 0.001*XSIZE)
C                   CALL GARROW(XP(IC)+PERIOD, YP(IC), 3,  
C     +                    XP(ICM), YP(ICM), 0, 27, 27, 0.001*XSIZE)
C                 ELSE
C                   CALL GARROW(XP(IC), YP(IC), 3, XP(ICM)+PERIOD, 
C     +                    YP(ICM), 0, 27, 27, 0.001*XSIZE)
C                   CALL GARROW(XP(IC)-PERIOD, YP(IC), 3,  
C     +                    XP(ICM), YP(ICM), 0, 27, 27, 0.001*XSIZE)
C                 ENDIF

               ELSE
        
                 CALL GARROW(XP(IC), YP(IC), IATYP1, 
     +             XP(ICM), YP(ICM), IATYP2, NTC, NTC, ARRSCL*SSC*XSIZE)

               ENDIF


             ENDIF

             ICM = ICM+1
             SWTCH1 = SWTCH2
             GOTO 300
         ELSE
	 
            ICM = ICM+1
            GOTO 200
         ENDIF
       ENDIF
      ENDIF

  500 CONTINUE

      IF(IZM .EQ. 1) THEN

C
C-----------------------------------------------------------------------
C
C Plot orientation vectors.
C

         MAXVEC = 0.
         YEXAG = 1.

         IF(IPLVEC .GT. 0) THEN

            DO 600 I=1, DFIL
               ZOR(2*I-1) = SQRT(XOR(I) * XOR(I) + YOR(I) * YOR(I) 
     +                           * YEXAG * YEXAG)

               MAXVEC = MAX(MAXVEC, ZOR(2*I-1))

               ZOR(2*I) = ATAN2(YOR(I) * YEXAG, XOR(I)) * PII

               IF(ZOR(2*I-1) .LT. SUPVEC) THEN
                 ZOR(2*I-1) = RDEF
               ENDIF


  600       CONTINUE



            CALL GSYMDA(1)
            CALL GSYMSY(2, 2)
            CALL GSYMCO(3, 0.003*SSC*XSIZE)
            CALL GSYMCO(4, 0.008*SSC*XSIZE)
            CALL GSYMCO(5, 0.008*SSC*XSIZE)
            CALL GSYMCO(7, 32)


            XSC = MAXVEC / VECSCL

            YSC = XSC * YSIZE / (XSIZE * YEXAG)

            CALL GSYMSC(1, 1, XSC)

            CALL GSYMDI(XP, YP, ZOR, 2, DFIL)

         ENDIF
C
C-----------------------------------------------------------------------
C
C Plot feature points.
C

        IF(NPNT .GT. 0) THEN
           CALL GWICOL(PTSCL*SSC*XSIZE, 9999)
           CALL RCLASS(ZL,ICOL,-5)
           CALL RSHADE(ITAB,0)
           CALL GPOI2Z(XP, YP, ZP, DFIL)

           IF(NCSC .GE. 0) THEN
             CALL GCOSCL(XMAX+0.06*(XMAX-XMIN), YMIN)
           ENDIF
	ENDIF

      ENDIF

C
C-----------------------------------------------------------------------
C
C PLOT AND LABEL AXES
C
      CALL RAXDIS(6, 1, 0)
      IF(NAL .EQ. 0)CALL RAXDIS(4, 0, 0)
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
C      CALL RTX(-1, 'COST FUNCTION WEIGHTS: W1 = ', 
C     +         XMIN, YMAX+0.1*(YMAX-YMIN))
C      CALL RTXNBO(5, 'E')
C      CALL RTXNC(W1, 0)
C      CALL RTXC(-1, ', W2 = ')
C      CALL RTXNC(W2, 0)
C      CALL RTXC(-1,', OPT. PARAM.S: DMAX = ')
C      CALL RTXNC(DM, 0)
C      CALL RTXC(-1,', PHIMAX = ')
C      CALL RTXNC(PHIMAX, 0)

C
C For Title
C
      CALL RTXHEI(3.0)
      CALL RTX(-2, TITLE, XMIN, YMAX + 0.15*(YMAX - YMIN))

      CALL GEMPTY
      WRITE(NUT, 900)
      READ(NIN, *) ANS
      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
         WRITE(NUT, 910)
         READ(NIN, '(A100)') TITSUB
         CALL RTXHEI(2.5)
         CALL RTX(-2, TITSUB, XMIN, YMAX+0.08*(YMAX-YMIN))
      ENDIF 
      CALL GEMPTY

      CALL GSEGCL(K)
      CALL GCLOSE
      RETURN
  790 FORMAT('There are currently ', I2, ' bins starting at ',
     +       F10.4, ' with step size ', F10.4)
  792 FORMAT('Do you want a unique output directory for the graphics?')
  795 FORMAT('Process Id is ', I8.8);
  800 FORMAT('do you wish to change the colour scale parameters')
  810 FORMAT('input the new colour scale parameters, ICNUM,'
     +       ' ZINIT and STEP')
  850 FORMAT('do you wish to change the field scaling')
  860 FORMAT('input the new field scaling') 
  900 FORMAT('do you want to use a subsidary title to the default')
  910 FORMAT('input the subsidary title')
      END

#endif




   
