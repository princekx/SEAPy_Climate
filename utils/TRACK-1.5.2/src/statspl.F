
/* fortran preprocessor file for stats. plot using uniras */

#ifndef     HAVEUNIRAS

      SUBROUTINE STATSPL()

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                     C
C  Dummy routine for STATSPL if no UNIRAS available.                  C
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

      SUBROUTINE STATSPL(XP, YP, SP, ZP, IAR, IBR, IM, ISN, IS, ISS,  
     +                   ISTY,IFTY, CMI, CMXG, CMYG, DCM, XA1, XA2,  
     +                   YA1, YA2, CT, GRID, NPX, NPY, SMTH, KERN,
     +                   IPGR, TDGRD)

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                                                                     C
C Routine to plot contour maps of stat. data using UNIRAS             C
C                                                                     C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C

      INTEGER ISN, IS, ISS, ISTY(IS), IFTY(IS), DCM, CMI(DCM)
      INTEGER NPX, NPY,  IM, IB, IV, IAR, IBR, IPGR, ICOL, ITCOL
      INTEGER IANS, I, J, K, N, M, L
      INTEGER INPS, INPR, NSEG, INCL
      INTEGER ISEG, IPR, IDREG, IREG, ICDIM
      INTEGER DCTY, AST, IL, NCHAR(3)
      INTEGER KERN(IS)
      INTEGER IVP
      INTEGER CCOL, IADD
      INTEGER IVLN
      INTEGER MCOL
      INTEGER IFMNX
      INTEGER ITD, IFD
      INTEGER ITPT, IFPT, ITOT
      INTEGER ILNSTY, INSTY
      INTEGER ISHTY
      PARAMETER(MCOL=20)

      REAL XP(IBR), YP(IBR), ZP(ISN * IM), SP(ISN*ISS), SMTH(IS)
      REAL CMXG(DCM), CMYG(DCM), GRID(NPX, NPY), TDGRD(NPX, NPY)
      REAL XA1, XA2, YA1, YA2, CT
      REAL ZL(MCOL), ZLT(MCOL), RAT
      REAL XSI, YSI
      REAL XOF, YOF, XOFF, YOFF, XCEN, YCEN
      REAL XSIZE, YSIZE, LPX, LPY
      REAL HEIGHT, XS, ZMAX, ZMIN
      REAL YEXAG, PII, PI2
      REAL XSC, YSC, MAXX
      REAL FSCALE, FSS
      REAL VVSC
      REAL DTHR

      LOGICAL FSTAT, PGRTY
      CHARACTER ANS*3, CONT*1, TITSUB*100
      CHARACTER*5 TT(3)
      CHARACTER NEWFIL*70, UDIR*70, CPID*8

      PARAMETER(PII=57.295828, PI2=1.570796)
      PARAMETER(ICDIM=100)

      INTEGER SDARRY(ICDIM), ICNW
      INTEGER IUCOL
      INTEGER IPID

      INTEGER GETPID, SYSTEM

C
C---------------------------------------------------------------------
C
C Includes
C
      INCLUDE 'paramc.ptr'
      INCLUDE 'params.ptr'
      INCLUDE 'params_lab.ptr'

#ifdef    NAMELISTS

      READ(NAMCHN, NML=PARAMC)
      READ(NAMCHN, NML=PARAMS)
      READ(NAMCHN, NML=PARAMS_LAB)
      REWIND(NAMCHN)

#endif

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
         NEWFIL=PATH // '/' // CPID // '/' // FILNM
         ipid=system(UDIR)
         FILNM = NEWFIL

      ELSE
         FILNM = PATH // '/' // FILNM
      ENDIF

      IUCOL = 0
      ITOT = ISN * ISS

      IVLN = 1
      VVSC = VSCALE

      ITD = 0
      IFD = 0

      FSS = 1.0

      NCHAR(1) = 5
      TT(1) = "below"
      NCHAR(2) = 5
      TT(2) = "above"
      NCHAR(3) = 0

      FSTAT = .TRUE.
      PGRTY = .FALSE.

      IF(IPGR .GT. 0) PGRTY = .TRUE.

      RAT = (YA2 - YA1) / (XA2 - XA1)
      ICOL = CNUM
      ZL(1) = CINIT
      ZL(2) = CSTEP

      WRITE(NUT, 790)CNUM, CINIT, CSTEP
      WRITE(NUT, 800)
      READ(NIN, *)ANS
      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
        WRITE(NUT, 810)
        READ(NIN, *)ICOL, ZL(1), ZL(2)
      ENDIF
C
C--------------------------------------------------------------------
C
C Check if output file exists and how many segments are in use
C
      INQUIRE(FILE=FILNM, EXIST=FSTAT)
      K = 1
      CALL GROUTE(DUMMY)
      CALL GFILE(10, NOUT, FILNM)
      CALL GOPEN

    1 CONTINUE
      CALL GSEGEX(K, IANS)
      IF(IANS .GT. 0) THEN
         K = K + 1
         GOTO 1
      ENDIF
      CALL GCLOSE

C
C--------------------------------------------------------------------
C
C Produce required plots.
C

      CALL GROUTE(DEVICE)
      CALL GFILE(10, NOUT, FILNM)
      CALL GOPEN
      IF(FSTAT) CALL GSEGWK(0)

      CALL GRPSIZ(XSI, YSI)
      CALL RUNDEF(RDEF, IDEF)
      CALL GCONSM(CSMTH)

      DO 2 I=1, CDIM

         CWI(I) = CONLW
         CLSTY(I) = CONLS
	 CAHGT(I) = 0.0
c         IF(AMOD(REAL(I), CNB) .LT. 0.00001) THEN
c           CWI(I) = 2.0 * CWI(I)
c           CLSTY(I) = 0
c         ENDIF

    2 CONTINUE 

      CALL RCLASS(ZL, ICOL, -5)

      IF(ITAB .GT. 6) CALL COLDEF(ITAB)

      XSI = XSI * SCALEX
      YSI = YSI * SCALEY


      WRITE(NUT, 900)
      READ(NIN, *) INPS

      NSEG = IS / INPS

      IF(NSEG*INPS .LT. IS) NSEG = NSEG + 1

      WRITE(NUT, 910)
      READ(NIN, *) INPR

      INCL = INPS / INPR

      IF(INPR*INCL .LT. INPS) INCL = INCL + 1


      IF(INPR .GE. INCL) THEN        
         IF(RAT .LT. 1.) THEN

            LPX = XSI / INPR
            LPY = RAT * LPX

         ELSE

            LPY = YSI / INCL
            LPX = LPY / RAT

         ENDIF

      ELSE

         IF(RAT .LT. 1.) THEN

           LPY = YSI / INCL
           LPX = LPY / RAT

         ELSE

           LPX = XSI / INPR
           LPY = LPX * RAT

         ENDIF

      ENDIF

      XCEN = (XSI - INPR * LPX) / 2.
      YCEN = (YSI - INCL * LPY) / 2.

      IF(NTEXT .EQ. 0 .AND. (NCSC .LT. 0 .OR. NCSC .GT. 2)) THEN
        XSIZE = 0.98 * LPX
        YSIZE = 0.98 * LPY
        XOF = 0.01 * AMIN1(LPX, LPY)
        YOF = XOF
      ELSE
        XSIZE = 0.8 * LPX
        YSIZE = 0.8 * LPY
        XOF = 0.06 * AMIN1(LPX, LPY)
        YOF = XOF
      ENDIF

      IF(PGRTY) THEN
        YEXAG = 1.0
      ELSE
        YEXAG = YSIZE / (XSIZE * RAT)
      ENDIF

      XS = 0.0025 * AMAX1(XSIZE, YSIZE)


      HEIGHT = 0.02*AMIN1(XSIZE,YSIZE)
C
C------------------------------------------------------------------------
C
C Set character escape's for superscript and subscript etc.
C
      CALL RTXESC(1, '~')
      CALL RTXESC(2, '^')
      CALL RTXESC(5, '#')

C
C------------------------------------------------------------------------
C
C  Set default contouring type
C

      WRITE(NUT, 920)
      WRITE(NUT, 930)
      READ(NIN, *) DCTY

      WRITE(NUT, 940)
      READ(NIN, *) CONT

      IPR = 1
      AST = 0
C
C------------------------------------------------------------------------
C
C Set contour line widths, annotation and line styles
C

      CALL RCTMOD(1, 255, 3, 3, 100.) 

      CALL GWICOL(0.002*XSIZE, 1)
      CALL GCONWI(CWI, ICOL)
      CALL GCONA(CHS*HEIGHT, 2, -0.8*XSIZE, NANOT)
      CALL GCONDA(CLSTY, ICOL)
      CALL GCONAB(IFILL, IFR)
      CALL RSHADE(ITAB, 0)

C
C------------------------------------------------------------------------
C
C Create country mask.
C

      IF(IMASK) THEN

         WRITE(NUT,941)

         CALL MASK(DCM, CMI, CMXG, CMYG, XA1, XA2, YA1, 
     +             YA2, CT, XS, IRREG, IRN, IREG)


      ENDIF


C
C------------------------------------------------------------------------
C
C Check for feature and track densities.
C
      DO 5 I = 1, IS

         IF(ISTY(I) .EQ. 5 .OR. ISTY(I) .EQ. 16) THEN
            IF(ISTY(I) .GT. IFD) THEN
               IFD = ISTY(I)
               IFPT = I
            ENDIF
         ENDIF

         IF(ISTY(I) .EQ. 8 .OR. ISTY(I) .EQ. 17) THEN
            IF(ISTY(I) .GT. ITD) THEN
              ITD = ISTY(I)
              ITPT = I
            ENDIF
         ENDIF

   5  CONTINUE

      IF(IFD .GT. 0) THEN

         IADD = 0

         DO 6 M = 1, IS

            IF(IFTY(M) .GT. 1)THEN

               IF(IFPT .GT. M) IADD = IADD + ISN

            ENDIF

    6    CONTINUE

         IFPT = (IFPT - 1) * ISN + IAR + IADD

         IFD = 1

      ENDIF

      IF(ITD .GT. 0) THEN

         IADD = 0

         DO 8 M = 1, IS

             IF(IFTY(M) .GT. 1)THEN

                IF(ITPT .GT. M) IADD = IADD + ISN

             ENDIF

    8    CONTINUE

         ITPT = (ITPT - 1) * ISN + IAR + IADD

         ITD = 1

      ENDIF


      DO 20 I = 1, NSEG

         ISEG = 1

         CALL GSEGCR(K)

         IF(NOBND .GT. 0) 
     +      CALL RRECT(XCEN, YCEN, XCEN+INPR*LPX, YCEN+INCL*LPY, 0, 2)

         DO 30 J = 1, INCL

            YOFF = YCEN + (J - 1) * LPY

            DO 40 N = 1, INPR

               IB = -1

               IF(IPR .GT. IS) GOTO 25
               IF(ISEG .GT. INPS) GOTO 25

               XOFF = XCEN + (N - 1) * LPX

               IF(NOBND .GT. 0) 
     +            CALL RRECT(XOFF, YOFF, XOFF+LPX, YOFF+LPY, 0 , 2)

               CALL GLIMIT(XA1, XA2, YA1, YA2, 0., 0.)
               CALL GVPORT(XOFF+XOF, YOFF+YOF, XSIZE, YSIZE)
               CALL RRECT(XOFF+XOF, YOFF+YOF, XOFF+XOF+XSIZE, 
     +                    YOFF+YOF+YSIZE, 0 , 2)

               CALL GSCALE

               CALL GMAPFR(10)

C
C--------------------------------------------------------------------
C
C Set region
C
               IF(IBR .GT. 1)THEN
                  CALL GREGCR(XP, YP, IBR, IDREG)
                  CALL GREGSO(IDREG, 1)
                  CALL GREGSP(1)
                  CALL GREGSS(IDREG, 1)
               ENDIF

               IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                  WRITE(NUT, '("Current plot is ", A30)')
     +              SLAB(ISTY(IPR))(1:SLLT(ISTY(IPR)))

                  WRITE(NUT, 930)
                  READ(NIN, *) DCTY

               ENDIF

C
C-----------------------------------------------------------------------
C
C Put data into plotting array for gridding
C
               IF(IFTY(IPR) .EQ. 1) THEN

                 IF(DCTY .EQ. 3 .OR. DCTY .EQ. 4) THEN

                    WRITE(NUT, 945)
                    DCTY = 0
C                    AST = AST + 2 * ISN
C                    GOTO 100

                 ENDIF

                 IF(DCTY .EQ. 5 .OR. DCTY .EQ. 6) THEN


                    WRITE(NUT, 948)
                    DO 31 M = 1, IS

                       IF(IFTY(M) .EQ. 1) 
     +                   WRITE(NUT, '(A30, 2X, I2)')
     +                      SLAB(ISTY(M))(1:SLLT(ISTY(M))), M

   31               CONTINUE

                    WRITE(NUT, 949)
                    READ(NIN,*)IB
                    IF(IB .LE. 0 .OR. IB .GT. M) THEN

                       WRITE(NUT, 951)
                       GOTO 100

                    ENDIF

                    IADD = 0

                    DO 32 M = 1, IS

                       IF(IFTY(M) .GT. 1)THEN

                          IF(IB .GT. M) IADD = IADD + ISN

                       ENDIF

   32               CONTINUE

                    WRITE(NUT, 952) SLAB(ISTY(IB))(1:SLLT(ISTY(IB)))

                    IV = (IB - 1) * ISN + IADD

                    ZMIN = SP(IV + 1)
                    ZMAX = ZMIN

                    FSCALE = FSS

                    IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                      WRITE(NUT, 780)
                      READ(NIN, *)ANS
                      IF(ANS(1:1) . EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                THEN
                        WRITE(NUT, 785)
                        READ(NIN, *)FSCALE
                        IF(FSCALE .LT. 0) FSCALE = ABS(FSCALE)
                      ENDIF

                    ENDIF

                    ZMAX = FSCALE * ZMAX
                    ZMAX = ZMIN

                    DO 33 M = 1, IBR

                      ZP(M) = FSCALE * SP(IV + IAR + M)

                      IF(ZP(M) .GT. ZMAX) ZMAX = ZP(M)
                      IF(ZP(M) .LT. ZMIN) ZMIN = ZP(M)

   33               CONTINUE


                    WRITE(NUT, 960) ZMIN, ZMAX

                    IF(DCTY .EQ. 5) THEN

                       WRITE(NUT, 962)
                       READ(NIN, *) ITCOL
                       IF(ITCOL .GT. BSDIM) THEN
                         WRITE(NUT,963) ITCOL, BSDIM
                         ITCOL = BSDIM
                       ENDIF
                       DO 34 M =1, ITCOL
                          WRITE(NUT, '("Input Limit",2X, I2)') M
                          READ(NIN,*) ZLT(M)
   34                  CONTINUE


                       CALL RCLASS(ZLT, ITCOL, 0)

                       CALL RCSDEF(ITABB, BSCALE, NUMCOL, CMODE, 
     +                             COLMAX)

                       CALL RSHADE(ITABB, 0)

                    ELSE

                       WRITE(NUT, 800)
                       READ(NIN, *)ANS
                       IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                 THEN

                         WRITE(NUT, 805)
                         READ(NIN, *)ANS
                         IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                   THEN
                            WRITE(NUT,'("Input number of classes")')
                            READ(NIN, *)ITCOL
                            CALL RCLSCA(ITCOL, IPERC, 0)
                         ELSE
                            WRITE(NUT, 810)
                            READ(NIN, *) ITCOL, ZLT(1), ZLT(2)
                            CALL RCLASS(ZLT, ITCOL, -5)
                         ENDIF
                       ELSE
                         CALL RCLASS(ZLT, ITCOL, -5)
                       ENDIF

                       CALL RSHADE(ITAB, 0)

                    ENDIF

                    CALL GIMETH(ITY)
                    CALL GINTP1(XP, YP, ZP, IBR, GRID, NPX, NPY)

                    CALL SUPPRES(GRID, TDGRD, XP, YP, ZP, SP, IBR, NPX,   
     +                   NPY, IFD, ITD, ISTY(IB), IFPT, ITPT, NIN, NUT, 
     +                   ITOT, DCTY, RDEF)


                    WRITE(NUT, 815)
                    READ(NIN, *)ISHTY
                    IF(ISHTY .EQ. 1) THEN          
                      CALL GCONWI(CWI, ICOL)
                      CALL GCONDA(CLSTY, -1)
                      CALL GCNR2V(GRID, NPX, NPY)
                    ELSE IF(ISHTY .EQ. 2) THEN
                      CALL GCONWI(CWI, -1)
                      CALL GCNR2S(GRID, NPX, NPY)
                    ELSE
                      CALL GCONDA(CLSTY, -1)
                      CALL GCONWI(CWI, ICOL)
                      CALL GCNR2S(GRID, NPX, NPY)
                    ENDIF

                    CALL GCONWI(CWI, ICOL)

                    IF(IBBX .EQ. 0) THEN
                      CALL GREGDR(IDREG, -1, 0.0, 1, 0)
                    ELSE
                      CALL GREGDR(IDREG, -1, -1.0, 1, 0)
                    ENDIF

                    IF(NCSC .EQ. 0) THEN
                      CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                    ELSE IF(NCSC .EQ. 1) THEN
                      CALL GCOSCS(XA2+(XA2-XA1)*0.01, 
     +                            YA1-(YA2-YA1)*0.04)
                    ELSE IF(NCSC .EQ. 2) THEN
                      CALL GCOSLI(XA2+(XA2-XA1)*0.01, 
     +                            YA1-(YA2-YA1)*0.04)
                    ENDIF

                    WRITE(NUT, 953) SLAB(ISTY(IPR))(1:SLLT(ISTY(IPR)))

                 ENDIF

                 FSCALE = FSS

                 IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                    WRITE(NUT, 780)
                    READ(NIN, *)ANS
                    IF(ANS(1:1) . EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
                       WRITE(NUT, 785)
                       READ(NIN, *)FSCALE
                       IF(FSCALE .LT. 0) FSCALE = ABS(FSCALE)
                    ENDIF

                 ENDIF

                 IFMNX = 0

                 DO 50 M = 1, IBR

                   IF(ABS(SP(AST + IAR + M) - RDEF) .GT. 1.0E-4)THEN
                     IF(IFMNX .EQ. 0) THEN
                       ZMIN = FSCALE * SP(AST + IAR + M)
                       ZMAX = ZMIN
                       IFMNX = 1
                     ENDIF
                     ZP(M) = FSCALE * SP(AST + IAR + M)
                     IF(ZP(M) .GT. ZMAX) ZMAX = ZP(M)
                     IF(ZP(M) .LT. ZMIN) ZMIN = ZP(M)
                    ELSE
                     ZP(M) = RDEF
                    ENDIF

   50            CONTINUE

                 AST = AST + ISN

                 IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                    WRITE(NUT, 960) ZMIN, ZMAX

                    WRITE(NUT, 800)
                    READ(NIN, *)ANS
                    IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
                      WRITE(NUT, 805)
                      READ(NIN, *)ANS
                      IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                THEN
                         WRITE(NUT,'("Input number of classes")')
                         READ(NIN, *)ICOL
                         CALL RCLSCA(ICOL, IPERC, 0)
                      ELSE IF(ANS(1:1) .EQ. "U" .OR. ANS(1:1) .EQ. "u")
     +                THEN
                         IUCOL = 1
                         WRITE(NUT,'("Input number of classes")')
                         READ(NIN, *)ICOL
                         IF(ICOL+1 .GT. MCOL) THEN
                            WRITE(NUT,963) ICOL, MCOL-1
                            ICOL = MCOL-1
                         ENDIF
                         DO 51 M =1, ICOL+1
                            WRITE(NUT, '("Input Limit",2X, I2)') M
                            READ(NIN,*) ZL(M)
   51                    CONTINUE
                         CALL RCLASS(ZL, ICOL+1, 0)
                      ELSE
                         WRITE(NUT, 810)
                         READ(NIN, *) ICOL, ZL(1), ZL(2)
                         CALL RCLASS(ZL, ICOL, -5)
                      ENDIF
                    ELSE
                      CALL RCLASS(ZL, ICOL, -5)
                    ENDIF

                    DO 52 M=1, ICOL

                       CLSTY(M) = 0
                       CWI(M) = CONLW
                       IF((ZL(1) + (M-1)*ZL(2)) .LT. 0.0) 
     +                   CLSTY(M) = CONLS

  52                CONTINUE
                    CALL GCONWI(CWI, ICOL)
                    CALL GCONDA(CLSTY, ICOL)
                    CALL RSHADE(ITAB,0)

                 ENDIF

                 CALL GIMETH(ITY)

                 CALL GINTP1(XP, YP, ZP, IBR, GRID, NPX, NPY)

                 CALL SUPPRES(GRID, TDGRD, XP, YP, ZP, SP, IBR, NPX,   
     +                       NPY, IFD, ITD, ISTY(IPR), IFPT, ITPT, NIN,  
     +                       NUT, ITOT, DCTY, RDEF)


               ELSE

                 IF(DCTY .LT. 3 .OR. DCTY .GT. 4) THEN

                    WRITE(NUT, 946)
                    DCTY = 3
C                    AST = AST + ISN
C                    GOTO 100

                 ENDIF

                 CALL GSYMDA(1)

                 IF(DCTY .EQ. 4) THEN     

                    WRITE(NUT, 948)
                    DO 60 M = 1, IS

                       IF(IFTY(M) .EQ. 1) 
     +                   WRITE(NUT, '(A30, 2X, I2)')
     +                      SLAB(ISTY(M))(1:SLLT(ISTY(M))), M

   60               CONTINUE

                    WRITE(NUT, 949)
                    READ(NIN,*)IB
                    IF(IB .LE. 0 .OR. IB .GT. M) THEN

                       WRITE(NUT, 951)
                       GOTO 100

                    ENDIF

                    WRITE(NUT, 952) SLAB(ISTY(IB))(1:SLLT(ISTY(IB)))

                    IADD = 0

                    DO 65 M = 1, IS

                       IF(IFTY(M) .GT. 1)THEN

                          IF(IB .GT. M) IADD = IADD + ISN

                       ENDIF

   65               CONTINUE

                    IV = (IB - 1) * ISN + IADD

                    FSCALE = FSS

                    IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                      WRITE(NUT, 780)
                      READ(NIN, *)ANS
                      IF(ANS(1:1) . EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                THEN
                        WRITE(NUT, 785)
                        READ(NIN, *)FSCALE
                        IF(FSCALE .LT. 0) FSCALE = ABS(FSCALE)
                      ENDIF

                    ENDIF

                    IFMNX = 0

                    DO 66 M = 1, IBR

                      IF(ABS(SP(IV + IAR + M) - RDEF) .GT. 1.0E-4)THEN
                        IF(IFMNX .EQ. 0) THEN
                          ZMIN = FSCALE * SP(IV + IAR + M)
                          ZMAX = ZMIN
                          IFMNX = 1
                        ENDIF

                        ZP(M) = FSCALE * SP(IV + IAR + M)
                        IF(ZP(M) .GT. ZMAX) ZMAX = ZP(M)
                        IF(ZP(M) .LT. ZMIN) ZMIN = ZP(M)
                      ELSE
                        ZP(M) = RDEF
                      ENDIF
                
   66               CONTINUE

                
                    IF(CONT .EQ. "F" .OR. CONT .EQ. "f") THEN

                       WRITE(NUT, 960) ZMIN, ZMAX

                       WRITE(NUT, 800)
                       READ(NIN, *)ANS
                       IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
                         WRITE(NUT, 805)
                         READ(NIN, *)ANS
                         IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") 
     +                   THEN
                            WRITE(NUT,'("Input number of classes")')
                            READ(NIN, *)ICOL
                            CALL RCLSCA(ICOL, IPERC, 0)
                         ELSE
                            WRITE(NUT, 810)
                            READ(NIN, *) ICOL, ZL(1), ZL(2)
                            CALL RCLASS(ZL, ICOL, -5)
                         ENDIF
                       ELSE
                         CALL RCLASS(ZL, ICOL, -5)
                       ENDIF

                       DO 67 M=1, ICOL

                          CLSTY(M) = 0
                          CWI(M) = CONLW
                          IF((ZL(1) + (M-1)*ZL(2)) .LT. 0.0) 
     +                      CLSTY(I) = CONLS

  67                   CONTINUE
                       CALL GCONWI(CWI, ICOL)
                       CALL GCONDA(CLSTY, ICOL)
                       CALL RSHADE(ITAB,0)

                    ENDIF


                    CALL GIMETH(ITY)
                    CALL GINTP1(XP, YP, ZP, IBR, GRID, NPX, NPY)

                    CALL SUPPRES(GRID, TDGRD, XP, YP, ZP, SP, IBR, NPX,   
     +                   NPY, IFD, ITD, ISTY(IB), IFPT, ITPT, NIN, NUT, 
     +                   ITOT, DCTY, RDEF)


                    WRITE(NUT, 964)
                    READ(NIN, *)IVLN

                    IF(IVLN .EQ. 1) THEN
                       CALL GCNR2S(GRID, NPX, NPY)
                    ELSE

                       CALL GCNR2V(GRID, NPX, NPY)

                       IF(ICOLL .EQ. 1) THEN

                         CALL RQSHAD(SDARRY, ICNW)
                         CALL GCONCO(SDARRY, ICNW)
                         CALL GCNR2V(GRID, NPX, NPY)

                       ENDIF

                    ENDIF

                    IF(IBBX .EQ. 0) THEN
                      CALL GREGDR(IDREG, -1, 0.0, 1, 0)
                    ELSE
                      CALL GREGDR(IDREG, -1, -1.0, 1, 0)
                    ENDIF

                    CALL GCLOPT(NCHAR, TT, 1.0*HEIGHT, NDEC, 0.0, 0)

                    IF(NCSC .EQ. 0) THEN
                      CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                    ELSE IF(NCSC .EQ. 1) THEN
                      CALL GCOSCS(XA2+(XA2-XA1)*0.01, 
     +                            YA1-(YA2-YA1)*0.04)
                    ELSE IF(NCSC .EQ. 2) THEN
                      CALL GCOSLI(XA2+(XA2-XA1)*0.01, 
     +                            YA1-(YA2-YA1)*0.04)
                    ENDIF      

                 ENDIF

C
C-----------------------------------------------------------------------
C
C Put vector atributes into array, grouped data 
C  i.e. p1.atr1, p1.atr2, p2.atr1, p2.atr2, .....
C

                 IV = AST + ISN
                 IL = 1

                 IVP = 0
                 MAXX = 0.

                 WRITE(NUT, 770)
                 READ(NIN, *) IVP

                 IF(IVP .EQ. 1) THEN
C                   WRITE(NUT, 775)
C                   READ(NIN, *)NPX, NPY
C                   IF(NVX .GT. NPX) NVX = NPX
C                   IF(NVY .GT. NPY) NVY = NPY


                   CALL GIMETH(ITY)
                   CALL GINTP1(XP, YP,SP(AST+IAR), IBR, GRDV, NVX, NVY)

                   DO 72 M=1, NVY
                      DO 73 L=1, NVX
                         ZP(IL) = GRDV(L, M) * GRDV(L, M)
                         ZP(IL+1) = GRDV(L,M)
                         IL = IL + 2

                         IF(((IL-1)/2)+1 .GT. ISN) GOTO 73 
   73                 CONTINUE
   72              CONTINUE
   74              CONTINUE

                   IL = 1

                   CALL GIMETH(ITY)
                   CALL GINTP1(XP, YP,SP(IV+IAR), IBR, GRDV, NVX, NVY)
                   DO 77 M=1, NVY
                      DO 78 L=1, NVX
                         IF(ZP(IL+1) .NE. RDEF)THEN
                         ZP(IL) = SQRT(ZP(IL) + GRDV(L, M) * GRDV(L, M)
     +                            *YEXAG*YEXAG)
                         ZP(IL+1) = ATAN2(GRDV(L,M)*YEXAG, ZP(IL+1))
                         ZP(IL+1) = ZP(IL+1) * PII
                         MAXX = MAX(MAXX, ZP(IL))

                         ELSE
                           ZP(IL) = RDEF
                         ENDIF

                         IF(ARMINS) THEN
                          IF(ZP(IL) .LT. SUPLEV) ZP(IL) = RDEF
                         ENDIF

                         IL = IL + 2
                         IF(((IL-1)/2)+1 .GT. ISN) THEN
                           WRITE(NUT, 776)
                           GOTO 79 
                         ENDIF
   78                 CONTINUE
   77              CONTINUE
   79              CONTINUE
                  
                 ELSE

                    DO 80 M = 1, IBR

                       ZP(IL) = SQRT(SP(AST+IAR+M) * SP(AST+IAR+M) + 
     +                          SP(IV+IAR+M) * SP(IV+IAR+M)*YEXAG*YEXAG)


                       ZP(IL+1) = ATAN2(SP(IV+IAR+M)*YEXAG, 
     +                            SP(AST+IAR+M))


                       ZP(IL+1) = ZP(IL+1) * PII


                       MAXX = MAX(MAXX, ZP(IL))

                       IF(ARMINS) THEN
                         IF(ZP(IL) .LT. SUPLEV) ZP(IL) = RDEF
                       ENDIF


                       IL = IL + 2

   80               CONTINUE

                 ENDIF

C
C--------------------------------------------------------------------
C
C Check for vectors out off bounds.
C
                 IL = 1
                 DO 81 M=1, IBR
                    IF(((XA2 - XP(M))*(XP(M) - XA1) .LT. 0) .OR.
     +                 ((YA2 - YP(M))*(YP(M) - YA1) .LT. 0)) THEN
                       ZP(IL) = RDEF
                       ZP(I+1) = RDEF
                    ENDIF
                    IL = IL + 2
   81            CONTINUE


C
C--------------------------------------------------------------------
C
C Scaling for flow field arrows.
C

                 WRITE(NUT, 996) VVSC
                 READ(NIN, *)ANS
                 IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
                    WRITE(NUT, 997)
                    READ(NIN,*)VVSC
                 ENDIF

                 XSC = MAXX / VVSC

                 YSC = XSC * YSIZE / (XSIZE * YEXAG)

                 CALL GSYMSY(2, 2)
                 CALL GSYMSC(1, 1, XSC) 
                 AST = AST + 2 * ISN
                 
               ENDIF

C
C-----------------------------------------------------------------------
C
C Perform statistic surpression 
C

               CALL GCLOPT(NCHAR, TT, 1.0*HEIGHT, NDEC, 0.0, 0)

               IF(IBBX .EQ. 0) THEN
                 CALL GREGDR(IDREG, -1, 0.0, 1, 0)
               ELSE
                 CALL GREGDR(IDREG, -1, -1.0, 1, 0)
               ENDIF

C
C-----------------------------------------------------------------------
C
C Perform contour and flow plotting
C

               IF(DCTY .EQ. 0) THEN

C                  CALL RFMODE('MIX')

                  WRITE(NUT, 961)
                  READ(NIN, *)ANS
                  IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
                     WRITE(NUT, 962)
                     READ(NIN, *) ITCOL
                     IF(ITCOL .GT. MCOL) THEN
                       WRITE(NUT,963) ITCOL, MCOL
                       ITCOL = MCOL
                     ENDIF
                     DO 86 M =1, ITCOL
                        WRITE(NUT, '("Input Limit",2X, I2)') M
                        READ(NIN,*) ZLT(M)
   86                CONTINUE

                     CALL RCLASS(ZLT, ITCOL, 0)

                     CALL RCSDEF(ITABB, BSCALE, NUMCOL, CMODE, COLMAX)

                     CALL RSHADE(ITABB, 0)

                     CALL GCNR2S(GRID, NPX, NPY)

                     IF(NCSC .EQ. 0) THEN
                       CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                     ELSE IF(NCSC .EQ. 1) THEN
                       CALL GCOSCS(XA2+(XA2-XA1)*0.01, 
     +                             YA1-(YA2-YA1)*0.04)
                     ENDIF

                  ENDIF

                  IF(IUCOL .EQ. 1) THEN
                     CALL RCLASS(ZL, ICOL+1, 0)
                  ELSE
                     CALL RCLASS(ZL, ICOL, -5)
                  ENDIF
                  CALL RSHADE(ITAB, 0)


                  WRITE(NUT, 849)
                  READ(NIN, *)INSTY
                  WRITE(NUT, 850)
                  DO 87 M=1, INSTY
                    READ(NIN, *) ILNSTY, CLSTY(ILNSTY), CWI(ILNSTY), 
     +                           CAHGT(ILNSTY)
  87              CONTINUE 
                  IF(INSTY .GT. 0) THEN
                    CALL GCONDA(CLSTY, ICOL)
                    CALL GCONWI(CWI, ICOL)
		    CALL GCONAH(CAHGT, ICOL)
		  ENDIF
  
                  SDARRY(1) = CNTCOL
                  CALL GCONCO(SDARRY(1), 1)
                  CALL GCNR2V(GRID, NPX, NPY)

                  IF(ICOLL .EQ. 1) THEN

                    CALL RQSHAD(SDARRY, ICNW)
                    CALL GCONCO(SDARRY, ICNW)
                    CALL GCNR2V(GRID, NPX, NPY)

                  ENDIF


                  IF(NCSC .EQ. 2) THEN
                    CALL GCOSLI(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ENDIF


               ELSE IF(DCTY .EQ. 1) THEN


                  CALL GCNR2S(GRID, NPX, NPY)

                  IF(NCSC .EQ. 0) THEN
                    CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                  ELSE IF(NCSC .EQ. 1) THEN
                    CALL GCOSCS(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ELSE IF(NCSC .EQ. 2) THEN
                    CALL GCOSLI(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ENDIF

               ELSE IF(DCTY .EQ. 2) THEN

                  CALL GCNR2S(GRID, NPX, NPY)

                  CALL GCNR2V(GRID, NPX, NPY)

                  
                  IF(NCSC .EQ. 0) THEN
                    CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                  ELSE IF(NCSC .EQ. 1) THEN
                    CALL GCOSCS(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ELSE IF(NCSC .EQ. 2) THEN
                    CALL GCOSLI(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ENDIF

               ELSE IF (DCTY .EQ. 3 .OR. DCTY .EQ. 4) THEN

                  IF(DCTY .EQ. 4) CALL GSYMCO(7, 33, 0)

                  IF(IVP .EQ. 1) THEN
                     CALL GSYMDR(ZP, 2, NVX, NVY)
                  ELSE
                     CALL GSYMDI(XP, YP, ZP, 2, IBR)
                  ENDIF

               ELSE IF (DCTY .EQ. 5 .OR. DCTY .EQ. 6) THEN

                  CALL RCLASS(ZL, ICOL, -5)
                  CALL RSHADE(ITAB, 0)

                  WRITE(NUT, 849)
                  READ(NIN, *)INSTY
		  
                  WRITE(NUT, 850)
                  DO 88 M=1, INSTY
                    READ(NIN, *) ILNSTY, CLSTY(ILNSTY), CWI(ILNSTY), 
     +                           CAHGT(ILNSTY) 
  88              CONTINUE 
                  IF(INSTY .GT. 0) THEN
                    CALL GCONDA(CLSTY, ICOL)
                    CALL GCONWI(CWI, ICOL)
		    CALL GCONAH(CAHGT, ICOL)
		  ENDIF

                  CALL GCNR2V(GRID, NPX, NPY)

                  IF(DCTY .EQ. 5) THEN

                     IF(ICOLL .EQ. 1) THEN

                       CALL RQSHAD(SDARRY, ICNW)
                       CALL GCONCO(SDARRY, ICNW)
                       CALL GCNR2V(GRID, NPX, NPY)

                     ENDIF


                     IF(NCSC .EQ. 2) THEN
                       CALL GCOSLI(XA2+(XA2-XA1)*0.01, 
     +                             YA1-(YA2-YA1)*0.04)
                     ENDIF

                  ENDIF
                  
               ELSE IF (DCTY .EQ. 7) THEN
C ----------------------------------------------------------------------
C  Suppress points outside region and below threshold

                  WRITE(NUT, 860)
                  READ(NIN, *)DTHR

                  DO 90 M=1, IBR
                     IF(((XA2 - XP(M))*(XP(M) - XA1) .LT. 0) .OR.
     +                  ((YA2 - YP(M))*(YP(M) - YA1) .LT. 0) .OR.
     +                  ZP(M) .LT. DTHR                          ) THEN
                        ZP(M) = RDEF
                     ENDIF
   90             CONTINUE 
	       
	       
                  CALL GPOI2Z(XP, YP, ZP, IBR)
		  
		  IF(NCSC .EQ. 0) THEN
                    CALL GCOSCL(XA2+(XA2-XA1)*0.01, YA1)
                  ELSE IF(NCSC .EQ. 1) THEN
                    CALL GCOSCS(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ELSE IF(NCSC .EQ. 2) THEN
                    CALL GCOSLI(XA2+(XA2-XA1)*0.01, YA1-(YA2-YA1)*0.04)
                  ENDIF
          
               ELSE

                  WRITE(NUT, 950) DCTY

               ENDIF


C
C-----------------------------------------------------------------------
C
C Plot country map
C
                  

               IF((DCTY .EQ. 0 .OR. DCTY .EQ. 3 .OR. DCTY .EQ. 5) 
     +          .AND. CMCOL .EQ. 0)THEN
                  CCOL = 1
               ELSE
                  CCOL = CMCOL
               ENDIF


               CALL COUNTRY(DCM, CMI, CMXG, CMYG, XA1, XA2, YA1, YA2,
     +                      CT, XS, CCOL)

               IF(IMASK) THEN

                  DO 95 M=1, IREG
                     IF(IRREG(M) .GT. 0 ) THEN
                        CALL GREGSS(IRREG(M), 1)
                        CALL GREGDR(IRREG(M), 1, 0.0, 1, 0)
                     ENDIF

   95             CONTINUE

               ENDIF

C
C-----------------------------------------------------------------------
C
C PLOT AND LABEL AXES
C
               IF(IAXTY .EQ. 0)THEN
                  CALL RAXDIS(6, 1, 0)
                  IF(NAL .EQ. 0)CALL RAXDIS(4, 0, 0)
                  IF(NOAX .EQ. 0) THEN
                     CALL RAXTEA(6, 2, 2, RDEF, IDEF, IDEF, RDEF)
                     CALL RAXTEX(6, -1, XLAB, RDEF, RDEF, 0.8*HEIGHT)
                     CALL RAXIS(1, YA1, HEIGHT, 1)
                     CALL RAXTEX(6, -1, YLAB, RDEF, RDEF, 0.8*HEIGHT)
                     CALL RAXIS(2, XA1, HEIGHT, 1)
                  ENDIF
C                ELSE IF(IAXTY .EQ. 1) THEN
C                  IF(NOAX .EQ. 0) THEN
C                  ENDIF
                ELSE
                  WRITE(NUT, 995)
                ENDIF

C
C-----------------------------------------------------------------------
C
C Plot labeling 
C

               CALL RTXCOL(2,2)
               CALL RTXHEI(1.5 * THS2 * HEIGHT)
               IF(NTEXT .GT. 0) THEN
                 IF(IB .LE. 0) THEN
                    CALL RTX(-2, SLAB(ISTY(IPR))(1:SLLN(ISTY(IPR))), 
     +                       XA1+(XA2-XA1)*0.1,   
     +                       YA2+0.08*(YA2-YA1))
                    CALL RTXHEI(1.5 * THS2* HEIGHT)
                    CALL RTX(-2, "kernel="// DENS(KERN(IPR))//"  ABW= ",
     +                       XA1+(XA2-XA1)*0.02,   
     +                       YA2+0.02*(YA2-YA1))
                    CALL RTXNC(SMTH(IPR), 2)
                    CALL RTXC(9, " degrees")
                 ELSE
                    CALL RTXHEI(1.5 * THS2* HEIGHT)                 
                    CALL RTX(-2,
     +                     SLAB(ISTY(IPR))(1:SLLN(ISTY(IPR)))//" and "
     +                     //SLAB(ISTY(IB))(1:SLLN(ISTY(IB))), 
     +                     XA1+(XA2-XA1)*0.0, YA2+0.02*(YA2-YA1))               
                 ENDIF
               ENDIF

  100          CONTINUE
               CALL RTXCOL(1, 1)
               CALL GSCAMM
               ISEG = ISEG + 1
               IPR = IPR + 1


   40       CONTINUE

   30    CONTINUE

   25    CONTINUE

C
C-----------------------------------------------------------------------
C
C Plot main title
C

         CALL RTXHEI(1.5 * THS1 * HEIGHT)
         CALL RTXCOL(2, 2)
         CALL RTX(-2, TITLE, 0.2*INPR*LPX, 1.06*INCL*LPY+YCEN)

         CALL GEMPTY

         WRITE(NUT, 970)
         READ(NIN, *)ANS
         IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
            WRITE(NUT, 980)
            READ(NIN, '(A100)') TITSUB
            CALL RTXHEI(1.5 * THS1 * HEIGHT)
            CALL RTX(-2, TITSUB, 0.2*INPR*LPX, 1.01*INCL*LPY+YCEN)
         ENDIF 

         CALL RTXCOL(1, 1)

         CALL GEMPTY

         CALL GSEGCL(K)
         CALL GCLEAR

         WRITE(NUT, 990)
         READ(NIN, *)

         K = K + 1

   20 CONTINUE

      CALL GCLOSE


      RETURN
  770 FORMAT('Do you want to interpolate the vectors before plotting,'/
     +       'input "0" for no or "1" for yes.                       ')
  775 FORMAT('How many grid points in X and Y are required for vectors')
  776 FORMAT('***WARNING***, too many points requested, results may be'/
     +       'ill-defined and unreliable.                             ')
  780 FORMAT('Do you want to re-scale the field before plotting.')
  785 FORMAT('What scaling do you want, must be positive.')
  790 FORMAT('There are currently ', I2, ' bins starting at ',
     +       F10.4, ' with step size ', F10.4)
  792 FORMAT('Do you want a unique output directory for the graphics?')
  795 FORMAT('Process Id is ', I8.8);
  800 FORMAT('do you wish to change the colour scale parameters')
  805 FORMAT('do you want automatic colour scale parameters, or user'/
     +       'chosen intervals, y, u or n')
  810 FORMAT('input the new colour scale parameters, CNUM, '
     +       'CINIT and CSTEP')
  815 FORMAT('Plot contour lines for background, input 1,            '/
     +       'shaded map without lines,          input 2,            '/
     +       'otherwize shaded map with lines.   input 3,            ')
  835 FORMAT('Input a suppresion threshold for this statistic.       ')
  849 FORMAT('How many orientation contour levels?                   ')
  850 FORMAT('Specify contour level IDs, styles, widths and          '/
     +       'annotation heights.                                    ')
  860 FORMAT('What plotting threshold is required?                   ')
  900 FORMAT('how many plots per segment are required')
  910 FORMAT('how many plots per row are required')
  920 FORMAT('SET DEFAULT CONTOURING TYPE')
  930 FORMAT('what is the contouring type required?                  '/
     +       '  input "0" for isoline map                            '/
     +       '  input "1" for colour-shaded map                      '/
     +       '  input "2" for composite map (shaded and isoline),    '/
     +       '            same statsitic.                            '/
     +       '  input "3" for simple flow map                        '/
     +       '  input "4" for flow map with colour contour background'/
     +       '  input "5" for isoline map with grey scale background,'/
     +       '            different statistics.                      '/
     +       '  input "6" for isoline map with colour background,    '/
     +       '            different statistics.                      '/
     +       '  input "7" for stippling.                             ')
  940 FORMAT('do you want free control over contouring, type "f",   '/
     +       '                  or use default contour, type "d".   ')
  941 FORMAT('CREATING COUNTRY MASK, PLEASE WAIT ................   ')
  945 FORMAT('***WARNING***, vector plots are not possible with '/
     +       '               scaler data. Using default.')
  946 FORMAT('***WARNING***, scaler plots are not possible with '/
     +       '               vector data, except as background. '/
     +       '               Using default.                   ')
  948 FORMAT('AVAILABLE SCALER FILDS ARE:-')
  949 FORMAT('input relvant integer descriptor')
  950 FORMAT('***ERROR***, no contour type available for value=', I2)
  951 FORMAT('***ERROR***, no data for colour coding')
  952 FORMAT('Plotting background field:- ', A30)
  953 FORMAT('Plotting foreground field:- ', A30)
  960 FORMAT('the min. and max. values of the field are', 2(F10.5,2X))
  961 FORMAT('Do you want contour shading?')
  962 FORMAT('Input the number of class limits, n, and their values')
  963 FORMAT('Maximum number of classes allowed for shading is', I2/
     +       'resetting number of classes limits to ', I2)
  964 FORMAT('What contouring type is required for composite vector'/
     +       'plots, Line=0, Shaded=1')
  970 FORMAT('do you want to use a subsidary title to the default')
  980 FORMAT('input the subsidary title')
  990 FORMAT('press RETURN to continue')
  995 FORMAT('***WARNING***, Axes not yet supported for this type of'/
     +       '               projection, no axes plotted.           ')
  996 FORMAT('The current scaling for vector plots is ', F8.4, ' do'/
     +       'you wish to use a different scaling.                  ')
  997 FORMAT('Input the new scaling factor.')
  998 FORMAT('The countour levels can now be changed for the line'/
     +       'contours.')
      END

C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C                                                                      C
C                         SUPPORT ROUTINES                             C
C                                                                      C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C

      SUBROUTINE COUNTRY(DCM, CMI, CMXG, CMYG, XA1, XA2, YA1, YA2, CTOL, 
     +                   XS, CCOL)

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                      C
C Routine to plot a country map.                                       C
C                                                                      C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C

      INTEGER DCM, CMI(DCM)
      INTEGER ICM, IC, CCOL

      REAL CMXG(DCM), CMYG(DCM), XA1, XA2, YA1, YA2, CTOL, XS
      REAL PTX, PTY

      LOGICAL SWTCH1, SWTCH2, SWTCH

      ICM = 1

   40 CONTINUE
      IF(ICM+1 .LE. DCM) THEN 
       PTX = (CMXG(ICM) - XA1)*(XA2 - CMXG(ICM))
       PTY = (CMYG(ICM) - YA1)*(YA2 - CMYG(ICM))
       IF((PTX .GE. 0.0) .AND. (PTY .GE. 0.0)) THEN
         SWTCH1 = .TRUE.
       ELSE
         SWTCH1 = .FALSE.
       ENDIF
   50  CONTINUE 
         
       IC = ICM + 1
       IF(IC .GT. DCM) GOTO 60
       PTX = (CMXG(IC) - XA1)*(XA2 - CMXG(IC))
       PTY = (CMYG(IC) - YA1)*(YA2 - CMYG(IC))
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

          IF(SWTCH)THEN
             CALL GWICOL(XS, CCOL)
             CALL GVECT(CMXG(ICM), CMYG(ICM), 0)
             CALL GVECT(CMXG(IC), CMYG(IC), 1)
          ENDIF

          ICM = ICM+1
          SWTCH1 = SWTCH2
          GOTO 50
       ELSE
          ICM = ICM+1
          GOTO 40
        ENDIF
      ENDIF
   60 CONTINUE
      RETURN
      END



      SUBROUTINE MASK(DCM, CMI, CMXG, CMYG, XA1, XA2, YA1, YA2, CTOL, 
     +                XS, IDREG, IRN, IREG)

      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                      C
C Routine to plot a country mask.                                      C
C                                                                      C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C

      INTEGER DCM, CMI(DCM)
      INTEGER ICM, IC, IREG, IRN, IDREG(IRN), INREG, I, J, IR
      INTEGER IRR, IBP, ICR, IRREG

      REAL CMXG(DCM), CMYG(DCM), XA1, XA2, YA1, YA2, CTOL, XS
      REAL XT, YT, XT1, YT1


      ICM = 1
      IREG = 0
      INREG = 0

      CALL GREGLO(XA1, YA1, -1, IRREG)
      CALL GREGLO(XA2, YA1, -1, IRREG)
      CALL GREGLO(XA2, YA2, -1, IRREG)
      CALL GREGLO(XA1, YA2, -1, IRREG)
      CALL GREGLO(XA1, YA1, -1, IRREG)

   40 CONTINUE
      IF(ICM+1 .LE. DCM) THEN

       IC = ICM + 1

       IF(CMI(ICM) .EQ. CMI(IC) ) THEN
           INREG = INREG + 1

           IF(INREG .EQ. 1) THEN

              IREG = IREG + 1
              IRR = IREG
              ICR = IREG
       
              CALL GREGLO(CMXG(ICM), CMYG(ICM), 1, IDREG(ICR))
              XT = CMXG(ICM)
              YT = CMYG(ICM)
              IBP = IRR


           ELSE
              CALL GREGLO(CMXG(ICM), CMYG(ICM), -1, IDREG(ICR))
           ENDIF

           IF(ABS(CMXG(ICM)-CMXG(IC)) .GT. CTOL) THEN

             IF(IBP .EQ. IRR) THEN
                IREG = IREG + 1
                IBP = IREG
                ICR = IREG
                XT1 = CMXG(IC)
                YT1 = CMYG(IC)

             ELSE
                IBP = IRR
                CALL GREGLO(XT1, YT1, -1, IDREG(ICR))
                ICR = IRR
             ENDIF

           ENDIF

       ELSE
          CALL GREGLO(CMXG(ICM), CMYG(ICM), -1, IDREG(ICR))
          CALL GREGLO(XT, YT, -1, IDREG(ICR))
          INREG = 0
       ENDIF

       ICM = IC
       GOTO 40
      ENDIF

      CALL GREGLO(XT, YT, -1, IDREG(ICR))



      DO 100 I=1, IREG

         DO 150 J=1, IREG

            IF(I .NE. J) THEN

               IR = 0

               IF(IDREG(I) .GT. 0 .AND. IDREG(J) .GT. 0)
     +            CALL GREGOP(IDREG(I), IDREG(J), 2, IR)

               IF(IR .GT. 0) THEN
                  CALL GREGOP(IDREG(I), IR, 3, IDREG(I))
                  CALL GREGOP(IDREG(J), IR, 3, IDREG(J))
               ENDIF


            ENDIF

  150    CONTINUE

  100 CONTINUE

      DO 160 I=1, IREG

         IF(IDREG(I) .GT. 0)
     +      CALL GREGOP(IDREG(I), IRREG, 2, IDREG(I))

  160 CONTINUE


      RETURN
      END


      SUBROUTINE SUPPRES(GRID, TDGRD, XP, YP, ZP, SP, IBR, NPX, NPY,  
     +                   IFD, ITD, ISTY, IFPT, ITPT, NIN, NUT, ITOT,
     +                   DCTY, RDEF)
      IMPLICIT NONE

C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                                      C
C Statistic supression routine.                                        C
C                                                                      C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C

      INTEGER L, M
      INTEGER NPX, NPY
      INTEGER IBR, IFD, ITD, ISTY
      INTEGER ITDD
      INTEGER NIN, NUT
      INTEGER IFPT, ITPT
      INTEGER ITOT, DCTY

      REAL GRID(NPX, NPY), TDGRD(NPX, NPY)
      REAL XP(IBR), YP(IBR), ZP(IBR)
      REAL SP(ITOT)
      REAL DSUP, RDEF

      CHARACTER ANS*3

      ITDD = 0

C      IF(ISTY .GE. 5 .AND. ISTY .LE. 8) RETURN

      IF(IFD .EQ. 1 .OR. ITD .EQ. 1) THEN

        WRITE(NUT, 815)
        READ(NIN, *)ANS
        IF(ANS(1:1) .EQ. "Y" .OR. ANS(1:1) .EQ. "y") THEN
           IF(IFD .eq. 1 .AND. ITD .EQ. 1) THEN
              WRITE(NUT, 816)
              READ(NIN, *)ITDD
              IF(ITDD .LT. 1 .OR. ITDD .GT. 2) THEN
                 WRITE(NUT, 817)ITDD
              ENDIF
           ELSE IF(IFD .EQ. 1) THEN
              ITDD = 1
           ELSE IF(ITD .EQ. 1) THEN
              ITDD = 2
           ELSE
           ENDIF 

           WRITE(NUT,820)
           READ(NIN,*) DSUP
             
           IF(ITDD .EQ. 1) THEN

              IF(DCTY .EQ. 7) THEN
                DO 70 M=1, IBR
                   IF(SP(IFPT + M) .LT. DSUP) ZP(M) = RDEF
   70           CONTINUE
              ELSE
                DO 80 M=1, IBR
                   ZP(M) = SP(IFPT + M)
   80           CONTINUE

                CALL GINTP1(XP, YP, ZP, IBR, TDGRD, NPX, NPY)

                DO 82 M = 1, NPX
                   DO 83 L = 1, NPY
                     IF(TDGRD(M,L) .LT. DSUP)GRID(M, L) = RDEF
   83              CONTINUE
   82           CONTINUE 
              ENDIF

           ELSE IF(ITDD .EQ. 2) THEN  

              IF(DCTY .EQ. 7) THEN
                DO 85 M=1, IBR
                   IF(SP(ITPT + M) .LT. DSUP) ZP(M) = RDEF
   85           CONTINUE
              ELSE
                DO 90 M=1, IBR
                   ZP(M) = SP(ITPT + M)
   90           CONTINUE

                CALL GINTP1(XP, YP, ZP, IBR, TDGRD, NPX, NPY)               
                DO 95 M = 1, NPX
                   DO 96 L = 1, NPY
                     IF(TDGRD(M,L) .LT. DSUP)GRID(M, L) = RDEF
   96              CONTINUE
   95           CONTINUE 
              ENDIF
           ENDIF

        ENDIF

      ENDIF

      RETURN

  815 FORMAT('statistic suppresion is available do you want to       '/
     +       'suppress the statistic in a low data density region.   '/
     +       'Yes or No                                              ')
  816 FORMAT('Do you want FEATURE, "1" or TRACK, "2" density         '/
     +       'suppression                                            ')
  817 FORMAT('***WARNING***, no suppression available for value=', I2 )
  820 FORMAT('What suppresion threshold is required.                 ')

      END


      SUBROUTINE COLDEF(IDC)

      INTEGER IDC

      INCLUDE 'paramc.ptr'
      INCLUDE 'paramcol.ptr'

#ifdef    NAMELISTS
      READ(NAMCHN, PARAMCOL)
      REWIND(NAMCHN)
#endif

      CALL RCSDEF(IDC,CSCALE,NCOL,CMODE,100.0)

      RETURN
      END

#endif
