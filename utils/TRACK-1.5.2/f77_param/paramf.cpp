C
C Parameters required for plotting processed data.
C
C FILNMF == the file for plotted, thresholded data and feature points.
C FILFEAT == file for plotted tracks for a single run.
C SPLICE == file for plotted tracks for a combined set of runs.
C SCALEX, SCALEY == scaling for plotting purposes.
C NDIN == input channel for reading input data to add to track plots
C TITLE1, TITLE2, TITLE3== titles of plots.
C DOVLAY == switch for plotting nodes.
C NU == switch for re-gridding non-uniform nodes
C
C
      INCLUDE '../f77_param/path.ptr'

      INTEGER MCTE
      PARAMETER(MCTE=20)

      CHARACTER*70 FILNAM, FILNMF, FILFET, SPLICE
      CHARACTER*50 TITLE1, TITLE2, TITLE3

      INTEGER NU, DOVLAY
      INTEGER INCT
      INTEGER IFR
      INTEGER NVEC, NPNT
      INTEGER IACT, NTC
      INTEGER ARRHTP
      INTEGER CONLS

      LOGICAL IZ

      REAL SCALEX, SCALEY, SOFFX, SOFFY, SSC
      REAL SUPVEC, VECSCL 
      REAL ARRSCL, PTSCL
      REAL ZNCT(MCTE)
      REAL CHS
      REAL CONLW, CSMTH

#ifdef NAMELISTS

      NAMELIST /PARAMF/ FILNAM, FILNMF, FILFET, SPLICE, TITLE1, TITLE2,
     +                  TITLE3, DOVLAY, NU, IZ, SCALEX, SCALEY, SOFFX,
     +                  SOFFY, SSC, SUPVEC, VECSCL, ARRSCL, ARRHTP, 
     +                  PTSCL, INCT, ZNCT, CHS, IFR, NVEC, NPNT, IACT, 
     +                  NTC, CONLW, CONLS, CSMTH

#else


      PARAMETER(FILNAM=PATH//'threshpl.pic')
      PARAMETER(FILNMF=PATH//'filterpl.pic')
      PARAMETER(FILFET=PATH//'featpl.pic')
      PARAMETER(SPLICE=PATH//'splicepl.pic')

      PARAMETER(TITLE1 = 'RELATIVE VORTICITY AT 850mb')

C      PARAMETER(TITLE1 ='')

      PARAMETER(TITLE2 = 
     + 'TRACKS OF RELATIVE VORTICITY AT 850mb')

      PARAMETER(TITLE3 = 
     +          'RELATIVE VORTICITY TRACKS AT 850mb')

C      PARAMETER(TITLE3 ='')

      PARAMETER(DOVLAY = 0, NU=1, IZ=.TRUE.)
C
C SCALEX=2.6, SCALEY=1.4 suitable scaling for output for Postscript
C
      PARAMETER(SCALEX=1.0, SCALEY=-1.0)
      PARAMETER(SOFFX=0.15, SOFFY=0.4)
      PARAMETER(SSC=2.0)
      PARAMETER(SUPVEC=0.001, VECSCL=0.05)
      PARAMETER(ARRSCL=0.001, PTSCL=0.003)
      PARAMETER(INCT=3)
      PARAMETER(CONLW=0.0, CONLS=7, CSMTH=5.0)

C ARRHTP controls arror head type, 1 for solid
      PARAMETER(ARRHTP=1)
C
C NVEC controls plotting of arrow heads
C NPNT controls plotting of feature points
C IACT controls using different color tables
C NTC  controls track color, default is 27.

      PARAMETER(IFR=0, CHS=1.0, NVEC=1, NPNT=1, IACT=1, NTC=27)

      DATA ZNCT /240.0, 268.0, 320.0/

#endif

      REAL CWI(1)

      INTEGER CLSTY(1)
