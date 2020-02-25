C
C Parameters for color table definition for plots of statistics.
C

      INTEGER CSDIM, CMODE
      INTEGER NCOL

      PARAMETER (CSDIM=10)

      REAL CSCALE(4, CSDIM)
      REAL MAXLEV

#ifdef  NAMELISTS

      NAMELIST /PARAMCOL/ CMODE, CSCALE, NCOL, MAXLEV

#else

      PARAMETER(CMODE=1, MAXLEV=100.0, NCOL=4)

      DATA CSCALE /1.0, 100.0, 100.0,  0.0,
     .             1.0, 100.0,   0.0, 100.0,
     .             1.0,   0.0, 100.0, 100.0
     .             1.0,   0.0,   0.0,  100.0/

#endif
