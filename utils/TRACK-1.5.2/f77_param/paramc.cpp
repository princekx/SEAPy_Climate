C
C Parameters common to plotting routines
C
C ZINIT, STEP, COL, ICNUM == UNIRAS parameters for colour plotting.
C XLAB, YLAB == labels for axes.
C SCALE == Suitable scaling of the data for plotting.
C IX, IY == full data dimensions.
C NOUT == output channel for writing the plotting file.
C SCALE == Suitable scaling of the data for plotting.
C VALUE == data value for undefined data point.
C NIN, NUT == standard I/O.
C NAL == switch for axis labels
C
      CHARACTER*30 XLAB, YLAB, DEVICE
      CHARACTER*20 DUMMY

      INTEGER COL, ICNUM, IDEF
      INTEGER NIN, NOUT, NUT, NAL, NOAX, NCSC
      INTEGER IAXTY, ICONT
      INTEGER ITAB
      INTEGER NAMCHN

      REAL ZINIT, STEP, VALUE, SCALE, RDEF 
      REAL TOLB

      PARAMETER(NAMCHN=20)

C
C Device for UNIRAS
C
C mx11      -- X display driver.
C hcposta4  -- A4 postcript driver.
C ldummy    -- dummy driver.
C

#ifdef NAMELISTS

      NAMELIST /PARAMC/ DEVICE, DUMMY, ZINIT, STEP, VALUE, COL, ICNUM,
     +                  RDEF, IDEF, XLAB, YLAB, TOLB, NOUT, SCALE, NIN,
     +                  NUT, NAL, NOAX, NCSC, IAXTY, ICONT, ITAB

#else

      PARAMETER(DEVICE='sel mx11;e')
      PARAMETER(DUMMY=' sel ldummy;e')
      PARAMETER(ZINIT=0.0, STEP=0.2, VALUE = -1.0, COL=32, ICNUM=20)
      PARAMETER(RDEF = 999.999, IDEF = 9999)
      PARAMETER(XLAB='LONGITUDE', YLAB='LATITUDE')

      PARAMETER(TOLB=5.0)
      PARAMETER(NOUT=10, SCALE=1.)
      PARAMETER(NIN=5, NUT=6)
      PARAMETER(NAL=1, NOAX=1, NCSC=1)
      PARAMETER(IAXTY=0, ICONT=1)
      PARAMETER(ITAB=1)

#endif
