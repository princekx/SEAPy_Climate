C
C-----------------------------------------------------------------------
C
C Data statment of labels for plotting of statistical data.
C
C
C NLAB == no. of statistics labels
C SLAB == array of labels
C


      INTEGER NLAB, DTY

      PARAMETER(NLAB=17, DTY=6)

      INTEGER SLLN(NLAB), SLLT(NLAB)

      CHARACTER*40 SLAB(NLAB)
      CHARACTER*15 DENS(DTY)

#ifdef NAMELISTS

      NAMELIST /PARAMS_LAB/ SLLN, SLLT, SLAB, DENS

#else


      DATA SLLN / 14,30,10,27,15,15,13,13,13,14,17,10,11,8,9,6,6/
      DATA SLLT / 13,30,10,27,15,15,13,13,13,13,17,10,11,8,9,6,6/

      DATA SLAB / "Mean Strength", "Standard Deviation of Strength",
     +            "Mean Speed", "Standard Deviation of Speed",
     +            "Feature Density", "Genesis Density", "Lysis Density",
     +            "Track Density", "Mean Velocity", "Mean Lifetime",
     +            "Growth/Decay Rate", "Anisotropy", "Orientation" ,
     +            "Tendency", "Mean Area", "Spare1", "Spare2" /

      DATA DENS / "Fisher (exp)", "Constant", "Linear", "Quadratic",
     +            "Non-Isotropic", "Power (SQT)" /


#endif
