*        SUBROUTINE 'GPFA'
*        SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
*
*        *** THIS IS THE ALL-FORTRAN VERSION ***
*            -------------------------------
*
*        CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
*
*        A IS FIRST REAL INPUT/OUTPUT VECTOR
*        B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
*        TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
*              BY CALLING SUBROUTINE 'SETGPFA'
*        INC IS THE INCREMENT WITHIN EACH DATA VECTOR
*        JUMP IS THE INCREMENT BETWEEN DATA VECTORS
*        N IS THE LENGTH OF THE TRANSFORMS:
*          -----------------------------------
*            N = (2**IP) * (3**IQ) * (5**IR)
*          -----------------------------------
*        LOT IS THE NUMBER OF TRANSFORMS
*        ISIGN = +1 FOR FORWARD TRANSFORM
*              = -1 FOR INVERSE TRANSFORM
*
*        WRITTEN BY CLIVE TEMPERTON
*        RECHERCHE EN PREVISION NUMERIQUE
*        ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
*
*----------------------------------------------------------------------
*
*        DEFINITION OF TRANSFORM
*        -----------------------
*
*        X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
*
*---------------------------------------------------------------------
*
*        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
*        SEE:
*
*        C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
*          FOR ANY N = (2**P)(3**Q)(5**R)",
*          SIAM J. SCI. STAT. COMP., MAY 1992.
*
*----------------------------------------------------------------------
*
      SUBROUTINE GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN,NJ)
      implicit double precision (a-h, o-z)
      implicit integer(i-n)
*
      DIMENSION A(*), B(*), TRIGS(*)
      DIMENSION NJ(3)

      IP = NJ(1)
      IQ = NJ(2)
      IR = NJ(3)
*
*     COMPUTE THE TRANSFORM
*     ---------------------
      I = 1
      IF (NJ(1).GT.0) THEN
         CALL GPFA2F(A,B,TRIGS,INC,JUMP,N,IP,LOT,ISIGN)
         I = I + 2 * ( 2**IP)
      ENDIF
      IF (NJ(2).GT.0) THEN
         CALL GPFA3F(A,B,TRIGS(I),INC,JUMP,N,IQ,LOT,ISIGN)
         I = I + 2 * (3**IQ)
      ENDIF
      IF (NJ(3).GT.0) THEN
         CALL GPFA5F(A,B,TRIGS(I),INC,JUMP,N,IR,LOT,ISIGN)
      ENDIF
*
      RETURN
      END
