C=======================================================================
C  Denit_DayCent, Subroutine
C
C  Determines denitrification based on DayCent model

C-----------------------------------------------------------------------
C  Revision history
C  06/12/2014 PG / CHP Written
C-----------------------------------------------------------------------
C  Called : SOIL
C  Calls  : Fert_Place, IPSOIL, NCHECK, NFLUX, RPLACE,
C           SOILNI, YR_DOY, FLOOD_CHEM, OXLAYER
C=======================================================================

      SUBROUTINE Denit_DayCent (DYNAMIC, ISWNIT, NSWITCH, 
     &    BD, DUL, KG2PPM, NO3, SAT, ST, SW     !Input
     &    DENITRIF, DLTSNO3)                        !Output

!-----------------------------------------------------------------------
      USE ModuleDefs 
      USE ModuleData
      IMPLICIT  NONE
      SAVE
!-----------------------------------------------------------------------
      CHARACTER*1 ISWNIT

!      LOGICAL IUON

      INTEGER DOY, DYNAMIC, L
      INTEGER NLAYR
      !INTEGER NSOURCE, YEAR, YRDOY   

      REAL ARNTRF, DLAG, FLOOD, TNOX, TNOXD, XMIN
      REAL TFDENIT, WFDENIT
      REAL DENITRIF(NL), NITRIF(NL), SNO3_AVAIL

      REAL DLTSNO3(NL)   
      REAL BD(NL), DUL(NL)
      REAL KG2PPM(NL) 
      REAL NO3(NL), SAT(NL)
      REAL SNO3(NL), SW(NL)
      
      INTEGER NSWITCH

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!     ------------------------------------------------------------------
!       Today's values
!       Seasonal cumulative vaules
        TNOX   = 0.0    !denitrification
        TN2O   = 0.0    ! N2O added        PG

!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS 
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!     ------------------------------------------------------------------
      IF (INDEX('N',ISWNIT) > 0) RETURN

!     ------------------------------------------------------------------
!     Loop through soil layers for rate calculations
!     ------------------------------------------------------------------
      TNOXD = 0.0
      TN2OD = 0.0     ! PG

!-----------------------------------------------------------------------
!       Denitrification section
!-----------------------------------------------------------------------
!       Denitrification only occurs if there is nitrate, SW > DUL and
!       soil temperature > 5.
        IF (NO3(L) .GT. 0.01 .AND. SW(L) .GT. DUL(L) .AND.
     &       ST(L) .GE. 5.0) THEN

!         Water extractable soil carbon: estimated according to
!         Rolston et al. 1980, as cited in Godwin & Jones 1991 (ASA
!         monograph #31). Calculate carbohydrate carbon as 40% of the
!         carbohydrate pool.
C-UPS     Corrected per e-mail 03/29/00
!         CW = 24.5 + 0.0031 * (SSOMC(L) + 0.4 * FPOOL(L,1)) * KG2PPM(L)

!     ----------------------------------------------------------------
!11/18/2003 UPS: THE NEW NTRANS SHOULD READ: 
!         CW = 24.5 + {0.0031 * SSOMC(L) + 0.4 * FPOOL(L,1)} * KG2PPM(L) 
!            = 24.5 + AVAILABLE CARBON FROM HUMIC FRACTION + FRESH C from  CARBOHYDRATE POOL 

!NOTES: 1. ONLY THE HUMIC C IS MULTIPLIED BY 0.0031 
!       2. SOILC IN GODWIN&JONES INCLUDED BOTH HUMIC C AND FRESH (LITTER POOL C) 
!       3. WE ARE USING ONLY THE CARBOHYDRATE POOL (*0.4 TO C) = ALL AVAILABLE 
!       4. FPOOL is still kg of Organic matter/ha (and not kg C/ha)?? 

!     SO WE NEED TO FIX BOTH DSSAT4 AND GODWIN AND SINGH EQN TO THE ABOVE. 

!          CW = 24.5 + (0.0031 * SSOMC(L) + 0.4 * FPOOL(L,1)) * KG2PPM(L) 

!     The above removed on 1/14/2004 as per email from AJG and UPS

!     ----------------------------------------------------------------
!     DENITRIFICATION - CORRECTIONS - 13 Jan 2004 (US)
!         From NTRANS:
!          CW = 24.5 + 0.0031 * (HUMC(L) + 0.4 * FPOOL(L,1)) * KG2PPM(L)
!     ----------------------------------------------------------------

!        From Century:
         !CHP changed 1/14/2004 per email from UPS / AJG
          CW = 24.5 + 0.0031 * (SSOMC(L) + 0.2 * LITC(L)) * KG2PPM(L)
!     ----------------------------------------------------------------
!
!     The DENITRIF or DNRATE calculations are identical in NTRANS 
!             (DSSAT4) and Godwin and Singh:
!
!     DENITRIF = {6.0 * 1.E-04 * CW * NO3(L) * WFDENIT * TFDENIT * DLAG }/ KG2PPM(L)
!     -- AS in NTRANS (DSSAT4)
!
!         = {in concentration unit}/KG2PPM
!         = kg N (denitrified)
!         = {in conc unit }/(10/BD*DLAYR)
!
!     DENITRIF = 6.0*1.E-05 * CW * NO3(L) * WFDENIT * TFDENIT * DLAG*BD(L)*DLAYR(L)
!     -- AS in GODWIN & SINGH
!
!     NOTE: CW is in concentration unit!!  Extractable C concentration.
!
!     CW = (SOILC*KG2PPM(L))*0.0031 + 24.5
!
!       where SOILC = SSOMC(L) + 0.4 * FOM(L) in kg/ha - Origianl definition
!        Later corrected to:
!        SOILC = SSOMC(L) + 0.4 * FPOOL(L,1)  -- because only carbohydrate 
!             pool from the FOM is assumed to provide labile/extractable C.
! 
!     CW   = ({SSOMC(L) + 0.4 * FPOOL(L,1)} * KG2PPM(L))*0.0031 + 24.5
!
!     The equations in Godwin and Jones as well as Godwin and Singh are incorrect!!

!     ----------------------------------------------------------------
!         Temperature factor for denitrification.
          TFDENIT = 0.1 * EXP (0.046 * ST(L))
          TFDENIT = AMAX1 (AMIN1 (TFDENIT, 1.), 0.)

!         Water factor for denitrification: only if SW > DUL.
          WFDENIT = 1. - (SAT(L) - SW(L)) / (SAT(L) - DUL(L))
          WFDENIT = AMAX1 (AMIN1 (WFDENIT, 1.), 0.)

          IF (WFDENIT .GT. 0.0) THEN
            DLAG(L) = DLAG(L) + 1
          ELSE
            DLAG(L) = 0
          ENDIF

          IF (DLAG(L) .LT. 5) THEN
            WFDENIT = 0.0
          ENDIF

!         Denitrification rate
C-UPS     Corrected per e-mail 03/29/00
!         DLAG REMOVED REVISED-US 4/20/2004
          DENITRIF(L) = 6.0 * 1.E-04 * CW * NO3(L) * WFDENIT * 
     &                 TFDENIT / KG2PPM(L)       
          DENITRIF(L) = AMAX1 (DENITRIF(L), 0.0)

!         The minimum amount of NO3 that stays behind in the soil and 
!         cannot denitrify is XMIN.
!         XMIN    = 0.25 / KG2PPM(L)
          XMIN    = 0.       !AJG

!         Check that no more NO3 denitrifies than there is, taking
!         into account what has already been removed by other
!         processes (thus use only negative DLTSNO3 values). This is a
!         protection against negative values at the integration step.
          SNO3_AVAIL = SNO3(L) + AMIN1 (DLTSNO3(L), 0.) - XMIN

!         Take the minimum of the calculated denitrification and the
!         amount of NO3 available for denitrification. 
          DENITRIF(L)  = AMIN1 (DENITRIF(L), SNO3_AVAIL)

C         If flooded, lose all nitrate --------REVISED-US
!          IF (FLOOD .GT. 0.0) THEN
!            !DNFRATE = SNO3(L) - 0.5/KG2PPM(L)        !XMIN?, SNO3_AVAIL?
!            DNFRATE = SNO3_AVAIL - 0.5/KG2PPM(L)        !XMIN?, SNO3_AVAIL?
!          ELSE
!            DNFRATE = 0.0
!          ENDIF

!         chp/us 4/21/2006
          IF (FLOOD .GT. 0.0 .AND. WFDENIT > 0.0) THEN
!            DENITRIF(L) = SNO3_AVAIL
!           chp 9/6/2011 remove 50% NO3/d = 97% removed in 5 days
!           previously removed 100% NO3/d
            DENITRIF(L) = SNO3_AVAIL * 0.5
          ENDIF

!chp 4/20/2004   DENITRIF = AMAX1 (DENITRIF, DNFRATE)
          DENITRIF(L) = AMAX1 (DENITRIF(L), 0.0)
          IF (NSWITCH .EQ. 6) THEN
            DENITRIF(L) = 0.0
          ENDIF

!         Reduce soil NO3 by the amount denitrified and add this to
!         the NOx pool
          DLTSNO3(L) = DLTSNO3(L) - DENITRIF(L)
          TNOX       = TNOX       + DENITRIF(L)
          TNOXD      = TNOXD      + DENITRIF(L)

        ELSE
!         IF SW, ST OR NO3 FALL BELOW CRITICAL IN ANY LAYER RESET LAG EFFECT.
          DLAG(L) = 0      !REVISED-US
        ENDIF   !End of IF block on denitrification.

      END DO   !End of soil layer loop.

      END DO   !End of soil layer loop.

      CALL PUT('NITR','TNOXD',ARNTRF) 
      CALL PUT('NITR','TN2OD',TN2OD)

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF

C-----------------------------------------------------------------------
  
      RETURN
      END SUBROUTINE Denit_DayCent

!=======================================================================
! Denit_DayCent Variables 
!-----------------------------------------------------------------------

!***********************************************************************
