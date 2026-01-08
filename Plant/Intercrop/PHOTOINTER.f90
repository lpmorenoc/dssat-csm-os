!C=======================================================================
!C  PHOTOINTER, Subroutine
!C  Compute daily intercepted PAR for intercrops based on Kropff et al.,
!C  1992 method and the PHOTO (canopy method available in DSSAT)
!C-----------------------------------------------------------------------
!C  REVISION HISTORY
!! 07/06/2022 Created from the PHOTO subroutine
!!-----------------------------------------------------------------------
!!  Called from:   INTERCROP
!C=======================================================================

      SUBROUTINE PHOTOINTER(CONTROL,                      &
         CROPS, KCANM, XLAIM, CANHTM, PAR,                & !Input
         FracIntRadM )                                      !Output

!-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      USE ModuleData
      IMPLICIT NONE
      SAVE

      CHARACTER*30 FILEIO

      INTEGER DYNAMIC
      INTEGER Canopylayers, I, J, Lcount

      REAL PAR, CANHTMAX, LAITemp,BaseHeightLayer, TopHeightLayer

      
      REAL, DIMENSION(2)  :: PREFCM, CANHTM, KCANM, XLAIM, IPARM, FracIntRadM
      REAL, DIMENSION(20,2) :: IPARML, LAIByLayer, CumLAIByLayer
      REAL, DIMENSION(20)   ::KLAIProd, KCumLAIProd
      CHARACTER*2, CROPS(2)

!-----------------------------------------------------------------------
!     Define constructed variable types based on definitions in
!     ModuleDefs.for.

!     The variable "CONTROL" is of type "ControlType".
      TYPE (ControlType) CONTROL

!     Transfer values from constructed data types into local variables.
      DYNAMIC = CONTROL % DYNAMIC
      FILEIO  = CONTROL % FILEIO

!C***********************************************************************
!C***********************************************************************
!C     Run Initialization - Called once per simulation
!C***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!C-----------------------------------------------------------------------
!      CALL PHOTIP(FILEIO,  
!     &  CCEFF, CCMAX, CCMP, FNPGN, FNPGT, LMXSTD, LNREF, PARMAX,   
!     &  PGREF, PHTHRS10, PHTMAX, ROWSPC, TYPPGN, TYPPGT, XPGSLW, YPGSLW)
!
!C-----------------------------------------------------------------------
!C     Adjust canopy photosynthesis for GENETIC input value of
!C     maximum leaf photosyntheses (LMXSTD).  Exponential curve from
!C     from the hedgerow photosynthesis model (Boote et al, 199?).
!C-----------------------------------------------------------------------
!      IF (PGREF .GT. 0.) THEN
!        PGLFMX = (1. - EXP(-1.6 * LMXSTD)) / (1. - EXP(-1.6 * PGREF))
!      ELSE
!        PGLFMX = 1.0
!      ENDIF

!     LPM reflection coefficient of canopy (equal both crops to test)   
      PREFCM = 0.0125

!C***********************************************************************
!C***********************************************************************
!C     Seasonal initialization - run once per season
!C***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!C-----------------------------------------------------------------------
      CANHTMAX = 0.0
      IPARM = 0.0

!C***********************************************************************
!C***********************************************************************
!C     Daily rate calculations
!C***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!C-----------------------------------------------------------------------
!C     Calculate intercepted radiation by 10 cm layers?
!C-----------------------------------------------------------------------
      !CANHTM(2) = CANHTM(1) !To test same height
      CANHTMAX = MAX(0.0,MAXVAL(CANHTM))
      Canopylayers = MIN(20, MAX(1,(INT(CANHTMAX * 10.0)+1)))
      LAIByLayer = 0.0
      KLAIProd = 0.0
      KCumLAIProd = 0.0
      CumLAIByLayer = 0.0
      DO I=1,2
          LAITemp = 0.0
          BaseHeightLayer = 0.0 
          TopHeightLayer = 0.0
          DO Lcount=0, Canopylayers-1
              J= Canopylayers - Lcount
              BaseHeightLayer = MAX(0.0 ,((J * 0.1) - 0.1))
              TopHeightLayer = MAX(0.0 ,(J * 0.1))
               IF (BaseHeightLayer < CANHTM(I) .AND. CANHTM(I) > 0.0) THEN
                   IF (TopHeightLayer > CANHTM(I)) THEN
                       LAIByLayer(J,I) = (XLAIM(I)/ CANHTM(I)) * MAX(0.0,(CANHTM(I)-BaseHeightLayer))
                   ELSE
                       LAIByLayer(J,I) = (XLAIM(I)/ CANHTM(I)) * 0.1
                   ENDIF
                   CumLAIByLayer(J,I) = LAITemp + LAIByLayer(J,I)
                   LAITemp = CumLAIByLayer(J,I)
               ENDIF
              KCumLAIProd(J) = KCumLAIProd(J) +(KCANM(I) * CumLAIByLayer(J,I))
              KLAIProd(J) = KLAIProd(J) +(KCANM(I) * LAIByLayer(J,I))
          ENDDO  
      ENDDO
      
      IPARM = 0.0 
      IPARML = 0.0
      DO I=1,2
          BaseHeightLayer = 0.0 
          DO Lcount=0, Canopylayers-1
              J= Canopylayers - Lcount
              BaseHeightLayer = MAX(0.0 ,((J * 0.1) - 0.1))
               IF (BaseHeightLayer < CANHTM(I) .AND. CANHTM(I) > 0.0 .AND. KLAIProd(J) > 0.0) THEN
                   IF (J == Canopylayers) THEN
                       IPARML(J,I) = (KCANM(I) * LAIByLayer(J,I) / KLAIProd(J)) * (1.0 - PREFCM(I)) * PAR * (1 - EXP(-KLAIProd(J)))
                   ELSE
                       IPARML(J,I) = (KCANM(I) * LAIByLayer(J,I) / KLAIProd(J)) * (1.0 - PREFCM(I)) * PAR * (1 - EXP(-KLAIProd(J))) * EXP(-KCumLAIProd(J+1))
                   ENDIF
                   
                   IPARM(I) = IPARM(I) + IPARML(J,I)
               ENDIF
          ENDDO
          IF (PAR > 0.0) FracIntRadM(I) = MIN(1.0, IPARM(I)/PAR)
      ENDDO
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END !SUBROUTINE PHOTOINTER

!=======================================================================
! Variable definitions for PHOTOINTER
!=======================================================================
! ERR      Error code for file operation 
! ERRKEY   Subroutine name for error file 
! FILEC    Filename for SPE file (e.g., SBGRO980.SPE) 
! FILECC   Path plus filename for species file (*.spe) 
! FILEIO   Filename for input file (e.g., IBSNAT35.INP) 
! KCAN     Canopy light extinction coefficient for daily PAR, for 
!            equidistant plant spacing, modified when in-row and between 
!            row spacing are not equal 
! KCANR    Canopy light extinction coefficient, reduced for incomplete 
!            canopy 
! LMXSTD   Maximum leaf photosyntheses for standard cultivar 
! LNREF    Value of leaf N above which canopy PG is maximum (for standard 
!            cultivar) 
! LNUM     Current line number of input file 
! LUNCRP   Logical unit number for FILEC (*.spe file) 
! LUNIO    Logical unit number for FILEIO 
! PAR      Daily photosynthetically active radiation or photon flux density
!            (moles[quanta]/m2-d)

!=======================================================================
! END SUBROUTINE PHOTO
!=======================================================================
