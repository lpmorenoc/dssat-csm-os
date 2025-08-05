C=======================================================================
C COPYRIGHT 1998-2022
C                     DSSAT Foundation
C                     University of Florida, Gainesville, Florida
C                     International Fertilizer Development Center
C                     
C ALL RIGHTS RESERVED
C=======================================================================
C=======================================================================
C  INTERCROP, Subroutine
C
C  This routine calls intercropping growth routines.
C-----------------------------------------------------------------------
C  Revision history
C
!  06/24/2022 LPM Create subroutine based on PLANT
C=======================================================================

      SUBROUTINE INTERCROP(CONTROL, ISWITCH,
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    IRRAMT, NH4, NO3, SKi_Avail, SPi_AVAIL,         !Input
     &    SNOW, SOILPROP, SRFTEMP, ST, SW,                !Input
     &    TRWUP, WEATHER, YREND, YRPLT,                   !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output

C-----------------------------------------------------------------------
!     The following models are currently supported:
!         'CRGRO' - CROPGRO
!         'MZCER' - CERES-Maize
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
! Each plant module must compute SATFAC, SWFAC, and TURFAC
C-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
      USE FloodModule

      IMPLICIT NONE
      EXTERNAL WARNING, CROPGRO,CSCERES_Interface, MZ_CERES,
     &  PHOTOINTER, NUPTAKINTER
      SAVE

      CHARACTER*1  MEEVP, RNMODE, ISWNIT, ISWWAT
      CHARACTER*6  ERRKEY
      PARAMETER (ERRKEY = 'INTERCROP')
      CHARACTER*8  MODEL
      CHARACTER*78 MESSAGE(10)    !Up to 10 lines of text to be output

      INTEGER DYNAMIC
      INTEGER RUN !, NVALP0
      INTEGER YREND, MDATE, YRPLT  !, YRSIM, YREMRG
      INTEGER STGDOY(20)

      REAL CANHT, CO2, DAYL, EO, EOP, EORATIO, EOS, EP, ES
      REAL KCAN, KEP, KSEVAP, KTRANS, LAI, NSTRES
      REAL PORMIN, RWUEP1, RWUMX, SRFTEMP, SNOW, IRRAMT
      REAL TMAX, TMIN, TRWU
      REAL TRWUP, TWILEN, XLAI, XHLAI
      
      REAL, DIMENSION(2)  :: HARVFRAC
      REAL, DIMENSION(NL) :: NH4, NO3, RLV, UPPM  !, RWU
      REAL, DIMENSION(NL) :: ST, SW, UNO3, UNH4, UH2O

      LOGICAL FixCanht, BUNDED    !, CRGRO
c-----------------------------------------------------------------------
C         Variables needed to run ceres maize.....W.D.B. 12-20-01
      CHARACTER*2 CROP
      REAL    SRAD

!-----------------------------------------------------------------------
C         Variables to run CASUPRO from Alt_PLANT.  FSR 07-23-03
      REAL PAR, TAVG, TGROAV  !CHP 7/26/04 , TDAY
      REAL TGRO(TS)

	INTEGER, PARAMETER :: CanopyLayers=3
	REAL, DIMENSION(1:NumOfStalks,CanopyLayers) :: LFmntDEF
!     P model
      REAL, DIMENSION(NL) :: PUptake, SPi_AVAIL, FracRts
      REAL PSTRES1

!     K model
      REAL, DIMENSION(NL) :: KUptake, SKi_Avail

!     ORYZA Rice model
      REAL, DIMENSION(0:NL) :: SomLitC
      REAL, DIMENSION(0:NL,NELEM) :: SomLitE
      LOGICAL, PARAMETER :: OR_OUTPUT = .FALSE.
      
!     Variables intercrop      
      INTEGER, DIMENSION(NumOfCrops)  :: MDATEM
      INTEGER STGDOYM(20,NumOfCrops)
      INTEGER      I, J, L
      REAL TotIntRad, EORATIOC
      REAL RTNO3, RTNH4
      REAL, DIMENSION(NL) :: RWU, KG2PPM
      REAL, DIMENSION(NumOfCrops)  :: CANHTM, EORATIOM, KCANM, KEPM
      REAL, DIMENSION(NumOfCrops)  :: KSEVAPM, KTRANSM, EOPM
      REAL, DIMENSION(NumOfCrops)  :: NSTRESM, PORMINM, PSTRES1M, RWUMXM
      REAL, DIMENSION(NumOfCrops)  :: XLAIM, XHLAIM, IPARM, FracIntRadM
      REAL, DIMENSION(NumOfCrops)  :: TRWUPM, FracIntRadMTrans
      REAL, DIMENSION(NumOfCrops)  :: TRNUM
      
      REAL, DIMENSION(NL,NumOfCrops) :: PUptakeM, FracRtsM, RLVM, UNO3M
      REAL, DIMENSION(NL,NumOfCrops) ::  UNH4M, KUptakeM, RWUM
      REAL, DIMENSION(NL,NumOfCrops) ::  PUNH4M, PUNO3M
      CHARACTER*2, CROPS(NumOfCrops)
      CHARACTER*8  MODELS(NumOfCrops)
      CHARACTER*30 FILEIOM(NumOfCrops)
      INTEGER :: test, DAS, CropStatus
!-----------------------------------------------------------------------
!     Constructed variables are defined in ModuleDefs.
      TYPE (ControlType)  CONTROL
      TYPE (SwitchType)   ISWITCH
      TYPE (SoilType)     SOILPROP
      TYPE (ResidueType)  HARVRES
      TYPE (ResidueType)  HARVRESM(NumOfCrops)
      TYPE (ResidueType)  SENESCE
      TYPE (ResidueType)  SENESCEM(NumOfCrops)
      TYPE (FloodWatType) FLOODWAT
      TYPE (FloodNType)   FLOODN
      TYPE (WeatherType)  WEATHER

!     Transfer values from constructed data types into local variables.
      CROP    = CONTROL % CROP
      DYNAMIC = CONTROL % DYNAMIC
      MODEL   = CONTROL % MODEL
      RNMODE  = CONTROL % RNMODE
      DAS     = CONTROL % DAS
      RUN     = CONTROL % RUN
      MODELS(1) = MODEL
      MODELS(2) = 'CRGRO'
      CROPS(1) = 'RY'
      CROPS(2) = 'CV'
      CONTROL % INTERCROP = CROPS
      FILEIOM(1) = 'DSSAT48.INP'
      FILEIOM(2) = 'DSSAT48_CV.INP'
      !LPM 09/14/2022 Added this temporary to allow N fixation
      ISWITCH % ISWSYM = 'Y'

      MEEVP  = ISWITCH % MEEVP
      ISWNIT = ISWITCH % ISWNIT
      ISWWAT = ISWITCH % ISWWAT
      BUNDED = FLOODWAT % BUNDED
      CO2    = WEATHER % CO2
      DAYL   = WEATHER % DAYL
      PAR    = WEATHER % PAR
      SRAD   = WEATHER % SRAD
      TAVG   = WEATHER % TAVG
      TGRO   = WEATHER % TGRO
      TGROAV = WEATHER % TGROAV
      TMAX   = WEATHER % TMAX
      TMIN   = WEATHER % TMIN
      TWILEN = WEATHER % TWILEN

!***********************************************************************
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
!     Non-CROPGRO crops can not use MEPHO = 'L' or MEEVP = 'Z' at
!       this time.  If Species files are modified for these options
!       in the future, we need to make this check on a crop by crop basis.
!     The plant routines do not use these codes, but the SPAM module
!       does and it will bomb when species parameters are not found.
      IF (INDEX(MODEL,'CRGRO') <= 0 .and. index(model,'PRFRM') <= 0
     &  .AND. ISWITCH % MEPHO .EQ. 'L') THEN
        ISWITCH % MEPHO = 'C'
!       Put ISWITCH data where it can be retreived
!         by other modules as needed.
        CALL PUT(ISWITCH)

!       Write message to WARNING.OUT file
        WRITE(MESSAGE(1),110)
        WRITE(MESSAGE(2),120) CROP
        WRITE(MESSAGE(3),130)
        CALL WARNING(3, ERRKEY, MESSAGE)
      ENDIF

  110 FORMAT('You have specified use of the Leaf-level photosynthesis')
  120 FORMAT('option, which is not available for crop ', A2, '.')
  130 FORMAT('Canopy photosynthesis option will be used.')

      IF (INDEX(MODEL,'CRGRO') <= 0 .and. index(model,'PRFRM') <= 0
     &  .AND. ISWITCH % MEEVP .EQ. 'Z') THEN
!       Default to Priestly-Taylor potential evapotranspiration
        ISWITCH % MEEVP = 'R'
!       Put ISWITCH data where it can be retreived
!         by other modules as needed.
        CALL PUT(ISWITCH)

!       Write message to WARNING.OUT file
        WRITE(MESSAGE(1),210)
        WRITE(MESSAGE(2),220) CROP
        WRITE(MESSAGE(3),230)
        WRITE(MESSAGE(4),240)
        CALL WARNING(4, ERRKEY, MESSAGE)
      ENDIF

  210 FORMAT('You have specified use of the Zonal evapotranspiration')
  220 FORMAT('option, which is not available for crop ', A2, '.')
  230 FORMAT('The Priestly-Taylor potential evapo-transpiration ')
  240 FORMAT('option will be used.')

!---------------------------------------------------------------------
!     Print warning if "dynamic ET" routine is used - need canopy height
!     Need to modify this crop code list as canopy height routines are added.
      IF (MEEVP .EQ. 'D' .AND. INDEX('RIWHMLMZSGPTBA',CROP) .GT. 0) THEN
        WRITE(MESSAGE(1),310)
        WRITE(MESSAGE(2),320)
        WRITE(MESSAGE(3),330) CROP
        WRITE(MESSAGE(4),340)
        CALL WARNING(4, ERRKEY, MESSAGE)
        !Trigger to set canopy height upon emergence
        FixCanht = .TRUE.
      ENDIF

  310 FORMAT('You are using the Penman-Monteith potential evapo-')
  320 FORMAT('transpiration method, which requires canopy height.')
  330 FORMAT('The ',A2,' crop routine does not calculate this ')
  340 FORMAT('value.  A default value will be used.')

!     Initialize output variables.
!     Each plant routine may or may not re-compute these values.
      CANHT    = 0.0
      CANHTM    = 0.0
      EOPM = 0.0
      EORATIO  = 1.0
      EORATIOC  = 0.0
      EORATIOM  = 1.0
      FracIntRadM = 0.0 
      FracIntRadMTrans = 0.0 
      KCAN     = 0.85
      KCANM     = 0.85
      KEP      = 1.0
      KEPM      = 1.0
      KSEVAP   = -99.
      KSEVAPM   = -99.
      KTRANS   = 1.0
      KTRANSM   = 1.0
      MDATE    = -99
      MDATEM    = -99
      NSTRES   = 1.0
      NSTRESM   = 1.0
      !NVALP0   = 10000
      PORMIN   = 0.02
      PORMINM   = 0.02
      PUNH4M    = 0.0
      PUNO3M    = 0.0
      RLV      = 0.0
      RLVM     = 0.0
      RTNH4   = 0.0
      RTNO3    = 0.0 
      RWUEP1   = 1.5
      RWUMX    = 0.03
      RWUMXM    = 0.03
      !SENESCE % CumResWt= 0.0
      !SENESCE % CumResE = 0.0
      STGDOY   = 9999999
      STGDOYM   = 9999999
      XHLAI    = 0.0
      XHLAIM    = 0.0
      XLAI     = 0.0
      XLAIM     = 0.0
      !YREMRG   = -99
      SENESCE % ResWt  = 0.0
      SENESCE % ResLig = 0.0
      SENESCE % ResE   = 0.0
      TotIntRad = 0.0 
      UNH4     = 0.0
      UNH4M    = 0.0
      UNO3     = 0.0
      UNO3M    = 0.0
      UH2O     = 0.0

      CALL READ_ASCE_KT(CONTROL, MEEVP)

!***********************************************************************
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
!     If this is not a sequenced run, don't use any previously calculated
!       harvest residue.
!     This should be done by each plant routine, but in case not:
      IF (RUN .EQ. 1 .OR. INDEX('QF',RNMODE) .LE. 0) THEN
        HARVRES % RESWT  = 0.0
        HARVRES % RESLig = 0.0
        HARVRES % RESE   = 0.0
        DO I=1, NumOfCrops
          HARVRESM(I) % RESWT  = 0.0
          HARVRESM(I) % RESLig = 0.0
          HARVRESM(I) % RESE   = 0.0
        ENDDO
!        HARVRES % CumResWt= 0.0
!        HARVRES % CumResE = 0.0
      ENDIF

!     Initialize output variables.
!     Each plant routine may or may not re-compute these values.
      CANHT    = 0.0
      CANHTM    = 0.0
      EOPM = 0.0 
!      EORATIO  = 1.0
      FracRts  = 0.0
      FracRtsM  = 0.0
      FracIntRadM = 0.0
      FracIntRadMTrans = 0.0 
!      KCAN     = 0.85
!      KEP      = 1.0
!      KSEVAP   = -99.
!      KTRANS   = 1.0
      KUptake = 0.0
      KUptakeM = 0.0
      NSTRES   = 1.0
      NSTRESM   = 1.0
!      PORMIN   = 0.02
      PUNH4M    = 0.0
      PUNO3M    = 0.0
      PSTRES1  = 1.0
      PSTRES1M  = 1.0
      PUPTAKE  = 0.0
      PUPTAKEM  = 0.0
      RLV      = 0.0
      RLVM     = 0.0
!      RWUEP1   = 1.5
!      RWUMX    = 0.03
      TotIntRad = 0.0
      UH2O     = 0.0
      UNH4     = 0.0
      UNH4M    = 0.0
      UNO3     = 0.0
      UNO3M    = 0.0
      XHLAI    = 0.0
      XHLAIM    = 0.0
      XLAI     = 0.0
      XLAIM     = 0.0
      SENESCE % ResWt  = 0.0
      SENESCE % ResLig = 0.0
      SENESCE % ResE   = 0.0
      DO I=1, NumOfCrops
        SENESCEM(I) % RESWT  = 0.0
        SENESCEM(I) % RESLig = 0.0
        SENESCEM(I) % RESE   = 0.0
      ENDDO

!***********************************************************************
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
      IF (CROP .NE. 'FA' .AND. 
     &    CONTROL % YRDOY .GE. YRPLT .AND. YRPLT .NE. -99) THEN

        SENESCE % ResWt  = 0.0
        SENESCE % ResLig = 0.0
        SENESCE % ResE   = 0.0
      DO I=1, NumOfCrops
        SENESCEM(I) % RESWT  = 0.0
        SENESCEM(I) % RESLig = 0.0
        SENESCEM(I) % RESE   = 0.0
      ENDDO

      ELSE
        CANHT = 0.0
        CANHTM = 0.0
        RLV   = 0.0
        RLVM   = 0.0
        XHLAI = 0.0
        XHLAIM = 0.0
        XLAI  = 0.0
        XLAIM = 0.0
        RETURN
      ENDIF

!-----------------------------------------------------------------------
      ENDIF   !End of dynamic loop prior to calls to crop models
!-----------------------------------------------------------------------

!***********************************************************************
!***********************************************************************
!     Call subroutine to define intercepted radiation by each crop:
       CALL PHOTOINTER(CONTROL,                      
     &    CROPS, KCANM, XLAIM, CANHTM, PAR,                 !Input
     &    FracIntRadM )                                     !Output
       !Estimate RWU by species
       
       CALL GET('SPAM','UH2O',RWU)
       Call GET('SPAM', 'EOPM',  EOPM, NumOfCrops)
       IF (ISWWAT .EQ. 'Y') THEN
           TRWUPM = 0.0
           DO I=1, NumOfCrops
               DO J=1, NL
                   IF (RLV(J) > 0.0) THEN
                       RWUM(J,I) = RWU(J) * RLVM(J,I)/RLV(J)
                   ELSE
                       RWUM(J,I) = 0.0
                   ENDIF
                   TRWUPM(I) = TRWUPM(I) + RWUM(J,I)
               ENDDO
           ENDDO
       ENDIF
       
      IF (ISWNIT .EQ. 'Y') THEN
          CALL NUPTAKINTER(DYNAMIC, SOILPROP,
!     &  NDMSDR, NDMTOT,                                  !Input
     &    NH4, NO3, RLV,RLVM, SW, !RTNH4M, RTNO3M,          !Input
     &   TRNUM, PUNH4M, PUNO3M)                            !Output
      ENDIF
      IF (DYNAMIC == INTEGR) THEN
!        OPEN (UNIT = test,FILE = 'pot_Nuptake.txt',POSITION="APPEND")
!         write (test, '(1I,44F8.4)') DAS, RLV(1:6), RLVM(1:6,:), 
!     &          PUNH4M(1:6,:), PUNO3M(1:6,:), TRNUM 
!        CLOSE (UNIT=test)
      ENDIF
            
      XLAI = 0.0
      XHLAI = 0.0
      RLV = 0.0
      HARVRES % RESWT  = 0.0
      HARVRES % RESLig = 0.0
      HARVRES % RESE   = 0.0
      SENESCE % ResWt  = 0.0
      SENESCE % ResLig = 0.0
      SENESCE % ResE   = 0.0
      EORATIOC = 0.0
      KUptake = 0.0
      UNH4 = 0.0
      UNO3 = 0.0 

!     Call crop models for all values of DYNAMIC:         
      DO I=1, NumOfCrops
        CONTROL % CROP = CROPS(I)
        CONTROL % FILEIO = FILEIOM(I)
        CROP = CONTROL % CROP
        Call PUT('PLANT', 'FracIntRadM',  FracIntRadM(I))
        IF (ISWNIT .EQ. 'Y') THEN
            Call PUT('SPAM', 'PUNO3M',  PUNO3M(:,I))
            Call PUT('SPAM', 'PUNH4M',  PUNH4M(:,I))
            Call PUT('SPAM', 'TRNUM',  TRNUM(I))
        ENDIF
        
        
      SELECT CASE (MODELS(I)(1:5))
!-----------------------------------------------------------------------
!     CROPGRO model
      CASE('CRGRO')
        CALL CROPGRO(CONTROL, ISWITCH,
     &    EOPM(I), HARVFRAC, NH4, NO3, SOILPROP, SPi_AVAIL,          !Input
     &    ST, SW, TRWUPM(I), WEATHER, YREND, YRPLT,                  !Input
     &    CANHTM(I), CropStatus,EORATIOM(I), HARVRESM(I), KSEVAPM(I),!Output
     &    KTRANSM(I), MDATEM(I),NSTRESM(I), PSTRES1M(I),             !Output
     &    PUptakeM(:,I), PORMINM(I), RLVM(:,I), RWUMXM(I),           !Output
     &    SENESCEM(I), STGDOYM(:,I), FracRtsM(:,I), UNH4M(:,I),      !Output
     &    UNO3M(:,I), XHLAIM(I), XLAIM(I))                           !Output

!     -------------------------------------------------
!     Maize, Sweetcorn
      CASE('MZCER')
        CALL MZ_CERES (CONTROL, ISWITCH,                         !Input
     &     EOPM(I), HARVFRAC, NH4, NO3, SKi_Avail,               !Input
     &     SPi_AVAIL, SNOW,                                      !Input
     &     SOILPROP, SW, TRWUPM(I), WEATHER, YREND, YRPLT,       !Input
     &     CropStatus,                                           !Output
     &     CANHTM(I), HARVRESM(I), KCANM(I), KEPM(I),            !Output
     &     KUptakeM(:,I), MDATEM(I), NSTRESM(I), PORMINM(I),     !Output
     &     PUptakeM(:,I),RLVM(:,I), RWUMXM(I), SENESCEM(I),      !Output
     &     STGDOYM(:,I),FracRtsM(:,I), UNH4M(:,I), UNO3M(:,I),   !Output
     &     XLAIM(I), XHLAIM(I))                                  !Output

        IF (DYNAMIC < RATE) THEN
!          KTRANS = KCAN + 0.15        !Or use KEP here??
          KTRANSM(I) = KEPM(I)        !KJB/WDB/CHP 10/22/2003
          KSEVAPM(I) = KEPM(I)
        ENDIF

!     -------------------------------------------------
!     Wheat and Barley CSCER
      CASE('CSCER')
        CALL CSCERES_Interface (CONTROL, ISWITCH,                !Input
     &     EOPM(I), YREND, NH4, NO3, SNOW, SOILPROP,             !Input
     &     SRFTEMP, ST, SW, TRWUPM(I), WEATHER, YRPLT, HARVFRAC, !Input
     &     CANHTM(I), HARVRESM(I), KCANM(I), KEPM(I), MDATEM(I), !Output
     &     NSTRESM(I),PORMINM(I), RLVM(:,I), RWUMXM(I),          !Output
     &     SENESCEM(I), STGDOYM(:,I), UNH4M(:,I), UNO3M(:,I),    !Output
     &     XLAIM(I))                                             !Output

        IF (DYNAMIC .EQ. SEASINIT) THEN
          KTRANSM(I) = KEPM(I)
          KSEVAPM(I) = KEPM(I)
          XHLAIM(I) = XLAIM(I)
        ELSEIF (DYNAMIC .EQ. INTEGR) THEN
          XHLAIM(I) = XLAIM(I)
        ENDIF
        
      END SELECT
       XLAI = XLAI + XLAIM(I) 
       XHLAI = XHLAI + XHLAIM(I)
       RLV =  RLV + RLVM(:,I)
       KUptake = KUptake + KUptakeM(:,I)
       UNH4 = UNH4 + UNH4M(:,I)
       UNO3 = UNO3 + UNO3M(:,I)
       EORATIOC = (EORATIOC + EORATIOM(I) * XHLAIM(I)) 
       HARVRES % RESWT  = HARVRES % RESWT + HARVRESM(I) % RESWT
       HARVRES % RESLig = HARVRES % RESLig + HARVRESM(I) % RESLig
       HARVRES % RESE   = HARVRES % RESE + HARVRESM(I) % RESE
       SENESCE % ResWt  = SENESCE % ResWt +  SENESCEM(I) % ResWt 
       SENESCE % ResLig = SENESCE % ResLig + SENESCEM(I) % ResLig
       SENESCE % ResE   = SENESCE % ResE  +  SENESCEM(I) % ResE  

      END DO
      Call PUT('PLANT', 'KTRANSM',  KTRANSM, NumOfCrops)
      Call PUT('PLANT', 'XHLAIM',  XHLAIM, NumOfCrops)
      !Call PUT('PLANT', 'RTNH4M', RTNH4M, NumOfCrops)
      !Call PUT('PLANT', 'RTNO3M', RTNO3M, NumOfCrops)
      !Estimate intercepted total radiation by crop to define 
      !transpiration. Use ktrans instead of kcan
       CALL PHOTOINTER(CONTROL,                      
     &    CROPS, KTRANSM, XHLAIM, CANHTM, SRAD,             !Input
     &    FracIntRadMTrans )                                !Output
      
       TotIntRad = SUM(FracIntRadMTrans)
       Call PUT('PLANT', 'FracIntRadMTrans', 
     &      FracIntRadMTrans, NumOfCrops)
      
      !Use/transfer maximum RWUMX to estimate total 
      !root water uptake (TRWUP)
      RWUMX = MAXVAL(RWUMXM)
      !Transfer maximum CANHT for LAND/SPAM
      CANHT = MAXVAL(CANHTM)
      !LPM 07/29/2022 Transfer maximum values of KTRANS, KSEVAP for LAND/SPAM. 
      !KTRANS is not used for intercropping in TRANS_Inter. 
      !KSEVAP is defined as -99 for CROPGRO and it is used in subroutine PSE
      !Should we change to the average or weighted average (by LAI) per species?
      KTRANS = MAXVAL(KTRANSM)
      KSEVAP = MAXVAL(KSEVAPM)
      IF (XHLAI> 0.0) EORATIO = EORATIOC/XHLAI
      
      !LPM 08/04/2022 Keep maximum value for N stress
      NSTRES = MAXVAL(NSTRESM)
      !LPM 08/04/2022 Keep maximum value for PORMIN
      PORMIN = MAXVAL(PORMINM)
      
      !LPM 08/05/2020 Modify KUptake in case is greater than SKi_Avail
      ! for both crops
      !LPM 08/30/2022 K uptake needs to be modified as done with N
      DO J = 1, NL
         IF (KUptake(J) > SKi_Avail(J)) THEN
          DO I=1, NumOfCrops
             KUptakeM(J,I) = KUptakeM(J,I) * SKi_Avail(J)/ KUptake(J)
          ENDDO
          KUptake(J) = SKi_Avail(J)
         ENDIF
      ENDDO
      
!      OPEN (UNIT = test,FILE = 'water_uptake.txt',POSITION="APPEND")
!       write (test, '(1I,8F8.3)') DAS, (TRWUPM/10.), (TRWUP/10.), 
!     & TRWU, (EP/10.), (EOP/10.), (EOPM/10.) 
!           CLOSE (UNIT=test)



!***********************************************************************
!***********************************************************************
!     Processing after calls to crop models:
!-----------------------------------------------------------------------
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
! Zero the value of HARVRES composite variable here
!!!NOTE: At this time, the variable has already been used to
!     initialize soil properties for this season.
!  This should be done by each plant routine, but in case not:
        HARVRES % RESWT  = 0.0
        HARVRES % RESLig = 0.0
        HARVRES % RESE   = 0.0

!***********************************************************************
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
!     Set default canopy height upon emergence (or first day with
!       LAI.  Should actually set these defaults within each
!       crop routine.
        IF (FixCanht .AND. (XLAI .GT. 0.0 .OR. XHLAI .GT. 0.0)) THEN
          CANHT = 0.5
          CANHTM = 0.5
          FixCanht = .FALSE.
        ENDIF

!***********************************************************************
      ENDIF

!***********************************************************************
      RETURN
      END SUBROUTINE INTERCROP

!===========================================================================
! Variable listing for Alt_Plant - updated 08/18/2003
! --------------------------------------------------------------------------
! CANHT     Canopy height (m)
! CO2       Atmospheric carbon dioxide concentration (µmol[CO2] / mol[air])
! CONTROL   Composite variable containing variables related to control
!             and/or timing of simulation.  The structure of the variable
!             (ControlType) is defined in ModuleDefs.for.
! CROP      Crop identification code
! DAYL      Day length on day of simulation (from sunrise to sunset) (hr)
! EOP       Potential plant transpiration rate (mm/d)
! EORATIO   Ratio of increase in EO with increase in LAI (up to LAI=6.0)
!             for use with FAO-56 Penman reference EO.
! ERRKEY    Subroutine name for error file
! FIXCANHT  Logical variable, =TRUE if default canopy height is to be set
!             by Alt_Plant routine upon emergence
! FLOODN    Composite variable which contains flood nitrogen mass and
!             concentrations. Structure of variable is defined in
!             ModuleDefs.for. (var.)
! FLOODWAT  Composite variable containing information related to bund
!             management. Structure of variable is defined in ModuleDefs.for.
! HARVFRAC  Two-element array containing fractions of (1) yield harvested
!             and (2) by-product harvested (fraction)
! HARVRES   Composite variable containing harvest residue amounts for total 
!             dry matter, lignin, and N amounts.  Structure of variable is 
!             defined in ModuleDefs.for. 
! ISWITCH   Composite variable containing switches which control flow of 
!             execution for model.  The structure of the variable 
!             (SwitchType) is defined in ModuleDefs.for. 
! IRRAMT    Irrigation amount (mm)
! KCAN      Canopy light extinction coefficient for daily PAR, for 
!             equidistant plant spacing, modified when in-row and between 
!             row spacing are not equal 
! KEP       Energy extinction coefficient for partitioning EO to EP 
! KSEVAP    Light extinction coefficient used for computation of soil 
!             evaporation 
! KTRANS    Light extinction coefficient used for computation of plant 
!             transpiration 
! MDATE     Harvest maturity date (YYYYDDD)
! MEEVP     Method of evapotranspiration ('P'=Penman, 'R'=Priestly-Taylor,
!             'Z'=Zonal)
! MESSAGE   Text array containing information to be written to WARNING.OUT
!             file.
! MODEL     Name of CROPGRO executable file
! NH4(L)    Ammonium N in soil layer L (µg[N] / g[soil])
! NL        Maximum number of soil layers = 20
! NO3(L)    Nitrate in soil layer L (µg[N] / g[soil])
! NSTRES    Nitrogen stress factor (1=no stress, 0=max stress)
! NVALP0    Set to 100,000 in PHENOLOG, used for comparison of times of
!             plant stages (d)
! PORMIN    Minimum pore space required for supplying oxygen to roots for
!             optimal growth and function (cm3/cm3)
! RLV(L)    Root length density for soil layer L (cm[root] / cm3[soil])
! RNMODE    Simulation run mode (I=Interactive, A=All treatments,
!             B=Batch mode, E=Sensitivity, D=Debug, N=Seasonal, Q=Sequence)
! RUN       Change in date between two observations for linear
!             interpolation
! RWUEP1    Threshold for reducing leaf expansion compared w/ ratio of
!             TRWU/EP1 (total potential daily root water uptake/ actual
!             transpiration)
! RWUMX     Maximum water uptake per unit root length, constrained by soil
!             water (cm3[water] / cm [root])
! SENESCE   Composite variable containing data about daily senesced plant
!             matter. Structure of variable is defined in ModuleDefs.for
! SNOW      Snow accumulation (mm)
! SOILPROP  Composite variable containing soil properties including bulk
!             density, drained upper limit, lower limit, pH, saturation
!             water content.  Structure defined in ModuleDefs.
! SRAD      Solar radiation (MJ/m2-d)
! ST(L)     Soil temperature in soil layer L (°C)
! STGDOY(I) Day when plant stage I occurred (YYYYDDD)
! SW(L)     Volumetric soil water content in layer L
!            (cm3 [water] / cm3 [soil])
! TMAX      Maximum daily temperature (°C)
! TMIN      Minimum daily temperature (°C)
! TRWUP     Potential daily root water uptake over soil profile (cm/d)
! TWILEN    Daylength from twilight to twilight (h)
! UNH4(L)   Rate of root uptake of NH4, computed in NUPTAK
!            (kg [N] / ha - d)
! UNO3(L)   Rate of root uptake of NO3, computed in NUPTAK (kg [N] / ha -d)
! XHLAI     Healthy leaf area index (m2[leaf] / m2[ground])
! XLAI      Leaf area (one side) per unit of ground area
!            (m2[leaf] / m2[ground])
! YREMRG    Day of emergence (YYYYDDD)
! YREND     Date for end of season (usually harvest date) (YYYYDDD)
! YRPLT     Planting date (YYYYDDD)
!===========================================================================
