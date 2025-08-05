!***********************************************************************
!     TRANSPIRATION MODULE Intercrop - File TRANS_Inter.for
!***********************************************************************
!     Includes subroutines:
!         TRANS_Inter - Calculates actual transpiration rate for
!                       intercrops.
!***********************************************************************
C=======================================================================
C  TRANS_Inter, Subroutine modified from TRANS
C  Calculates potential transpiration rate of intercrop, EOP (mm/d).
C-----------------------------------------------------------------------
C  REVISION HISTORY
!  07/19/2022 LPM Modified from TRANS
!-----------------------------------------------------------------------
!  Called by: WATBAL
!  Calls:     None
C=======================================================================
      SUBROUTINE TRANS_Inter(CONTROL, MEEVP,
     &    CO2, EO, ET0, EVAP,                             !Input
     &    WINDSP, XHLAI,                                  !Input
     &    WEATHER,                                        !Input
     &    EOP)                                            !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
      USE YCA_Growth_VPD
      IMPLICIT NONE
      EXTERNAL TRATIO, WARNING
      
      TYPE (WeatherType) WEATHER
      TYPE (ControlType) CONTROL

      CHARACTER*1  MEEVP
      CHARACTER(len=6), PARAMETER :: ERRKEY = 'IPECO'
      CHARACTER(len=78)  MSG(2)

      INTEGER DYNAMIC
      INTEGER hour, I

      REAL CO2, EO, EVAP, FDINT, TAVG, WINDSP, XHLAI
      REAL EOP, TRAT, EOP_reduc, EOP_max
      REAL KCB, REFET
      REAL PHSV, PHTV, TDEW, TMIN
      REAL, DIMENSION(TS)    ::TAIRHR ,ET0
      REAL, DIMENSION(NumOfCrops)  :: KTRANSM, XHLAIM, TRATM
      REAL, DIMENSION(NumOfCrops)  :: EOPM, FDINTM, EOP_reducM
      REAL, DIMENSION(NumOfCrops)  :: FracIntRadMTrans
      CHARACTER*2, CROPS(NumOfCrops)

!     FUNCTION SUBROUTINES:
      REAL TRATIO
       
      DYNAMIC = CONTROL % DYNAMIC
      CROPS = CONTROL % INTERCROP
      TAVG   = WEATHER % TAVG  
      TDEW   = WEATHER % TDEW   
      TMIN   = WEATHER % TMIN 
      TAIRHR = WEATHER % TAIRHR
      
      CALL GET('SPAM', 'KCB', KCB)
      CALL GET('SPAM', 'REFET', REFET)
      Call GET('PLANT', 'KTRANSM',  KTRANSM, NumOfCrops)
      Call GET('PLANT', 'XHLAIM',  XHLAIM, NumOfCrops)
      Call GET('PLANT', 'FracIntRadMTrans',FracIntRadMTrans, NumOfCrops)

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
        EOP = 0.0
        EOPM = 0.0

!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
        EOP = 0.0
        EOPM = 0.0
        TRATM = 0.0 
        EOP_reduc = 0.0
       DO I=1, NumOfCrops
            TRATM(I) = TRATIO(CROPS(I), CO2, TAVG, WINDSP, XHLAIM(I)) 

!-----------------------------------------------------------------------
C       Estimate light interception.  NOTE 01/15/03 We don't want PAR
C       Estimate ENERGY interception.  NOTE 01/15/03 We don't want PAR
C       interception.  Changed to reflect energy interception which
C       correctly considers both PAR and the infrared.  Villalobos
C       measured 0.52 for Transpiration component.  Theory of Goudriaan
C       supports combined interception coefficient of 0.5 for PAR and NIR

C       01/15/03 - Work of Sau et al, shows that a K of 0.5 was better in
C       all cases, for PT form as well as the Dynamic form for predicting
C       soil water balance and predicting measured ET.

!       LPM 07/27/2022 Currently KCB is a species coefficient and we would
!       need to modify the ASCE method so it can work for intercropping

        IF (KCB .GE. 0.0) THEN 
          EOPM(I) = KCB * REFET !KRT added for ASCE dual Kc ET approach
        ELSE  
          !FDINT = 1.0 - EXP(-(KTRANS) * XHLAI) 
            FDINTM(I) = FracIntRadMTrans(I)
            IF (meevp .NE.'H') THEN 
                EOPM(I) = EO * FDINTM(I)
            ELSE
              CALL GET('SPAM', 'PHSV' ,phsv)
              CALL GET('SPAM', 'PHTV' ,phtv)
      
              IF (phsv <= 0.0) THEN
                  MSG(1) = "VPD sensitivity parameter PHSV" //
     &              " is not defined for EVAPO method (H)."
                  MSG(2) = "Program will stop."
                  CALL WARNING(2, ERRKEY, MSG)
                  CALL ERROR(ERRKEY,4,"",0)
              ENDIF
              IF (phtv <= 0.0) THEN
                  MSG(1) = "VPD threshold parameter PHTV is" //
     &              " not defined for EVAPO method (H)."
                  MSG(2) = "Program will stop."
                  CALL WARNING(2, ERRKEY, MSG)
                  CALL ERROR(ERRKEY,4,"",0)
              ENDIF
              DO hour = 1,TS 
                  VPDFPHR(hour) =  get_Growth_VPDFPHR(PHSV, PHTV, TDEW, 
     &                     TMIN, TAIRHR, hour)
                  EOPH(hour) = (ET0(hour) * FDINTM(I)) * VPDFPHR(hour)
                  EOPM = EOPM + EOPH(hour)
              ENDDO
          ENDIF
          EOP_reducM(I) = EOPM(I) * (1. - TRATM(I))  
          EOPM(I) = EOPM(I) * TRATM(I)
          EOP_reduc = EOP_reduc + EOP_reducM(I) 
          ENDIF
          EOP = EOP + EOPM(I)
       ENDDO

C         01/15/03 KJB  I think the change to "Same" K for EOS and EOP
C         may cause next function to be less driving, but below still
C         will depend on whether actual soil evapo (EVAP) meets EOS

!         IF ((EOP + EVAP) .GT. (EO * TRAT)) EOP = EO * TRAT - EVAP

!         Need to limit EOP to no more than EO (reduced by TRAT effect on EOP) 
!         minus actual evaporation from soil, mulch and flood
        IF (KCB .GE. 0.0) THEN 
          EOP = EOP
        ELSE
          EOP_max = EO - EOP_reduc - EVAP
          EOP = MIN(EOP, EOP_max)
        ENDIF 
        
        EOP = MAX(EOP,0.0)
        
        IF (EOP < SUM(EOPM) .AND. SUM(EOPM) > 0.0) THEN
            DO I=1, NumOfCrops
                EOPM(I) = (EOPM(I)/SUM(EOPM)) * EOP
            ENDDO
        ENDIF
        Call PUT('SPAM', 'EOPM',  EOPM, NumOfCrops)
        

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE TRANS_Inter

!-----------------------------------------------------------------------
!     TRANS VARIABLE DEFINITIONS:
!-----------------------------------------------------------------------
! CO2     Atmospheric carbon dioxide concentration (ppm)
! EO      Potential evapotranspiration rate (mm/d)
! EOP     Potential plant transpiration  (mm/d)
! EVAP    Actual soil + mulch + flood evaporation rate (mm/d)
! FDINT   Fraction of direct solar radiation captured by canopy
! KCAN    Canopy light extinction coefficient for daily PAR, for
!           equidistant plant spacing, modified when in-row and between row
!           spacings are not equal
! KCB     Basal crop coefficient for ASCE dual Kc ET method
! LNUM    Current line number of input file
! REFET   ASCE Standardized Reference Evapotranspiration (alfalfa or grass)
! TAVG    Average daily temperature (�C)
! TRAT    Relative transpiration rate for CO2 values other than 330 ppm
! TRATIO  Function subroutine which calculates relative transpiration rate.
!
! WINDSP  Wind speed (km/d)
! XHLAI   Leaf area index (m2[leaf] / m2[ground])
!-----------------------------------------------------------------------
!     END SUBROUTINE TRANS_Inter
!-----------------------------------------------------------------------

