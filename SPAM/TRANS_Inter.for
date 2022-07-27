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
     &    CO2, EO, ET0, EVAP, KTRANS,                     !Input
     &    WINDSP, XHLAI,                                  !Input
     &    WEATHER,                                        !Input
     &    EOP)                                            !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
      USE YCA_Growth_VPD
      IMPLICIT NONE
      
      TYPE (WeatherType) WEATHER
      TYPE (ControlType) CONTROL

      CHARACTER*1  MEEVP
      CHARACTER(len=6), PARAMETER :: ERRKEY = 'IPECO'
      CHARACTER(len=78)  MSG(2)

      INTEGER DYNAMIC
      INTEGER hour, I

      REAL CO2, EO, EVAP, FDINT, KTRANS, TAVG, WINDSP, XHLAI
      REAL EOP, TRAT, EOP_reduc, EOP_max
      REAL KCB, REFET
      REAL PHSV, PHTV, TDEW, TMIN
      REAL, DIMENSION(TS)    ::TAIRHR ,ET0
      REAL, DIMENSION(NumOfCrops)  :: KTRANSM, XHLAIM, TRATM
      CHARACTER*2, CROPS(NumOfCrops)
      REAL KtransbyLAI

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

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
        EOP = 0.0

!***********************************************************************
!***********************************************************************
!     DAILY RATE CALCULATIONS
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
        EOP = 0.0
        TRAT = 0.0 
        KtransbyLAI = 0.0 
        DO I=1, NumOfCrops
            TRATM(I) = TRATIO(CROPS(I), CO2, TAVG, WINDSP, XHLAIM(I))
            TRATM(I) = TRATM(I) * XHLAIM(I) / XHLAI
            TRAT = TRAT + TRATM(I)
            KtransbyLAI = KtransbyLAI + (KTRANSM(I) * XHLAIM(I))
        ENDDO

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

        IF (KCB .GE. 0.0) THEN
          EOP = KCB * REFET !KRT added for ASCE dual Kc ET approach
        ELSE  
          !FDINT = 1.0 - EXP(-(KTRANS) * XHLAI) 
            FDINT = 1.0 - EXP(-KtransbyLAI)
            IF (meevp .NE.'H') THEN 
                EOP = EO * FDINT
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
                  EOPH(hour) = (ET0(hour) * FDINT) * VPDFPHR(hour)
                  EOP = EOP + EOPH(hour)
              ENDDO
          ENDIF
          EOP_reduc = EOP * (1. - TRAT)  
          EOP = EOP * TRAT

C         01/15/03 KJB  I think the change to "Same" K for EOS and EOP
C         may cause next function to be less driving, but below still
C         will depend on whether actual soil evapo (EVAP) meets EOS

!         IF ((EOP + EVAP) .GT. (EO * TRAT)) EOP = EO * TRAT - EVAP

!         Need to limit EOP to no more than EO (reduced by TRAT effect on EOP) 
!         minus actual evaporation from soil, mulch and flood
          EOP_max = EO - EOP_reduc - EVAP
          EOP = MIN(EOP, EOP_max)
        ENDIF  

        EOP = MAX(EOP,0.0)

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

