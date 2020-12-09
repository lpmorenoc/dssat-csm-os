Module Forecast
  USE ModuleDefs

  Type DailyWeatherType
    INTEGER YRDOY
    REAL RAIN, SRAD, TMAX, TMIN, PAR
  End Type
  TYPE (DailyWeatherType), Allocatable :: Obs_data(:)

  INTEGER FCOUNT, FODAT, FSTART

! Track years of ensemble weather data
  INTEGER EnsYearFirst, EnsYearLast, EnsYearCurrent

CONTAINS

!========================================================================
SUBROUTINE FCAST_STORE(FCOUNT)
! Determine start and end of forecast season (i.e., prior to ensemble)
!   and the first and last years for the historical ensemble weather years
! Store the forecast year weather data plus data from short term forecast file. 
!   These will be used in every year of the forecast ensemble. 

  USE ModuleData

  IMPLICIT NONE

  CHARACTER*1  MEWTH
  CHARACTER*12 FILEW, FILEWC, FILEWG
  CHARACTER*80 PATHWT, PATHWTC, PATHWTG, PATHWTW
  CHARACTER*92 FILEWW
  INTEGER, INTENT(OUT)  :: FCOUNT
  INTEGER DOY, I, Obs_YRDOY, RSEED1, YREND, YR, YRSIM
  INTEGER INCDAT, TIMDIF
  REAL CCO2, DCO2, OZON7, PAR, RAIN, REFHT, RHUM
  REAL SRAD, TAMP, TAV, TDEW, TMAX, TMIN, VAPR, WINDHT
  REAL WINDSP, XELEV, XLAT, XLONG
  LOGICAL SUCCESS
  TYPE (ControlType) CONTROL, CONTROL2

  EXTERNAL IPWTH, YR_DOY

  CALL GET(CONTROL)
  FODAT = CONTROL % FODAT
  YRSIM = CONTROL % YRSIM

  IF (CONTROL % ENDYRS .GT. 1) RETURN

! =======================================================================
! FODAT = date of forecast. Weather prior to this date come from observations. After, they come from ensemble.
! FCOUNT = number of days of "observed weather" to be used in forecast. 
  FSTART = YRSIM
  FCOUNT = TIMDIF(FSTART, FODAT)

! IF FODAT is before or equal to YRSIM, then all weather data are from ensembles.
  IF (FCOUNT .LT. 1) RETURN

  IF (ALLOCATED(Obs_data)) THEN
    DEALLOCATE(Obs_data)
  ENDIF

! Allocate the array size for observed weather data
  Allocate (Obs_data(0:FCOUNT))

! These physical years are needed for historical ensemble, but not for generated weather
! Determine the first and last years of the weather ensemble using simulation start date as the reference.
  CALL YR_DOY(YRSIM,YR,DOY)
! First year of historical ensemble is NYRS before start of simulation date
  EnsYearFirst = YR - CONTROL % NYRS
! Last year of historical ensemble is one year before start of simulation date
  EnsYearLast = YR - 1
! Start with first year of ensemble
  EnsYearCurrent = EnsYearFirst

! =======================================================================
! Get and store weather data between YRSIM and FODAT-1
  CONTROL2 = CONTROL
  CONTROL2 % DYNAMIC = RATE

! IPWTH was already initialized for the forecast year. Go ahead and retrieve 
!   data for this year and store for retrieval at the beginning of each 
!   ensemble.
  DO I = 0, FCOUNT
    Obs_YRDOY = INCDAT(FSTART, I-1)
    CONTROL2 % YRDOY = Obs_YRDOY

    CALL IPWTH(CONTROL2,                              &
        CCO2, DCO2, FILEW, FILEWC, FILEWG, FILEWW,    &    !Output
        MEWTH, OZON7, PAR,                            &    !Output
        PATHWTC, PATHWTG, PATHWTW,                    &    !Output
        RAIN, REFHT, RHUM, RSEED1, SRAD,              &    !Output
        TAMP, TAV, TDEW, TMAX, TMIN, VAPR, WINDHT,    &    !Output
        WINDSP, XELEV, XLAT, XLONG, YREND,            &    !Output
        RATE)                                           

    Obs_data(I) % YRDOY= Obs_YRDOY
    Obs_data(I) % SRAD = SRAD
    Obs_data(I) % TMAX = TMAX
    Obs_data(I) % TMIN = TMIN
    Obs_data(I) % RAIN = RAIN
    Obs_data(I) % PAR  = PAR
  ENDDO

  IF (MEWTH .EQ. 'M') THEN
    CONTROL2 % YRDOY = EnsYearFirst * 1000 + DOY
    CONTROL2 % DYNAMIC = RUNINIT

!   Re-initialize IPWTH for the historical ensemble.
    CALL IPWTH(CONTROL2,                              &
        CCO2, DCO2, FILEW, FILEWC, FILEWG, FILEWW,    &    !Output
        MEWTH, OZON7, PAR,                            &    !Output
        PATHWTC, PATHWTG, PATHWTW,                    &    !Output
        RAIN, REFHT, RHUM, RSEED1, SRAD,              &    !Output
        TAMP, TAV, TDEW, TMAX, TMIN, VAPR, WINDHT,    &    !Output
        WINDSP, XELEV, XLAT, XLONG, YREND,            &    !Output
        RUNINIT)                                           
  ENDIF

  RETURN
END SUBROUTINE FCAST_STORE

!========================================================================
SUBROUTINE FCAST_RETRIEVE(WDATE, RAIN, TMAX, TMIN, SRAD, PAR, YRDOY_F)
  USE ModuleData
  INTEGER DOY, YR, YRDOY, YRDOY_F
  INTEGER WDATE, FCODE, I, TIMDIF
  REAL RAIN, TMAX, TMIN, SRAD, PAR
  TYPE (ControlType) CONTROL
  EXTERNAL :: YR_DOY, TIMDIF

! YRDOY_F = 0 - No observed data for this date, use ensemble data
! YRDOY_F > 0 - This is the date for the historical ensemble (or generated weather)

  YRDOY_F = 0

! No forecast weather data
  IF (FCOUNT .LE. 0) RETURN

! Current date is after forecast date
  IF (WDATE .GE. FODAT) THEN
    CALL GET(CONTROL)
    CALL YR_DOY(WDATE, YR, DOY)
    EnsYearCurrent = EnsYearFirst + CONTROL % ENDYRS - 1
    YRDOY_F = EnsYearCurrent * 1000 + DOY

  ELSE
!   Find today's forecast weather data
    I = TIMDIF(FSTART, WDATE) + 1
    IF (WDATE .NE. Obs_data(I) % YRDOY) THEN
      PRINT *, WDATE, Obs_data(I) % YRDOY, "Error in forecast weather data"
    ENDIF
    
    RAIN = Obs_data(I) % RAIN
    SRAD = Obs_data(I) % SRAD
    TMAX = Obs_data(I) % TMAX
    TMIN = Obs_data(I) % TMIN
    PAR  = Obs_data(I) % PAR
    FCODE = 1
  ENDIF

  RETURN
END SUBROUTINE FCAST_RETRIEVE

!========================================================================
SUBROUTINE FCAST_FINISH()
  DEALLOCATE (Obs_data)
  RETURN
END SUBROUTINE FCAST_FINISH

!========================================================================
End Module Forecast
!========================================================================
