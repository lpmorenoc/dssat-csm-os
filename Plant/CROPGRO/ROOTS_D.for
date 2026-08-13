C=======================================================================
C  ROOTSD, Subroutine, K.J. BOOTE - DYNAMIC 1-D ROOTING
C-----------------------------------------------------------------------
C
C  Calculates root growth, extension, respiration, and senescence
C
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  01/09/1989 GH  Written.
C  09/29/1995 KJB Changed to lessen effect of water deficit on root depth
C                 increase.
C  01/19/1996 JWJ Added effects of excess water.
C  01/20/1996 KJB Increase root extension rate under drought stress
C  09/13/1998 CHP Modified for modular format
C  09/14/1998 CHP Changed TROOT to TRLV to match same variable in ROOTDM
C                 Changed SWDF1 to SWFAC to match variable name in WATBAL
C  05/11/1999 GH  Incorporated in CROPGRO
!  02/21/2005 SJR Moved ISWWAT condition here to allow computation of 
!                 root senescence even when water not simulated.
!  10/04/2005 SJR Include senescence due to water stress in total 
!                 daily senescence.
!  10/20/2005 CHP Added optional minimum root mass for senescence, 
!                 RTWTMIN, to species file
!  01/19/2006 CHP Fixed discrepancies between plant root senescence  
!                 calculated and that sent to soil routines for addition
!                 to organic matter.
!  12/14/2017 KJB Replaced with Dynamic 1-D rooting from Jones & Bland 
!  10/24/2024 CHP Added TRLV to PlantGro.OUT
!-----------------------------------------------------------------------
!  Called by  :  PLANT
!  Calls      :  IPROOTD, INROOTD
!=======================================================================

      SUBROUTINE ROOTS_D(DYNAMIC,
     &    AGRRT, CROP, DLAYR, DS, DTX, DUL, FILECC, FRRT, !Input
     &    BD,CLAY,SILT,STONES,TOTBAS, EXTAL, EXCA,        !Input
     &    ISWWAT, LL, NLAYR, PG, PLTPOP, RO, RP, RTWT,    !Input
     &    SAT, SW, SWFAC, VSTAGE, WR, WRDOTN, WTNEW,DAS,ST,      !Input
     &    RLV, RTDEP, SATFAC, SENRT, SRDOT, TRLV)               !Output
!kjb
C-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      IMPLICIT NONE
      SAVE
                
      CHARACTER*1 ISWWAT
      CHARACTER*2 CROP
      CHARACTER*92 FILECC
      CHARACTER*3   TYPRDT

      INTEGER L, L1, NLAYR

      INTEGER DYNAMIC

      REAL CUMDEP, DEPMAX, DTX, FRRT,
     &  PG, RFAC1, RFAC2, RFAC3,
     &  RLDSM, RLNEW, RO, RP,
     &  RTDEP, RTSDF, RTSEN, RTWT, SRDOT, SWDF, SWFAC,
     &  TRLDF, TRLV, WRDOTN   !, TRTDY
      REAL CGRRT, AGRRT
      REAL SWEXF, PORMIN, RTEXF, RTSURV
      REAL RTDEPI, SUMEX, SUMRL, SATFAC
      REAL PLTPOP, WTNEW, VSTAGE
      REAL TABEX              !Function subroutine located in UTILS.for
      REAL XRTFAC(4), YRTFAC(4), FNRDT(4)
      REAL DLAYR(NL), DS(NL), DUL(NL), ESW(NL), LL(NL), RLDF(NL)
      REAL RLGRW(NL), RLSEN(NL), RLV(NL), RLV_WS(NL), RRLF(NL)
      REAL SW(NL), SAT(NL), WR(NL)
      REAL GRESPR(NL), MRESPR(NL), RESPS(NL)
      REAL SENRT(NL)

!     Added 10/20/2005 for minimum RLV calculations
!     RTWTMIN = minimum root mass per layer; used to limit senescence
!                 (g/m2) (species file parameter)
!     TRLV_MIN  = conversion of RTWTMIN to RLV units per layer
      REAL TRLV_MIN, RLSENTOT, FACTOR, RTWTMIN
      REAL TotRootMass, CumRootMass
!kjb-
      TYPE (SoilType)    SOILPROP
C
      CHARACTER*56 TITLE
      CHARACTER*20 OUTFILE
      CHARACTER*1 SELECT
!      REAL Vstage, WRDOTN,LL
      REAL, DIMENSION(NL) :: NO3, NH4, ST
!      REAL, DIMENSION(10) :: NO3, NH4, ST, SW
      INTEGER DAY, DAS, YEAR, Rstage
      LOGICAL GOOD
C     
      INTEGER I, IJ, IR
!kjb
      REAL WCG,PD,GS,GSR,GSY
      REAL CAA,CAX,ALA,ALX,SFT
      REAL,DIMENSION(NL) :: BD,SAN,EXCA,ALS,BDO,
     ! BDX,LT,SST,LW,STP,SCD,PO,CWP,SCA,SAL,SBD,SCF,ZA,LWN,LWF,
     ! SAI,ASF1,ASF2,ASF3,LWR,CEC,CAS,WFP,ASF,CLAY,SILT, STONES,
     ! TOTBAS, EXTAL

!      IMPLICIT NONE
C      REAL DM
C
      REAL,DIMENSION(NL) :: WFL,DZ,GPL,LWA,DWL,
     1  RWL,DLL,RLL,GAL

      REAL LWS,LWM,DDI,DRD,WFT,TEMP,GPS,TRW
C
      REAL DM, DAC, DMA, DMC, DMD, RD
C
C      COMMON /BLK1/ CAA,CAX,ALA,ALX,TBS,TOP,SFT,BD,SAN,SIL,ROK,
C     1  SMB,EAL,CA,CLA,ALS,BDO,BDX
C      COMMON /BLK2/ UL,LT,SST,LW,STP,SCD
C      COMMON /BLK12/ PO,CWP,SCA,SAL,SBD,SCF
C      COMMON /BLK13/ Z,ZA
C      COMMON /BLK23/ LWF,SAI,ASF1,ASF2,ASF3
C      COMMON /BLK123/ IJ

C
C  INITIALIZE VARIABLES
C
C     DS(1) = 0.       
      IJ = 0
      IR = 0
      DAY = 0.
      DMA = 0.
      DAC = 0.
      GSY = 0.
C      DMC = 0.
!kjb
!***********************************************************************
!***********************************************************************
!     Run Initialization - Called once per simulation
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
      CALL IPROOTD(FILECC,                                    !Input
     &  GSR, LWM, LWS, WCG, SFT, CAA, CAX, ALA, ALX,          !Output
     & PORMIN, RFAC1, RLDSM, RTDEPI, RTEXF,                   !Output
     &  RTSEN, RTSDF, RTWTMIN, XRTFAC, YRTFAC)                !Output

      DEPMAX = DS(NLAYR)

!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      SRDOT = 0.0       
      RLV   = 0.0
      RTDEP = 0.0       
      SENRT = 0.0
      SUMEX = 0.0
      SUMRL = 0.0
      SATFAC = 0.0
      RTDEP = 0.

!kjb
!    DO 20 I=1,10
!       RWL(I) = 0.
!       RLL(I) = 0.
!       RLV(I) = 0.
!   20 CONTINUE
C
!     Open input files
!      CALL DailyInput(1, 
!     &    YEAR, DOY,              !Date
!     &    WRDOTN, VStage, RStage, !New root dry matter (g/m2/d)
!     &    ST,                     !Soil temperature (oC)
!     &    SW,                     !Soil water (mm3/mm3)
!     &    NO3, NH4,               !Soil N (ppm)
!     &    GOOD)                   !Signal success
!kjb - prior 7-8 lines, should now be coming from CROPGRO, not read in
!kjb - do not know about the setting to zero
!
!kjb - next need to read in the static soil parameters, at least once per season
!      Need to consider what each of these does.
!
C KJB MOVED HERE
      OutFile = "Roots.OUT"
C      OPEN(3,FILE=OUTFILE,ACCESS='SEQUENTIAL',STATUS='NEW')  !DOES NOT HELP
      OPEN(3,FILE=OUTFILE,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      
 1000 FORMAT (A56)
C
C  READ AND WRITE THE GENETIC PARAMETERS
C

      WRITE (3,7000)
      WRITE (3,5000) GSR,LWM,LWS,WCG,CAA,CAX,ALA,ALX, SFT
 5000 FORMAT (2X,'GSR = ',F5.1,' LWM = ',F5.0,
     1 ' LWS = ',F5.0,' WCG = ',F5.1,
     2 ' CAA = ',F5.1,' CAX = ',F5.1,
     3 'ALA = ',F5.0,' ALX = ',F5.0,/,
     4 ' SFT = ',F5.1,//)
 7000 FORMAT (//,28X,'GENETIC CHARACTERISTICS',//)
C
C      RD = PD
!
!      DO 30 I=1,10
      DO 30 I=1,NL
       LWR(I) = LWS
       IF (WR(I).GT.0) THEN
           SCD(I) = 1.0
       ELSE
           SCD(I) = 0.0
       ENDIF
   30 CONTINUE
  110 CONTINUE

C  CALCULATE STATIC SOIL STRESS FACTORS
C
!kjb - used passing arguments, to and from, rather than common here
C    
      CALL SCALC(NLAYR,DS,CAA,CAX,ALA,ALX,BD,CLAY,SILT,STONES,   !Input
     1  TOTBAS,EXTAL,EXCA,                                       !Input
     1  ZA,CWP,PO,CEC,ALS,SAL,CAS,SCA,BDO,BDX,SBD,SCF)        !Output
C
C  HEADING AND OUTPUT FROM SUBROUTINE SCALC
C
C kjb, need to replace 11 with NLAYR in all places, but need to read in?
C kjb, maybe output NLAYR  ?
      WRITE (3,7010)
      WRITE (3,9000)
      DO 200 I = 1,11
      WRITE (3,3010) I,SAN(I),SILT(I),CLAY(I),STONES(I),BD(I),
     1   PO(I),DUL(I),LL(I),CWP(I),CEC(I)
  200 CONTINUE
      WRITE (3,9010)
      DO 210 I=1,11
      WRITE (3,3020) I,TOTBAS(I),EXTAL(I),ALS(I),EXCA(I),SAL(I),
     1   SCA(I),SBD(I),SCF(I),SCD(I),BDO(I),BDX(I)
  210 CONTINUE
C
! 1000 FORMAT (A56) !twice
 3010 FORMAT (2X,I3,4X,1X,3(3X,F4.1),2X,F4.2,3X,F4.2,
     1 3X,F4.2,3(3X,F4.2),3X,F4.1)
 3020 FORMAT (2X,I3,2X,4(3X,F4.1),7(3X,F4.2))
! 5000 FORMAT (2X,'RDX = ',F5.1,' GSR = ',F5.1,' LWM = ',F5.0, !twice
!     1 'LWS = ',F5.0,' WCG = ',F5.1,' TBS = ',F5.0,/,
!     2 'TOP = ',F5.0,' CAA = ',F5.1,' CAX = ',F5.1,
!     3 'ALA = ',F5.0,' ALX = ',F5.0,' PD = ',F5.2,/,
!     4 'GSD = ',F5.1,' SFT = ',F5.1,//)
! 7000 FORMAT (//,28X,'GENETIC CHARACTERISTICS',//)
 7010 FORMAT (//,29X,'INPUT CHARACTERISICS')
 9000 FORMAT (//,2X,'LAYER',5X,'SAN',4X,'SIL',4X,'CLA',4X,
     1 'ROK',4X,'BD',5X,'PO',5X,'DUL',5X,'LL',4X,'CWP',4X,'CEC'/)
 9010 FORMAT (//,2X,'LAYER',4X,'SMB',4X,'EAL',4X,'ALS',5X,'CA',4X,
     1 'SAL',4X,'SCA',4X,'SBD',4X,'SCF',4X,'SCD',4X,'BDO',4X,'BDX',/)
! 9020 FORMAT (1X,'DPTH',3X,'RLV',5X,'GAL',6X,'DWL',6X,'RWL',7X,
!     1 'LWR',5X,'STP',5X,'SST',5X,'SAI',5X,'ASF',/)
C      WRITE (3,7070)
C      WRITE (3,7080)
C 7070 FORMAT (//,29X,'STRESS OUTPUT CHARACTERISICS')
C 7080 FORMAT (//,'DOY',1X,'LAYER',4X,'LL',3X,'DUL',3X,'SW',3X,'ST',3X,
C     1 'LWF',3X,'WFP',3X,'SST',3X,'SAI',3X,'STP',3X,'ASF',3X,'ASF1',
C     1  3X,'ASF2',3X,'ASF3',/)
C
      WRITE (3,8010)
C      WRITE (3,8020)
      WRITE (3,8015)
 8010 FORMAT (//,29X,'ROOT DISTRIBUTION CHARACTERISICS')
 8015 FORMAT(//'    GS,   DDI,   DRD,  RTDEP, WFT, TEMP, DMD,  GPS, 
     1     TRW, DMA,   DMC')
C
 8020 FORMAT (//,'DOY',1X,'LAYER',4X,'LL',3X,'DUL',4X,'SW',3X,'ST',3X,
     1 'LWF',3X,'WFL',4X,'LWA',3X,'LWR',4X,'DWL',4X,'RWL',4X,'DLL',
     1  4X,'RLL',4X,'RLV',4X,'GPL',4X,'GAL',/)
       
C

!kjb
! Next is just species, and will be modified by the static and dynamic soil 
!-----------------------------------------------------------------------
C     ROOT DEPTH INCREASE RATE WITH TIME, cm/physiological day
C-----------------------------------------------------------------------
!      IF (CROP .NE. 'FA' .AND. ISWWAT .NE. 'N') THEN
      IF (CROP .NE. 'FA') THEN
        RFAC2 = TABEX(YRTFAC,XRTFAC,0.0,4)
      ENDIF
      
      CumRootMass = 0.0
C

!***********************************************************************
!***********************************************************************
!     EMERGENCE CALCULATIONS - Performed once per season upon emergence
!         or transplanting of plants
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. EMERG) THEN
!-----------------------------------------------------------------------
!       Call INROOT for initialization of root variables on
!       day of emergence.  (GROW emergence initialization
!       must preceed call to INROOT.)
!-----------------------------------------------------------------------
      CALL INROOTD(
     &  DLAYR, FRRT, NLAYR, PLTPOP, RFAC1, RTDEPI, WTNEW, !Input
     &  RLV, RTDEP)                                       !Output

      RFAC3 = RFAC1

      TRLV = 0.0
      DO L = 1,NLAYR
        TRLV = TRLV + RLV(L) * DLAYR(L) ! cm[root] / cm2[ground]
      ENDDO
C    
      CumRootMass = WTNEW * FRRT * PLTPOP * 10. 
!        kg[root]  g[tissue] g[root]    plants   kg/ha
!        -------- = ----- * --------- * ------ * ----- 
!           ha      plant   g[tissue]     m2      g/m2
C
      CALL DCALC(NLAYR,LL,DUL,SW,ST,TYPRDT,FNRDT,SFT,
     1  CWP,PO,SCA,SAL,SBD,SCF,SCD,              !inputs
     1  LWF,WFP,SST,SAI,STP,ASF,ASF1,ASF2,ASF3)  !outputs

      RTDEP = RTDEPI
      GS = VSTAGE   !Divide by max vstage for 0-1 fraction
      DMD = WTNEW * FRRT * PLTPOP         !WRDOTN is g/m2, and math below for GAL(i) use
						!LWA(i) which is in m/g, without adjusting
						!for wrong units.
C 
C
      CALL RTDIST(NLAYR,DS,ZA,GSR,GSY,GS, WRDOTN,   !inputs
     1  WCG,LWS,LWM,DMD,ASF,ASF1,ASF2,ASF3,SAI,LWF, RFAC2, DTX,     !inputs
     1  DDI,DRD,RTDEP,WFL,DZ,WFT,LWN,LWA,RLV,TEMP,GPL,GPS,  !outputs
     1  GAL,RWL,DMA,DWL,DLL,RLL,LWR,TRW,DMC,DAC)         !outputs
C
C
C  HEADING FOR OUTPUT FROM SUBROUTINE RTDIST
C      WRITE (3,7020)
C 7020 FORMAT (//,25X,'DAILY ROOT DISTRIBUTION OUTPUT',//)
C      WRITE(3,9030)
C 9030 FORMAT(" YEAR DOY    GS   RD    TRW    DMA    DMD    DAC",
C     &  "    DMC  DPTH    RLV     GAL      DWL      RWL       LWR",
C     &  "     STP     SST     SAI     ASF")
C 
!***********************************************************************
!***********************************************************************
!     DAILY RATE/INTEGRATION
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
C kjb -definitions
C      LW = SW             !soil water -- same units
C      UL = DUL            !soil water -- same units
C      LT = ST             !soil temperature -- same units
C      PD = RTDEPI * 0.01*.75
      GS = VSTAGE   !Divide by max vstage for 0-1 fraction
      DMD = WRDOTN        !WRDOTN is g/m2, and math below for GAL(i) use
						!LWA(i) which is in m/g, without adjusting
						!for wrong units.
C   
C   
C      
      CALL DCALC(NLAYR,LL,DUL,SW,ST,TYPRDT,FNRDT,SFT,
     1  CWP,PO,SCA,SAL,SBD,SCF,SCD,              !inputs
     1  LWF,WFP,SST,SAI,STP,ASF,ASF1,ASF2,ASF3)  !outputs

C      DO 777 I = 1,11      
C      WRITE (3,7090) DAS,I,LL(I),DUL(I),SW(I),ST(I),LWF(I),
C     1   WFP(I),SST(I),SAI(I),STP(I),ASF(I),ASF1(I),ASF2(I),
C     1   ASF3(I)
C  777 CONTINUE
C 7090 FORMAT (1X,I3,2X,I3,1X,3(2X,F4.2),2X,F4.1,9(2X,F4.2))
C 
C    WRDOTN SAME AS DMD
C    PD SAME AS RTDEPI (INPUT), AS MODEL INITIATES AT EMERGENCE
C      PD = RTDEP * 0.01
C      RD = PD
C
C     Calculate Root Depth Rate of Increase, Physiological Day (RFAC2)
C-----------------------------------------------------------------------
      RFAC2 = TABEX(YRTFAC, XRTFAC, VSTAGE, 4)
      
      CALL RTDIST(NLAYR,DS,ZA,GSR,GSY,GS, WRDOTN,    !inputs
     1  WCG,LWS,LWM,DMD,ASF,ASF1,ASF2,ASF3,SAI,LWF, RFAC2, DTX,     !inputs
     1  DDI,DRD,RTDEP,WFL,DZ,WFT,LWN,LWA,RLV,TEMP,GPL,GPS,  !outputs
     1  GAL,RWL,DMA,DWL,DLL,RLL,LWR,TRW,DMC,DAC)         !outputs
C
       WRITE (3,8111) GS,DDI,DRD,RTDEP,WFT,TEMP,DMD,GPS,TRW,DMA,DMC
 8111 FORMAT(3(1X,F6.2),1X,F6.1,3(1X,F4.2),2X,F5.2,1X,F6.1,1X,
     1     F6.2,1X,F6.1/)
C      DO 888 I = 1,NLAYR      
C      WRITE (3,8880) DAS,I,LL(I),DUL(I),SW(I),ST(I),LWF(I),
C     1   WFL(I),LWA(I),LWR(I),DWL(I),RWL(I),DLL(I),RLL(I),
C     1   RLV(I),GPL(I),GAL(I)
C  888 CONTINUE
 8880 FORMAT (1X,I3,2X,I3,3(2X,F4.2),2X,F4.1,2(2X,F4.2),
     1 2(1X,F5.1),3(1X,F6.2),1X,F7.2,3(1X,F6.2))
      
C      REAL,DIMENSION(NL) :: WFL,DZ,GPL,LWA,LWR,DWL,
C     1  RWL,DLL,RLL,RLV,GAL,Z,ZA,
C     1  LWF,SAI,ASF,ASF1,ASF2,ASF3
C      REAL LWS,LWM,LWN,DDI,DRD,WFT,TEMP,GPS,TRW
C      REAL WCG,PD,GSD,GS,GSR,GSY
C      REAL DM, DAC, DMA, DMC, DMD, RD, RDX
C      INTEGER I, IJ, IRC kjb, old next
C

      RLNEW = WRDOTN * RFAC1 / 10000.
      CGRRT = AGRRT * WRDOTN

C-----------------------------------------------------------------------
C     Calculate root length per cm2 soil and initiate growth,
C     respiration and senescence by layer
C-----------------------------------------------------------------------
!      TRTDY = 0.0
!     1/19/2006 Remove TRTDY and replace with TRLV -- RLV is only updated
!     once, so yesterday's value is stored in TRLV here.
      DO L = 1,NLAYR
!       TRTDY = TRTDY + RLV(L) * DLAYR(L) ! cm[root] / cm2[ground]
        RRLF(L)   = 0.0
        RLSEN(L)  = 0.0
        RLGRW(L)  = 0.0
        MRESPR(L) = 0.0
        GRESPR(L) = 0.0
        RESPS(L)  = 0.0
      ENDDO

!     Update RFAC3 based on yesterday's RTWT and TRLV
      IF (RTWT - WRDOTN .GE. 0.0001 .AND. TRLV .GE. 0.00001) THEN
!       RFAC3 = TRTDY * 10000.0 / (RTWT - WRDOTN)
!       RTWT has not yet been updated today, so use yesterday's
!       value and don't subtract out today's growth - chp 11/13/00
        RFAC3 = TRLV * 10000.0 / RTWT
      ELSE
        RFAC3 = RFAC1
      ENDIF

!     10/20/2005 Limit RLV decrease due to senscence to 
!       a minimum resulting root weight 
      IF (RTWTMIN > 0.0) THEN
        TRLV_MIN = RTWTMIN * RFAC3 / 1.E4   !same units as TRLV
!        cm/cm2  =  (g/m2) *(cm/g) / (cm2/m2)
      ELSE
!       Set TRLV_MIN to zero -- no minimum root mass
        TRLV_MIN = 0.0
      ENDIF
      
!-----------------------------------------------------------------------
      TRLDF  = 0.0
      CUMDEP = 0.0
      SUMEX  = 0.0
      SUMRL  = 0.0
      RLV_WS = 0.0
      RLSEN  = 0.0

      DO L = 1,NLAYR
        L1 = L
        CUMDEP = CUMDEP + DLAYR(L)
        SWDF = 1.0
        SWEXF = 1.0

C-----------------------------------------------------------------------
C     2/21/05 - SJR - move conditional call for water stress from CROPGRO 
C     to ROOTS.  Allows root senescence when Water dynamics option is 
C     turned off.  Water stress options set to no stress levels.  This 
C     also allows output of root growth dynamics without limimiting 
C     water or N uptake. 
C-----------------------------------------------------------------------
        IF (ISWWAT .EQ. 'Y') THEN
          IF (SAT(L)-SW(L) .LT. PORMIN) THEN
            SWEXF = (SAT(L) - SW(L)) / PORMIN
            SWEXF = MIN(SWEXF, 1.0)
          ENDIF

          SUMEX = SUMEX + DLAYR(L) * RLV(L) * (1.0 - SWEXF)
          SUMRL = SUMRL + DLAYR(L) * RLV(L)

          ESW(L) = DUL(L) - LL(L)
          IF (SW(L) - LL(L) .LT. 0.25*ESW(L)) THEN
            SWDF = (SW(L) - LL(L)) / (0.25*ESW(L))
            SWDF = MAX(SWDF, 0.0)
          ENDIF
        ENDIF
C-----------------------------------------------------------------------

        RTSURV = MIN(1.0,(1.-RTSDF*(1.-SWDF)),(1.-RTEXF*(1.-SWEXF)))
        IF (RLV(L) .GT. RLDSM .AND. TRLV + RLNEW > TRLV_MIN) THEN
!         1/14/2005 CHP Don't subtract water stress senescence 
!           yet - combine with natural senescence and check to see if 
!           enough RLV for senescence to occur (TRLV > TRLV_MIN)
          !RLV(L) = RLV(L) * RTSURV
          RLV_WS(L) = RLV(L) * (1.0 - RTSURV)
        ELSE
          RLV_WS(L) = 0.0
        ENDIF

C-----------------------------------------------------------------------
        RLDF(L) = WR(L) * DLAYR(L) * MIN(SWDF,SWEXF)
        IF (CUMDEP .LT. RTDEP) THEN
          TRLDF = TRLDF + RLDF(L)
        ELSE
          IF (WR(L) .GT. 0.0 .AND. RLNEW .GT. 0.0) THEN
C            IF (L .EQ. 1) THEN
C              RTDEP = RTDEP + DTX * RFAC2
C kjb Not using RTDEP from ROOTS main
C      RTDEP = RD * 100.
C            ELSE
C              RTDEP = RTDEP + DTX * RFAC2 * MIN(SWDF,SWEXF) *
C     &                (1. + 0.25 * (1. - MAX(SWFAC,0.40)))
C kjb Not using RTDEP from ROOTS main
C      RTDEP = RD * 100.
C-----------------------------------------------------------------------
C-KJB  DO NOT WANT TO DECREASE ROOT DEPTH WITH STRESS.  IF PG TO ROOTS
C IS LOW BECAUSE OF SEED GROWTH OR IF WATER DEFICIT CAUSES LOW PG TO ROOTS
C DESPITE INCREASED PARTITIONING TO ROOTS, THEN RLV WILL NOT INCREASE
C SO THERE WILL BE NO EFFECTIVE INCREASE IN WATER EXTRACTION.  IDEALLY THE
C DECISION SHOULD BE BASED ON AMOUNT OF ROOT GROWTH VS NORMAL UNSTRESS.
C ACCELERATE FROM 1.0 TO 0.5, STAY FLAT, SHOULD DROP AGAIN, 0.5 TO 0.0
C AS THE OTHER FUNCTION ACTS.  NOTE:  SWFAC*2.0 WAS NOT USED IN ALL 40 CASES
C EXCEPT 1985-RAINFED WHERE DELETING INCREASED YIELD 2764 TO 2770 KG/HA.
C NOW ACCELERATING ROOT GROWTH BY ABOUT 12-13% AT SWFAC=0.50.  THIS
C HELPS IOWA 88 AND VEG STRESS TRTS IN 1981 AND 1985. INCR SEED AND BIO.
C-----------------------------------------------------------------------
C            ENDIF
            IF (RTDEP .GT. DEPMAX) THEN
               RTDEP = DEPMAX
C kjb Not using RTDEP from ROOTS main
C      RTDEP = RD * 100.
            ENDIF
          ENDIF
          RLDF(L) = RLDF(L) * (1. - (CUMDEP - RTDEP) / DLAYR(L))
          TRLDF = TRLDF + RLDF(L)
          GO TO 2900
        ENDIF
       ENDDO
C-----------------------------------------------------------------------
C     Calculate root senescence, growth, maintenance and growth
C     respiration, and update root length density for each layer.
!-----------------------------------------------------------------------
 2900 CONTINUE

      IF (SUMRL .GT. 0.0) THEN
         SATFAC = SUMEX/SUMRL
      ELSE
         SATFAC = 0.0
      ENDIF

      SRDOT = 0.0
      RLSENTOT = 0.0

      DO L = 1,L1
        IF (TRLDF .LT. 0.00001) THEN
          RRLF(L) = 1.0
        ELSE
          RRLF(L) = RLDF(L)/TRLDF
        ENDIF
!-----------------------------------------------------------------------
!       MRESPR, GRESPR, and RESPS are not used anywhere
!                       chp 9/22/98
!-----------------------------------------------------------------------
        MRESPR(L) = (RLV(L)/RFAC1*RO*DLAYR(L)*100.0
     &    +RRLF(L)*FRRT*PG*RP) * 44.0 / 30.0
        GRESPR(L) = RRLF(L) * (CGRRT-WRDOTN) * 44.0 /30.0
        RESPS(L) = MRESPR(L) + GRESPR(L)
!-----------------------------------------------------------------------
        RLGRW(L) = RRLF(L) * RLNEW / DLAYR(L) !cm[root]/cm3[ground]

        IF (TRLV + RLNEW > TRLV_MIN) THEN
          RLSEN(L) = RLV(L) * RTSEN * DTX
        ELSE
          RLSEN(L) = 0.0
        ENDIF

!       Limit total senescence in each layer to existing RLV
        IF (RLV(L) - RLSEN(L) - RLV_WS(L) + RLGRW(L) < 0.0) THEN
          RLSEN(L) = RLV(L) + RLGRW(L) - RLV_WS(L)
        ENDIF 

!       RLSENTOT is profile senescence, water stress and natural cm/cm2
        RLSENTOT = RLSENTOT + (RLSEN(L) + RLV_WS(L)) * DLAYR(L)
      ENDDO

!     If senescence too high (results in TRLV < TRLV_MIN) then
!       reduce senescence in each layer by factor.
      IF (RLSENTOT > 1.E-6 .AND. TRLV + RLNEW - RLSENTOT < TRLV_MIN)THEN
        FACTOR = (TRLV + RLNEW - TRLV_MIN) / RLSENTOT
        FACTOR = MAX(0.0, MIN(1.0, FACTOR))
        RLSEN  = RLSEN  * FACTOR
        RLV_WS = RLV_WS * FACTOR
      ENDIF

!     Update RLV and TRLV based on today's growth and senescence
      TRLV = 0.0
      DO L = 1, NLAYR
C kjb - This turned off root growth and senescence from this subroutine
C        RLV(L) = RLV(L) + RLGRW(L) - RLSEN(L) - RLV_WS(L)
        TRLV = TRLV + RLV(L) * DLAYR(L)

!       Keep senescence in each layer for adding C and N to soil
        !SENRT(L) = RLSEN(L) * DLAYR(L) / RFAC1 * 10000. * 10. !kg/ha
!       1/14/2005 CHP - water stress senesence needs to be inlcuded.
        SENRT(L) = (RLSEN(L) + RLV_WS(L)) * DLAYR(L) / RFAC3 * 1.E5 
!                   cm[root]              g[root]   1000 cm2   10(kg/ha)
!         kg/ha  =  -------- * cm[soil] * ------- * -------- * ---------
!                  cm3[soil]             cm[root]     m2         (g/m2)

        SENRT(L) = AMAX1(SENRT(L), 0.0)
        SRDOT = SRDOT + SENRT(L)/10.        !g/m2

!       Not used:
        !TRLGRW = TRLGRW + RLGRW(L) * DLAYR(L)
        !TRLSEN = TRLSEN + RLSEN(L) * DLAYR(L)
      ENDDO

!     11/13/2000 CHP Sum RLSEN for total root senescence today.  
!     SRDOT = TRLSEN / RFAC3 * 10000.     !g/m2

!     Total root senescence = water stress + natural senescence
!     10/3/2005 SJR
!      SRDOT = (TRTDY + RLNEW - TRLV) * 10000.0 / RFAC3    !g/m2
      SRDOT = AMAX1(SRDOT, 0.0)

      TotRootMass = TRLV / RFAC3 * 1.E5
!                   cm[root]   g[root]   10000 cm2   10(kg/ha)
!          kg/ha  = -------- * ------- * -------- * ---------
!                  cm2[soil]   cm[root]     m2         (g/m2)

      CumRootMass = CumRootMass + WRDOTN * 10. - SRDOT * 10. 
           
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE ROOTS_D
!=======================================================================


!=======================================================================
!  IPROOTD Subroutine
!  Reads root parameters from input files.
!----------------------------------------------------------------------
!  REVISION HISTORY
!  09/13/1998 CHP Written
C  08/12/2003 CHP Added I/O error checking
!-----------------------------------------------------------------------
!  Called : ROOTS
!  Calls  : FIND, ERROR, IGNORE
C=======================================================================
      SUBROUTINE IPROOTD(
     &  FILECC,                                                     !Input
     &  GSR, LWM, LWS, WCG, SFT, CAA, CAX, ALA, ALX,                !Output
     & PORMIN, RFAC1, RLDSM, RTDEPI, RTEXF,                         !Output
     &  RTSEN, RTSDF, RTWTMIN, XRTFAC, YRTFAC)                      !Output

!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT NONE

      CHARACTER*6 ERRKEY
      PARAMETER (ERRKEY = 'ROOTS')

      CHARACTER*6 SECTION
      CHARACTER*80 CHAR
      CHARACTER*92 FILECC
      CHARACTER*3   TYPRDT
      
      INTEGER LUNCRP, ERR, LNUM, ISECT, FOUND, II

      REAL RTDEPI, RLDSM, PORMIN
      REAL RFAC1, RTSEN, RTSDF, RTEXF
      REAL XRTFAC(4), YRTFAC(4), FNRDT(4)
      REAL GSR, LWM, LWS, WCG, SFT
      REAL CAA, CAX, ALA, ALX

!     Added 10/20/2005 for minimum root mass for senescence
      CHARACTER (len=7) RWMTXT
      REAL RTWTMIN

!-----------------------------------------------------------------------
!     ***** READ ROOT GROWTH PARAMETERS *****************
!-----------------------------------------------------------------------
!     Read in values from input file, which were previously input
!       in Subroutine IPCROP.
!-----------------------------------------------------------------------
      CALL GETLUN('FILEC', LUNCRP)
      OPEN (LUNCRP,FILE = FILECC, STATUS = 'OLD',IOSTAT=ERR)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,0)

!-----------------------------------------------------------------------
!    Find and Read Photosynthesis Section
!-----------------------------------------------------------------------
      SECTION = '!*ROOT'
      CALL FIND(LUNCRP, SECTION, LNUM, FOUND)
      IF (FOUND .EQ. 0) THEN
        CALL ERROR(SECTION, 42, FILECC, LNUM)
      ELSE
        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(5F6.0)',IOSTAT=ERR) RTDEPI,RFAC1,RTSEN,RLDSM,RTSDF
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(8F6.0)',IOSTAT=ERR)(XRTFAC(II),YRTFAC(II),II = 1,4)
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)

        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(12X,2F6.0)',IOSTAT=ERR) PORMIN, RTEXF
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
        
        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(5F6.0)',IOSTAT=ERR) GSR, LWM, LWS, WCG, SFT
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
        
        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(4F6.0)',IOSTAT=ERR) CAA, CAX, ALA, ALX
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
        
        CALL IGNORE(LUNCRP,LNUM,ISECT,CHAR)
        READ(CHAR,'(F6.0,T45,A7)',IOSTAT=ERR) RTWTMIN, RWMTXT
        IF (ERR .NE. 0 .OR. RWMTXT .NE. 'RTWTMIN') THEN
          RTWTMIN = 0.0
        ENDIF
        
        READ(CHAR,'(4F6.0,3X,A3)',IOSTAT=ERR)(FNRDT(II),II=1,4),TYPRDT
        IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
           
      ENDIF

      CLOSE (LUNCRP)

!***********************************************************************
      RETURN
      END SUBROUTINE IPROOTD
!=======================================================================

C=======================================================================
C  INROOTD Subroutine
C  Initializes root variables at emergence.
C----------------------------------------------------------------------
C  REVISION HISTORY
C  04/01/91 GH  Adapted for CROPGRO
C  06/17/98 CHP Modified for modular format
C  05/11/99 GH  Incorporated in CROPGRO
C-----------------------------------------------------------------------
C  Called : CROPGRO
C  Calls  : None
C=======================================================================
      SUBROUTINE INROOTD(
     &  DLAYR, FRRT, NLAYR, PLTPOP, RFAC1, RTDEPI, WTNEW, !Input
     &  RLV, RTDEP)                                       !Output

!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT NONE

      INTEGER L, NLAYR

      REAL DEP,RLINIT
      REAL RTDEP,RTDEPI,CUMDEP
      REAL RFAC1
      REAL WTNEW, FRRT, PLTPOP
      REAL RLV(NL), DLAYR(NL)

!***********************************************************************
C     INITIALIZE ROOT DEPTH AT EMERGENCE
C-----------------------------------------------------------------------
      RTDEP = RTDEPI
C-----------------------------------------------------------------------
C     DISTRIBUTE ROOT LENGTH EVENLY IN ALL LAYERS TO A DEPTH OF
C     RTDEPTI (ROOT DEPTH AT EMERGENCE)
C-----------------------------------------------------------------------
      CUMDEP = 0.

      DO L = 1,NLAYR
        RLV(L) = 0.0
      ENDDO

      DO L = 1,NLAYR
           DEP = MIN(RTDEP - CUMDEP, DLAYR(L))
           RLINIT = WTNEW * FRRT * PLTPOP * RFAC1 * DEP / ( RTDEP *
     &          10000 )
!        cm[root]   g[root]    plants  cm[root]   m2
!        -------- = -------- * ------ * ------- * ---
!      cm2[ground]   plant       m2     g[root]   cm2

           CUMDEP = CUMDEP + DEP
           RLV(L) = RLINIT / DLAYR(L)
           IF (CUMDEP .GE. RTDEP) GO TO 300
      ENDDO

  300 CONTINUE
!***********************************************************************
      RETURN
      END SUBROUTINE INROOTD
!=======================================================================
!=======================================================================
C=======================================================================
C  SCALC Subroutine
C  Computes static factors on rooting following Jones et al. (1991)
C----------------------------------------------------------------------
C  REVISION HISTORY
C  08/18/10 KJB placed in CROPGRO
C-----------------------------------------------------------------------
C  Called : CROPGRO
C  Calls  : None
C=======================================================================
      SUBROUTINE SCALC(NLAYR,DS,CAA,CAX,ALA,ALX,BD,CLAY,SILT,STONES,  !Input
     1  TOTBAS,EXTAL,EXCA,                                            !Input
     1  ZA,CWP,PO,CEC,ALS,SAL,CAS,SCA,BDO,BDX,SBD,SCF)        !Output

!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT NONE

!***********************************************************************
C
C  THIS SUBROUTINE CALCULATES THE STATIC (CONSTANT OVER SEASON)
C  STRESS FACTORS WHICH LIMIT ROOT GROWTH IN PARTICULAR SOIL LAYERS
C kjb
C      REAL CLA, SAN
C
C      DIMENSION CEC(NL)
C
!      REAL DAC, DMA, DMC, DMD, GSR, GSY, RD, RDX, RLL
!      REAL DAC, DMA, DMC, DMD, GSR, GSY, RD, RDX, RLL, RLV, RWL, Z
      INTEGER I, IJ, IR, L, NLAYR
!kjb
      REAL LWS
      REAL WCG,RTDEP
      REAL CAA,CAX,ALA,ALX,SFT
      REAL,DIMENSION(NL) :: BD,SAN,STONES,TOTBAS,EXTAL,EXCA,ALS,BDO,
     ! BDX,PO,CWP,SCA,SAL,SBD,SCF,ZA,Z_CM,CEC,CAS,DS,CLAY,SILT
!
C kjb
C      COMMON /BLK1/ CAA,CAX,ALA,ALX,BD,SAN,SIL,ROK,
C     1  SMB,EAL,CA,CLA,ALS,BDO,BDX
C      COMMON /BLK12/ PO,CWP,SCA,SAL,SBD,SCF
C      COMMON /BLK13/ Z,ZA
C      COMMON /BLK123/ IJ

!      PARAM   = "PARAM.TXT"
!      OutFile = "Roots.OUT"

!      OPEN(1,FILE=PARAM,ACCESS='SEQUENTIAL',STATUS='OLD')
!      OPEN(3,FILE=OUTFILE,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      
C  READ AND WRITE THE GENETIC PARAMETERS
C
!      READ (1,'(14F8.0)')RDX,GSR,LWM,LWS,WCG,TBS,TOP,CAA,CAX,ALA,ALX,PD,
!     1  GSD,SFT
!      WRITE (3,7000)
!      WRITE (3,5000) RDX,GSR,LWM,LWS,WCG,TBS,TOP,CAA,CAX,ALA,ALX,PD,
!     1  GSD,SFT
C
C      RD = PD
C      DO 30 I=1,NL
C       LWR(I) = LWS
C   30 CONTINUE
C
C  READ THE FOLLOWING SOIL PROPERTIES FOR UP TO NL SOIL LAYERS
C
!      READ (1,1000) TITLE
!      DO 100 I=1,NL
!	! CHP changed Z(I+1) to Z(I):
!      READ (1,'(11F8.0)',END=110) Z(I),BD(I),SAN(I),SIL(I),ROK(I),
!     1    SMB(I),EAL(I),CA(I),UL(I),LL(I),SCD(I)
!       IJ = IJ+1
!  100 CONTINUE  !bad here, as another 100 below
C
C  CALCULATE STATIC SOIL STRESS FACTORS
C
      DO 100 I=1,NLAYR
C
C  CALCULATE PERCENT OF SAND AND POROSITY OF EACH LAYER
C
      SAN(I) = 100.-CLAY(I)-SILT(I)
      PO(I) = ((100.-STONES(I))/100.)*(1.-BD(I)/2.65)
      CEC(I) = TOTBAS(I)+EXTAL(I)
C
C
C  CALCULATE DEPTH TO CENTER OF EACH LAYER
C
      IF (I.EQ.1) THEN
          ZA(I)= DS(I)/2
      ELSE
          ZA(I)=DS(I-1)+(DS(I)-DS(I-1))/2
      ENDIF
      
C
C
C  CALCULATE CRITICAL WATER-FILLED PORE FRACTION
C
      CWP(I)=.4+0.004*CLAY(I)
C
C  CALCULATE CALCIUM STRESS FACTOR BY LAYER
C
C kjb CAS(I) had to be made layer-wise, original code may have been wrong.
C as CAX and CAA, as well as ALX and ALA are constants.
C if EXCA(I) or CEC(I) or EAL(I) are missing then SCA(I)=1.0 and SAL(I)=1.0
      SCA(I)=1.0
      IF (EXCA(I).GT.0.0.AND.CEC(I).GT.0.0)  THEN
       CAS(I)=EXCA(I)/CEC(I)*100.
      IF (CAS(I).LT.CAA) SCA(I)=(CAS(I)-CAX)/(CAA-CAX)
      ENDIF
C
C  CALCULATE ALUMINUM TOXICITY STRESS FACTOR BY LAYER
C
      SAL(I)=1.
      IF (EXTAL(I).GT.0.0.AND.CEC(I).GT.0.0) THEN
      ALS(I)=EXTAL(I)/CEC(I)*100.
      IF (ALS(I).LE.ALA) THEN
        SAL(I)=1.
       ELSEIF (ALS(I).GE.ALX) THEN
        SAL (I)=0.
       ELSE
        SAL(I)=(ALX-ALS(I))/(ALX-ALA)
       ENDIF
      ENDIF
      
C
C  CALCULATE BULK DENSITY STRESS FACTOR BY LAYER
C
      BDO(I)=1.1+0.005*SAN(I)
      BDX(I)=1.6+0.004*SAN(I)
      IF (BD(I).LE.BDO(I)) THEN
        SBD(I)=1.
       ELSEIF (BD(I).GE.BDX(I)) THEN
        SBD(I)=0.
       ELSE
        SBD(I)=(BDX(I)-BD(I))/(BDX(I)-BDO(I))
      ENDIF
C
C  CALCULATE COARSE FRAGMENT STRESS FACTOR BY LAYER
C
      SCF(I)=(100.-STONES(I))/100.
C
  100 CONTINUE

      RETURN
C
C      END
C
C!***********************************************************************
C KJB
C      RETURN
      END SUBROUTINE SCALC
C************************************************
C
C************************************************
C
C
      SUBROUTINE DCALC(NLAYR,LL,DUL,LW,LT,TYPRDT,FNRDT,SFT,
     1  CWP,PO,SCA,SAL,SBD,SCF,SCD,              !inputs
     1  LWF,WFP,SST,SAI,STP,ASF,ASF1,ASF2,ASF3)  !outputs
C
C  THIS SUBROUTINE CALCULATES DYNAMIC STRESS FACTORS THROUGHOUT
C  THE SOIL PROFILE ON A DAILY BASIS
C
!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT NONE
      INTEGER I, NLAYR
      REAL SFT, FNRDT(4), CURV
      CHARACTER*3   TYPRDT
      REAL,DIMENSION(NL) :: 
     1  LL,DUL,LW,LT,
     1  CWP,PO,SCA,SAL,SBD,SCF,SCD,                 
     1  LWF,WFP,SST,SAI,STP,ASF,ASF1,ASF2,ASF3
C
C      COMMON /BLK2/ LL(10),UL(10),LT(10),SST(10),LW(10),STP(10),
C     1  SCD(10),TBS,TOP,SFT
C      COMMON /BLK12/ PO(10),CWP(10),SCA(10),SAL(10),SBD(10),SCF(10)
C      COMMON /BLK23/ LWF(10),SAI(10),ASF(10),ASF1(10),ASF2(10),ASF3(10)
C      COMMON /BLK123/ IJ
C
C  CALCULATE DYNAMIC ROOT GROWTH STRESS FACTORS FOR SOIL LAYERS
C
      DO 100 I=1,NLAYR
C
C  CALCULATE LAYER STRENGTH FACTORS
C
      IF (LW(I).LT.LL(I)) THEN
        LWF(I)=0.
       ELSEIF (LW(I).GT.DUL(I)) THEN
        LWF(I)=1.0
       ELSE
        LWF(I)=(LW(I)-LL(I))/(DUL(I)-LL(I))
      ENDIF
      SST(I)=SBD(I)*SIN(1.57*LWF(I))
C
C  CALCULATE LAYER AERATION FACTORS
C
      WFP(I)=LW(I)/PO(I)
      SAI(I)=1.
      IF (WFP(I).GE.CWP(I)) THEN
        SAI(I)=SFT+(1-WFP(I))*((1-SFT)/(1-CWP(I)))
       ELSE
        SAI(I)=1.0
      ENDIF
      IF(SAI(I).LT.0.) SAI(I)=0.
C
C  CALCULATE LAYER TEMPERATURE FACTORS
C
      !IF(LT(I).GE.TBS)THEN
      !  STP(I)=SIN(1.5707*(LT(I)-TBS)/(TOP-TBS))
      ! ELSE
      !  STP(I)=0.
      ! ENDIF
       
       STP(I) = CURV(TYPRDT,FNRDT(1),FNRDT(2),FNRDT(3),FNRDT(4),LT(I))
C
C  CALCULATE MINIMUM OF STATIC AND DYNAMIC STRESS FACTORS
C
      ASF(I)=MIN(STP(I),SST(I),SAI(I),SCA(I),SAL(I),SCF(I),SCD(I))
      ASF1(I)=MIN(STP(I),SCA(I),SAL(I),SCD(I))
      ASF2(I)=(MIN(SST(I),SAI(I),SCF(I)))**.5
      ASF3(I)=MIN(SST(I),SAL(I),SCF(I))
      IF(ASF(I).LT.0.) ASF(I)=10E-10
      IF(ASF1(I).LT.0.) ASF1(I)=10E-10
      IF(ASF2(I).LT.0.) ASF2(I)=10E-10
      IF(ASF3(I).LT.0.) ASF3(I)=10E-10
C
  100 CONTINUE
C
      RETURN
C
      END SUBROUTINE DCALC
C
C*********************************************************
C*********************************************************
C    PD SAME AS RTDEPI
C
      SUBROUTINE RTDIST(NLAYR,DS,ZA,GSR,GSY,GS, WRDOTN,    !inputs
     1  WCG,LWS,LWM,DMD,ASF,ASF1,ASF2,ASF3,SAI,LWF, RFAC2, DTX,      !inputs
     1  DDI,DRD,RTDEP,WFL,DZ,WFT,LWN,LWA,RLV,TEMP,GPL,GPS,  !outputs
     1  GAL,RWL,DMA,DWL,DLL,RLL,LWR,TRW,DMC,DAC)         !outputs
C
C  THIS SUBROUTINE DISTRIBUTES ROOT GROWTH THROUGHOUT THE SOIL
C  PROFILE IN RESPONSE TO ROOT SYSTEM GROWTH, LAYER DEPTH,
C  STATIC STRESS FACTORS, AND DYNAMIC STRESS FACTORS
!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT NONE
C
      REAL,DIMENSION(NL) :: WFL,DZ,GPL,LWA,LWR,DWL,
     1  RWL,DLL,RLL,RLV,GAL,DS,ZA,
     1  LWN,LWF,SAI,ASF,ASF1,ASF2,ASF3
      REAL LWS,LWM,DDI,DRD,WFT,TEMP,GPS,TRW
      REAL WCG,GS,GSR,GSY, RFAC2, DTX, WRDOTN
      REAL DM, DAC, DMA, DMC, DMD, RTDEP
      INTEGER I, IR, NLAYR
!kjb


C
C      COMMON /BLK3/ GAL(10),DWL(10),RLL(10),RWL(10),RLV(10),LWR(10),
C     1  RDX,GSR,LWM,LWS,WCG,GSD,GS,DMD,GSY,DMC,DMA,TRW,RD,
C     2  IR,DAC
C      COMMON /BLK13/ Z(11),ZA(10)
C      COMMON /BLK23/ LWF(10),SAI(10),ASF(10),ASF1(10),ASF2(10),ASF3(10)
C      COMMON /BLK123/ IJ
C
C  INITIALIZE VARIABLES
C
      DMA=0.
      GPS=0.
      TRW=0.
      WFT=0.
C      RD = PD
C      DMC = 0.0
C
C  DETERMINE IF ROOT SYSTEM IS STILL GROWING DOWNWARD
C
C      GS = VSTAGE / 19.   !Divide by max vstage for 0-1 fraction

      IF(GS.LT.GSR.OR.WRDOTN.GT.0.0) THEN
C
C  DETERMINE POTENTIAL INCREASE IN ROOT DEPTH
C
C      RD = RTDEP
C      DDI=RDX*(GS-GSY)/GSR
      DDI =  RFAC2
      DO 100 I=1,NLAYR
       IF (RTDEP.LT.DS(I)) THEN
C          RTDEP = RTDEP + DTX * RFAC2 * MIN(SWDF,SWEXF) *
C     &                (1. + 0.25 * (1. - MAX(SWFAC,0.40)))
C       DRD=DDI * MIN(ASF1(I),ASF2(I))
        DRD=DDI * MIN(ASF1(I),ASF2(I))
        IF(DRD.GT.(DS(I)-RTDEP)) THEN
          DDI=DDI-(DS(I)-RTDEP)
          DRD = DS(I)-RTDEP
          RTDEP = RTDEP +DRD
C          RD=Z(I+1)
         ELSE
           RTDEP=RTDEP+DRD
           GO TO 200
        ENDIF
       ENDIF
  100 CONTINUE
C
C  DETERMINE NUMBER OF LAYERS CURRENTLY CONTAINING ROOTS
C
C      DO 222 I=1,NLAYR
  200   IF(I.LT.NLAYR)THEN
        IR=I
         ELSE
        IR=NLAYR
         ENDIF
C  222 CONTINUE  
C
C this setS IR equal NLAYR so OK
         ENDIF
  300 DO 400 I=1,IR
C
C  CALCULATE WEIGHTING FACTOR FOR GROWTH HABIT
C
      WFL(I)=(1.-ZA(I)/300.0)**WCG
      IF (I.EQ.1) THEN
          DZ(I)=DS(I)
      ELSE
          DZ(I)=DS(I)-DS(I-1)
      ENDIF
      
      WFT=WFT+WFL(I)
C
C  CALCULATE NORMAL LENGTH/WEIGHT RATIO AS AFFECTED BY GROWTH
C  STAGE AND LAYER DEPTH
C
      IF(ZA(I).GT.RTDEP)THEN
        LWN(I)=LWS
      ELSE
Ckjb, was just LWN.  I made LWM(I). MAYBE ERROR IN PRIOR CODE.
        LWN(I)=LWS-(GS/GSR)*(LWS-LWM)*(1.-ZA(I)/RTDEP)
      ENDIF
C
C  CALCULATE ACTUAL ROOT LENGTH/WEIGHT RATIOS AS AFFECTED BY
C  ALUMINUM, STRENGTH OR COARSE FRAGMENTS
C
      LWA(I)=LWN(I)/(1.+3.*(1.-ASF3(I)))
C
  400 CONTINUE
C
C  DETERMINE POTENTIAL FOR ROOT GROWTH IN EACH LAYER
C
C kjb The first part is just to create some RLV in a layer.  Modify later?
C The second part (TEMP) is a driver of increased RLV that assumes rapid growth
C in a layer requires a certain amount of existing root presence.  Reasonable.
C All of this along with ASF, depth factor, and present cm/g.  Math ends up with g
C Hopefully g/m2
C
C WFL is the weighting of assimilate to layers (sums to 1 over all layers)
      DO 500 I=1,IR
       IF(RLV(I).EQ.0.0)THEN
C         RLV(I)=0.001*((RD-Z(I))/DZ(I))
         IF(I.EQ.1) THEN
             RLV(I)=0.01*(RTDEP/DZ(I))
         ELSE
             RLV(I)=0.01*((RTDEP-DS(I-1))/DZ(I))
         ENDIF 
       ENDIF
       WFL(I)=WFL(I)/WFT
       TEMP= (5.0*RLV(I))/(0.025+RLV(I))
!       GPL(I)=(TEMP*ASF(I)*WFL(I)*DZ(I)*1E4)/LWA(I)
       GPL(I)=(TEMP*ASF(I)*WFL(I)*DZ(I)*1E2)/LWA(I)
       GPS=GPS+GPL(I)
  500 CONTINUE
C
C
      DO 600 I=1,IR
C
C  DISTRIBUTE ROOT GROWTH BY LAYERS
C   DMD (g) divided by LWA (m/g) leaves you with GAL in m.  GPL/GPS distributes to layers
C
      IF(GPS.GT.DMD)THEN
        GAL(I)=DMD*(GPL(I)/GPS)*LWA(I)
      ELSE
C kjb, the next better not happen, otherwise, dry matter to root does not balance
          GAL(I)=GPL(I)*LWA(I)
       ENDIF
C Note:  Math indicates that GAL is meters, and RWL is g (presumed g/m2)
C DMA is "new" root weight over all layers.  Should match DMD.
       RWL(I)=RWL(I)+(GAL(I)/LWA(I))
      DMA=DMA+GAL(I)/LWA(I)
C
C  CALCULATE ROOT DEATH BY LAYER, DWL(I) is today's root death (g) by layer
C  kjb - Do not see why LWF (distribution weighting should act here
C
C Maybe do senescence in main root or pass info to here? Do not like the GS and GSD
      DWL(I)=0.01*RWL(I)*(1+MAX((1.-LWF(I)),(1.-SAI(I))))
      IF(GS.GT.GSR)THEN
        DWL(I)=DWL(I)+RWL(I)*((GS-GSR)/(1.-GSR))**3.
      ENDIF
C RWL is actual root mass (g) in a layer, subtracting dead roots (g)
C DLL is today's dead root mass (g) converted to m
C RLL is actual root length (m) in a layer, where GAL is new root (m) and DLL is today's dead root
C RLV is root length density, where RLL is m, divided by layer thickness in m
C LWR is re-computed root length per unit mass (m/g) from root length (m) divided by root mass (g)
C TRW is summation of root mass over all layers
      RWL(I)=RWL(I)-DWL(I)
      DLL(I)=DWL(I)*LWR(I)
      RLL(I)=RLL(I)+GAL(I)-DLL(I)
      RLV(I)=RLL(I)/(DZ(I)*1E2)
      LWR(I)=RLL(I)/(RWL(I)+10E-10)
      TRW=TRW+RWL(I)
C
  600 CONTINUE
C
      GSY=GS
      DMC=DMC+DMD
C kjb deleted, as DM does not exist in RTDIST.  What is DAC supposed to be?
C      DAC=DAC+DM
C
!     Possible missing code here (OLD FROM ORIGINAL?, KJB)
      return
      END   SUBROUTINE RTDIST
C
CC
!-----------------------------------------------------------------------
!       Variable definitions
!-----------------------------------------------------------------------
! AGRRT     Mass of CH2O required for new root growth (g[CH2O] / g[root])
! CGRRT     Carbon demand for new root growth (g[CH2O] / m2 / d)
! CROP      Crop identification code 
! CUMDEP    Cumulative depth of soil profile (cm)
! DEP       Cumulative soil depth (cm)
! DEPMAX    Maximum depth of reported soil layers (cm)
! DLAYR(L)  Soil thickness in layer L (cm)
! DS(L)     Cumulative depth in soil layer L (cm)
! DTX       Thermal time that occurs in a real day based on vegetative 
!             development temperature function (thermal days / day)
! DUL(L)    Volumetric soil water content at Drained Upper Limit in soil 
!             layer L (cm3 [H2O] /cm3 [soil])
! ESW(L)    Plant extractable soil water by layer (= DUL - LL) (cm3/cm3)
! FILECC    Path plus filename for species file (*.spe) 
! FRRT      Fraction of vegetative tissue growth that goes to roots on a 
!             day (g[root] / g[veg])
! GRESPR(L) Growth respiration for new root growth in layer L 
! LL(L)     Volumetric soil water content in soil layer L at lower limit
!             ( cm3/cm3)
! LUNCRP    Logical unit number for FILEC (*.spe file) 
! LUNIO     Logical unit number for FILEIO 
! MRESPR(L) Maintenance respiration for new root growth in layer L 
! NL        Maximum number of soil layers = 20 
! NLAYR     Number of soil layers 
! PG        Daily gross photosynthesis (g[CH2O] / m2 / d)
! PLTPOP    Plant population (# plants / m2)
! PORMIN    Minimum pore space required for supplying oxygen to roots for 
!             optimal growth and function (cm3/cm3)
! RESPS(L)  Total respiration for new root growth in layer L 
! RFAC1     Root length per unit  root weight. (cm/g)
! RFAC2     Root depth increase rate with time (cm / physiol. day)
! RFAC3     Ratio of root length to root weight at the current time (cm/g)
! RLDF(L)   Combined weighting factor to determine root growth distribution
! RLDSM     Minimum root length density in a given layer, below which 
!             drought-induced senescence is not allowed.
!             (cm [root ]/ cm3 [soil])
! RLGRW(L)  Incremental root length density in soil layer L
!             (cm[root] / cm3[soil])
! RLINIT    Initial root density (cm[root]/cm2[ground])
! RLNEW     New root growth added (cm[root]/cm2[ground]/d)
! RLSEN(L)  Root length density senesced today (cm[root]/ cm3[soil])
! RLV(L)    Root length density for soil layer L (cm[root] / cm3[soil])
! RO        Respiration coefficient that depends on total plant mass
!             (g[CH2O] / g[tissue])
! RP        proportion of the day's photosynthesis which is respired in the 
!             maintenance process 
! RRLF(L)   Root length density factor ratio (RLDF(L) / TRLDF) 
! RTDEP     Root depth (cm)
! RTDEPI    Depth of roots on day of plant emergence. (cm)
! RTEXF     Fraction root death per day under oxygen depleted soil 
! RTSDF     Maximum fraction of root length senesced in a given layer per 
!             physiological day when water content in a given layer falls 
!             below 25 % of extractable soil water. 
! RTSEN     Fraction of existing root length which can be senesced per 
!             physiological day. (fraction / ptd)
! RTSURV(L) Fraction survival of roots on a given day, taking into account 
!             death due to excess or deficit water conditions 
! RTWT      Dry mass of root tissue, including C and N
!             (g[root] / m2[ground])
! SAT(L)    Volumetric soil water content in layer L at saturation
!             (cm3 [water] / cm3 [soil])
! SRDOT     Daily root senescence (g / m2 / d)
! SW(L)     Volumetric soil water content in layer L
!             (cm3 [water] / cm3 [soil])
! SWDF      Soil water deficit factor for layer with deepest roots (0-1) 
! SWEXF     Excess water stress factor for layer with deepest roots (0-1) 
! SWFAC     Effect of soil-water stress on photosynthesis, 1.0=no stress, 
!             0.0=max stress 
! TABEX     Function subroutine - Lookup utility 
! TRLDF     Total root length density factor for root depth (cm)
! TRLGRW    Total new root length density in soil layer L
!             (cm[root] / cm2[soil])
! TRLSEN    Total root length density senesced today (cm[root]/ cm2[soil])
! TRLV      Total root length per square cm soil today 
!             (cm[root]/cm2[soil])
! TRTDY     Total root length per square cm soil yesterday 
!             (cm[root]/cm2[soil])
! VSTAGE    Number of nodes on main stem of plant 
! WR(L)     Root hospitality factor, used to compute root growth 
! WRDOTN    Dry weight growth rate of new root tissue including N but not C 
!             reserves (g[root] / m2[ground]-d)
! WTNEW     Initial mass of seedling or seed (g / plant)
! XRTFAC(I) V-stage at which rate of increase in root depth per 
!             physiological day is YRTFAC(I). (# leaf nodes)
! YRTFAC(I) Rate of increase in root depth per degree day at V-stage 
!             XRTFAC(I). (cm / (physiol. day))
!***********************************************************************
!      END SUBROUTINES ROOTS, IPROOT, and INROOT
!=======================================================================
