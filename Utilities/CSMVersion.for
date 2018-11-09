!=======================================================================
!  MODULE CSMVersion
!  07/05/2018 CV  Written
!=======================================================================

      MODULE CSMVersion
!     Contains defintion of version and branch which will be
!     used throughout the model.

!=======================================================================

!     Global CSM Version Number
      TYPE VersionType
        INTEGER :: Major = 4
        INTEGER :: Minor = 7
        INTEGER :: Model = 2
        INTEGER :: Build = 3
      END TYPE VersionType
      TYPE (VersionType) Version
      CHARACTER(len=*), PARAMETER :: VBranch = '-ASCE_kc_chp'
!     Implement Kelly Thorp's ASCE PET routines for soybean
!     GitHash: 6928dc31442e29f6c1b00f6f612818362b57a951

!     Version history:  
!       4.7.2.3  chp 10/25/2018 Rice - G5 parameter   
!       4.7.2.2  chp 10/25/2018 Sugarbeet, Forage, Cassava,    
!                               ASCE ET added (on hold),
!                               Fixed bug with maize byproduct harvest percentage 
!       4.7.2.1  chp 07/08/2018 Use CSMVersion, minor bug fixes, Linux hacks, 
!                                 VSCode files
!       4.7.2.0  chp 05/07/2018 v4.7.2 Release 2018 Workshop
!       4.7.1.0  chp 10/27/2017 v4.7.1 Release
!       4.7.0.0  chp 08/09/2017 v4.7.0 Release
!       4.6.5.1  chp 05/10/2017 v4.6.5 Release 2017 Workshop  
!       4.6.0.1  chp 06/28/2011 v4.6.0 Release
!       4.5.1.0  chp 10/10/2010 V4.5.1 Release
!       4.0.2.0  chp 08/11/2005 v4.0.2 Release
!       4.0.1.0  chp 01/28/2004 v4.0.1 Release 

!=======================================================================
      END MODULE CSMVersion

