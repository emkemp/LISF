!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
!BOP
!
! !MODULE: APHROPRCP_obsMod
! \label(APHROPRCP_obsMod)
!
! !INTERFACE:
module APHROPRCP_obsMod
!
! !USES:
  use ESMF

  implicit none
  PRIVATE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This module handles the observation plugin for the
!  APHRODITE (Asian Precipitaton - Highly-Resolved
!  Observational Data Integration Towards Evaluation)
!  daily gridded precipitation data.
!
!  https://climatedataguide.ucar.edu/climate-data/aphrodite-asian-precipitation-highly-resolved-observational-data-integration-towards
!
!  Coverage is from 1951 - 2007, daily at
!  0.25 deg lat/lon resolution
!
! !FILES USED:
!
! !REVISION HISTORY:
!  09 Aug 2017   Sujay Kumar  Initial Specification
!  14 Sep 2020   Eric Kemp    Added V1901 netCDF support, valid 1998-2015.
!                             Updated interpolation.
!  15 Sep 2020   Eric Kemp    Significant changes to mimic CHIRPSv2 reader.
!EOP

  PUBLIC :: APHROPRCP_obsinit
  PUBLIC :: APHROPRCPobs

  type, public :: aphroprcpobsdec
     character*100        :: odir
     character*100        :: loc
     character*100        :: version ! EMK
     integer              :: nc, nr

     real,    allocatable     :: rlat(:)
     real,    allocatable     :: rlon(:)

     ! EMK Used for upscale averaging
     integer, allocatable        :: n11(:)

     ! Used with budget interpolation (precip)
     integer, allocatable        :: n112(:,:)
     integer, allocatable        :: n122(:,:)
     integer, allocatable        :: n212(:,:)
     integer, allocatable        :: n222(:,:)
     real,    allocatable        :: w112(:,:)
     real,    allocatable        :: w122(:,:)
     real,    allocatable        :: w212(:,:)
     real,    allocatable        :: w222(:,:)

     real :: datares

  end type aphroprcpobsdec

  type(aphroprcpobsdec), allocatable :: aphroprcpobs(:)

contains

!BOP
!
! !ROUTINE: APHROPRCP_obsInit
! \label{APHROPRCP_obsInit}
!
! !INTERFACE:
  subroutine APHROPRCP_obsinit(i)
!
! !USES:
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS:
    integer,   intent(IN) :: i
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine initializes and sets up the data structures required
!  for reading APHRO PCP data.
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP

    real               :: gridDesci(50)
    integer            :: status

    if(.not.allocated(aphroprcpobs)) then
       allocate(aphroprcpobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, aphroprcpobs(i)%odir, &
         label='APHRO PCP data directory:', rc=status)
    call LVT_verify(status, 'APHRO PCP data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, aphroprcpobs(i)%loc, &
         label='APHRO PCP data region:', rc=status)
    if(status.ne.0) then
       write(LVT_logunit,*) '[ERR] APHRO PCP data region: not defined'
       write(LVT_logunit,*) '[ERR] options are: '
       write(LVT_logunit,*) "[ERR] 'MA' (for Monsoon Asia)' "
       write(LVT_logunit,*) "[ERR] 'ME' (for Middle East)' "
       write(LVT_logunit,*) "[ERR] 'RU' (for Northern Eurasia)' "
       write(LVT_logunit,*) "[ERR] 'PR' (for Combined Eurasia)' "
       call LVT_endrun()
    endif

    ! EMK Add version number
    call ESMF_ConfigGetAttribute(LVT_config, aphroprcpobs(i)%version, &
         label='APHRO PCP version:', rc=status)
    if(status.ne.0) then
       write(LVT_logunit,*) '[ERR] APHRO PCP version: not defined'
       write(LVT_logunit,*) '[ERR] options are: '
       write(LVT_logunit,*) "[ERR] 'V1101' Original Release (1951-2007)"
       write(LVT_logunit,*) "[ERR] 'V1901' (1998-2015) "
       call LVT_endrun()
    endif

    gridDesci = 0
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(aphroprcpobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    if(aphroprcpobs(i)%loc.eq."MA") then
       aphroprcpobs(i)%nc = 360
       aphroprcpobs(i)%nr = 280

       aphroprcpobs(i)%datares = 0.25 ! deg

       gridDesci(1) = 0
       gridDesci(2) = 360
       gridDesci(3) = 280
       gridDesci(4) = -14.875
       gridDesci(5) = 60.125
       gridDesci(6) = 128
       gridDesci(7) = 54.875
       gridDesci(8) = 149.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64

       ! EMK Use budget-bilinear interpolation if APHRODITE is at coarser
       ! resolution than the analysis grid; otherwise, use upscale averaging.
       if (LVT_isAtAFinerResolution(aphroprcpobs(i)%datares)) then

          allocate(aphroprcpobs(i)%n112(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%n122(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%n212(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%n222(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%w112(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%w122(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%w212(LVT_rc%lnc*LVT_rc%lnr, 25))
          allocate(aphroprcpobs(i)%w222(LVT_rc%lnc*LVT_rc%lnr, 25))

          aphroprcpobs(i)%n112 = 0
          aphroprcpobs(i)%n122 = 0
          aphroprcpobs(i)%n212 = 0
          aphroprcpobs(i)%n222 = 0
          aphroprcpobs(i)%w112 = 0
          aphroprcpobs(i)%w122 = 0
          aphroprcpobs(i)%w212 = 0
          aphroprcpobs(i)%w222 = 0

          call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
               LVT_rc%lnc*LVT_rc%lnr, &
               aphroprcpobs(i)%rlat, aphroprcpobs(i)%rlon,&
               aphroprcpobs(i)%n112, aphroprcpobs(i)%n122, &
               aphroprcpobs(i)%n212, aphroprcpobs(i)%n222, &
               aphroprcpobs(i)%w112, aphroprcpobs(i)%w122, &
               aphroprcpobs(i)%w212, aphroprcpobs(i)%w222)

       else

          allocate(aphroprcpobs(i)%n11(aphroprcpobs(i)%nc*aphroprcpobs(i)%nr))
          aphroprcpobs(i)%n11 = 0

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc,&
               aphroprcpobs(i)%nc*aphroprcpobs(i)%nr, &
               LVT_rc%lnc*LVT_rc%lnr, &
               aphroprcpobs(i)%n11)
       end if
    else
       write(LVT_logunit,*) "[ERR] The Aphrodite plugin only supports"
       write(LVT_logunit,*) "[ERR] the MA region currently"
       call LVT_endrun()
    endif

  end subroutine APHROPRCP_obsinit

end module APHROPRCP_obsMod
