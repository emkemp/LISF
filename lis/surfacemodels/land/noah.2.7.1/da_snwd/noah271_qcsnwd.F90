!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah271_qcsnwd
! \label{noah271_qcsnwd}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_qcsnwd(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah271_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  integer                :: t
  integer                :: status
  real, allocatable          :: swe(:)
  real, allocatable          :: snod(:)

  real                   :: swemax,snodmax
  real                   :: swemin,snodmin
 
!  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
!  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)
 
!  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
!  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

!  call ESMF_AttributeGet(sweField,"Max Value",swemax,rc=status)
!  call LIS_verify(status)
!  call ESMF_AttributeGet(sweField,"Min Value",swemin,rc=status)
!  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Max Value",snodmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Min Value",snodmin,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
!     if(swe(t).lt.swemin) then 
!        swe(t) = swemin
!     endif
     if(snod(t).lt.snodmin) then 
        snod(t) = snodmin
     endif
!Additional qc's to ensure consistency.. 
!     if(swe(t).le.0.0) snod(t) = 0.0
!     if(snod(t).lt.swe(t)) then 
!        snod(t) = swe(t)
!     endif
  enddo
end subroutine noah271_qcsnwd

