!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noah33_dasnow_Mod
!BOP
!
! !MODULE: noah33_dasnow_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah33_dasnow_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: noah33_dasnow_init
! \label{noah33_dasnow_init}
! 
! !INTERFACE:
  subroutine noah33_dasnow_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
  end subroutine noah33_dasnow_init
end module noah33_dasnow_Mod