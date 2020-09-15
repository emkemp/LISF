!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: readAPHROPRCPobs
! \label{readAPHROPRCPobs}
!
! !INTERFACE:
subroutine readAPHROPRCPobs(source)
!
! !USES:
   use ESMF
   use LVT_coreMod,    only : LVT_rc, LVT_domain, LVT_isAtAfinerResolution
   use LVT_histDataMod
   use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
        LVT_releaseUnitNumber, LVT_verify, LVT_endrun
   use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
   use APHROPRCP_obsMod,    only : aphroprcpobs
   use map_utils
   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   
   implicit none
!
! !INPUT PARAMETERS:
   integer,  intent(in) :: source
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine reads and processes the gridded APHRO precip
!  data. The data for a given year is read into memory at the
!  start of a year and indexed into during each day.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  09 AUG 2017: Sujay Kumar, Initial Specification
!  14 Sep 2020: Eric Kemp, Added V1901 support.  Updated interpolation, and
!                 added consistent check for missing data.
!  15 Sep 2020: Eric Kemp.  Significant rewrite.  Now mimics CHIRPSv2 reader.
!EOP

   character*100       :: filename
   character*4         :: fyr
   logical             :: file_exists
   real                :: &
        rainf(aphroprcpobs(source)%nc, aphroprcpobs(source)%nr)
   real                :: &
        prcp_in(aphroprcpobs(source)%nc*aphroprcpobs(source)%nr)
   logical*1           :: lb(aphroprcpobs(source)%nc*aphroprcpobs(source)%nr)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp_out(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp_final(LVT_rc%lnc, LVT_rc%lnr)
   integer             :: c, r
   integer             :: nid, varid
   integer             :: iret
   logical             :: alarmCheck
   real                :: fill_value
   integer             :: yr1, mo1, da1, hr1, mn1, ss1
   integer             :: yr2, mo2, da2, hr2, mn2, ss2
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   integer :: start(3), count(3)
   integer :: jda2
   real                :: currTime

   ! Initialize variables
   prcp_out = LVT_rc%udef
   prcp_final = LVT_rc%udef

   currTime = float(LVT_rc%dhr(source))*3600+ &
        60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   ! EMK...Only read at 00Z
   alarmCheck = (mod(currtime,86400.0).eq.0)

   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = 0
   ss1 = 0
   call ESMF_TimeSet(time1, yy=yr1, mm=mo1, dd=da1, &
        h=hr1, m=mn1, s=ss1, calendar=LVT_calendar, rc=iret)
   call LVT_verify(iret)

   ! If it is 00Z, use previous day's time level
   if (mod(currtime, 86400.0).eq.0) then
      call ESMF_TimeIntervalSet(lis_ts, s=86400, &
           rc=iret)
      call LVT_verify(iret)
   else
      call ESMF_TimeIntervalSet(lis_ts, s=0, &
           rc=iret)
      call LVT_verify(iret)
   end if
   time2 = time1 - lis_ts

   call ESMF_TimeGet(time2, yy=yr2, mm=mo2, dd=da2, &
        h=hr2, m=mn2, s=ss2, calendar=LVT_calendar, &
        dayOfYear=jda2, rc=iret)
   call LVT_verify(iret)

   if (alarmCheck) then

      ! Create filename
      ! FIXME...Put in separate subroutine
      if (aphroprcpobs(source)%loc .eq. "MA") then
         write(fyr, fmt='(i4.4)') yr2
         ! EMK...Add support for V1901
         if (trim(aphroprcpobs(source)%version) .eq. "V1901") then
            filename = trim(aphroprcpobs(source)%odir)//&
                 '/APHRO_MA_025deg_V1901.'//trim(fyr)//'.nc'
         else
            if( LVT_rc%dyr(source) >= 2007) then           ! Yeosang Yoon
               filename = trim(aphroprcpobs(source)%odir)//&
                    '/APHRO_MA_025deg_V1101_EXR1.'//trim(fyr)//'.nc'
            else
               filename = trim(aphroprcpobs(source)%odir)//&
                    '/APHRO_MA_025deg_V1101.'//trim(fyr)//'.nc'
            end if
         end if
      end if

      inquire(file=trim(filename), exist=file_exists)

      if (file_exists) then
         write(LVT_logunit,*) '[INFO] Reading APHRO data ', trim(filename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
         iret = nf90_open(path=trim(filename), mode=NF90_NOWRITE, ncid=nid)
         call LVT_verify(iret, 'Error opening file'//trim(filename))

         iret = nf90_inq_varid(nid, 'precip', varid)
         call LVT_verify(iret, 'Error nf90_inq_varid: precip')

         start(1) = 1
         start(2) = 1
         start(3) = jda2
         count(1) = aphroprcpobs(source)%nc
         count(2) = aphroprcpobs(source)%nr
         count(3) = 1

         iret = nf90_get_var(nid, varid, rainf, start=start, &
              count=count)
         call LVT_verify(iret, 'Error nf90_get_var: precip')

         ! EMK Get missing value from file
         iret = nf90_get_att(nid, varid, "_FillValue", fill_value)
         call LVT_verify(iret, 'Error nf90_get_att: precip')

         iret = nf90_close(nid)
         call LVT_verify(iret, 'Error nf90_close')
#endif

         ! EMK Make sure consistent missing values are used.
         do r = 1, aphroprcpobs(source)%nr
            do c = 1, aphroprcpobs(source)%nc
               if (rainf(c, r) .eq. fill_value) then
                  rainf(c, r) = LVT_rc%udef
               end if
            end do
         end do

         lb = .false.
         prcp_in = LVT_rc%udef
         do r = 1, aphroprcpobs(source)%nr
            do c = 1, aphroprcpobs(source)%nc
               if (rainf(c, r) .ge. 0) then
                  prcp_in(c+(r-1)*aphroprcpobs(source)%nc) = &
                       rainf(c, r)
                  lb(c+(r-1)*aphroprcpobs(source)%nc) = .true.
               endif
            enddo
         enddo

         ! Interpolate
         if (LVT_isAtAFinerResolution(aphroprcpobs(source)%datares)) then
           call conserv_interp(LVT_rc%gridDesc, lb, prcp_in, &
                lo, prcp_out, &
                aphroprcpobs(source)%nc * aphroprcpobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr, &
                aphroprcpobs(source)%rlat, aphroprcpobs(source)%rlon, &
                aphroprcpobs(source)%w112, aphroprcpobs(source)%w122, &
                aphroprcpobs(source)%w212, aphroprcpobs(source)%w222, &
                aphroprcpobs(source)%n112, aphroprcpobs(source)%n122, &
                aphroprcpobs(source)%n212, aphroprcpobs(source)%n222, &
                LVT_rc%udef, iret)
        else
           call upscaleByAveraging(&
                aphroprcpobs(source)%nc * aphroprcpobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
                aphroprcpobs(source)%n11, lb, &
                prcp_in, lo, prcp_out)
        end if

     endif ! File exists

     do r = 1, LVT_rc%lnr
        do c = 1, LVT_rc%lnc
           prcp_final(c, r) = prcp_out(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo

     write(LVT_logunit,*) '[INFO] Finished processing ', trim(filename)

  end if ! alarmCheck

  ! Convert mm/day to kg/m2s (note that 1 mm = 1 kg/m2)
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp_final(c,r).ge.0) then
           prcp_final(c,r) = prcp_final(c,r)/86400.0 !kg/m2s
        else
           prcp_final(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf, source, prcp_final, &
       vlevel=1, units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc, source, prcp_final, &
       vlevel=1, units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, prcp_final, &
       vlevel=1, units='kg/m2s')

  ! Now convert from kg/m2s to kg/m2
  do r = 1, LVT_rc%lnr
     do c = 1, LVT_rc%lnc
        if (prcp_final(c, r) .ge. 0) then
           prcp_final(c, r) = prcp_final(c, r)*86400.0 !kg/m2
        else
           prcp_final(c, r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf, source, prcp_final, &
       vlevel=1, units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc, source, prcp_final, &
       vlevel=1, units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, prcp_final, &
       vlevel=1, units='kg/m2')

end subroutine readAPHROPRCPobs


