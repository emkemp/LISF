!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_MERIT1K_slope
! \label{read_MERIT1K_slope}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 Aug  2012: KR Arsenault; Expanded for elevation tiling
!  30 May  2017: KR Arsenault; Expanded for Antarctica
!  03 Mar  2020: Yeosang Yoon; Modify codes for MERIT DEM
!  04 Apr  2024: Yeosang Yoon; Fix bug slope value along the coastal area
!
! !INTERFACE:
subroutine read_MERIT1K_slope( n, num_bins, fgrd, slopeave )

! !USES:
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use calc_SlopeAspect_module, only: calc_Slope_fromElev
  use LDT_paramTileInputMod,   only: param_1dbin_areacalc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n  
  integer, intent(in) :: num_bins
  real,    intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out):: slopeave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves static elevation data from the MERIT source
!  and reprojects it to the latlon projection. 
!
!  Source information:
!   http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/index.html
!   Native MERIT DEM has a resolution '0.0008333333' (~ 90 m). To consider
!   computational complexity and current usage of LIS requirement, LIS team
!   upscales the DEM to have a resolution '0.008333' (~ 1 km) and this
!   reader handle the upscaled DEM.
!   Upscaled method: Area averageing technique with 25% cutoff for land
!   determination.
!!
!  For the Antarctica files, currently use the GTOPO30 tiled files:
!   https://lta.cr.usgs.gov/GTOPO30
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[fgrd]
!   output grid fractions for elevation bins (or bands)
!  \item[slopeave]
!   elevation average for each bin (or band)
!  \end{description}
!EOP

   integer, parameter :: input_cols = 360*120
   integer, parameter :: input_rows = 180*120   ! MERIT elev tiles span: 90N to -90S
   real,    parameter :: input_xres = 1.0/120.0
   real,    parameter :: input_yres = 1.0/120.0

   integer, parameter :: tile_nc=3600, tile_nr=3600
   integer, parameter :: tile_nc_antarc=7200, tile_nr_antarc=3600

   integer, parameter :: tile_col=12, tile_row=6   ! Added 6th for Antarctica
   !**  60 tiles come as 5 rows, 12 in each for north of -60 S.
   integer, parameter :: tile_col_antarc=6
   !**  6 tiles for Antarctica 4th row, 6 cols
   integer :: tilecol_temp

   character*3 :: ns(tile_row)
   character*4 :: we(tile_col)
   character*4 :: we_antarc(tile_col_antarc)
   integer     :: col_cnt1, col_cnt2, row_cnt1, row_cnt2

   integer   :: ftn
   logical   :: file_exists
   integer   :: c, r, t, i, j, k, l, line
   integer   :: subpnc, subpnr, glpnc, glpnr
   integer   :: mi                          ! Total number of input param grid array points
   integer   :: mo                          ! Total number of output LIS grid array points
   real      :: param_gridDesc(20)          ! Input parameter grid desc array
   real      :: subparam_gridDesc(20)       ! Subsetted Input parameter grid desc array
   integer, allocatable :: lat_line(:,:), lon_line(:,:)
   integer, allocatable :: n11(:)           ! Maps each input grid point to output grid.
   real,    allocatable :: gi1(:)           ! input parameter 1d grid
   logical*1,allocatable:: li1(:)           ! input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go1)

   real, allocatable :: read_elevtile(:,:,:)               ! Tiled elevation fields
   integer*2, allocatable :: read_elevtile_antarc(:,:,:)   ! Tiled files - Antarctica

   real, allocatable :: mosaic_elev(:,:)    ! Mosaicked-tile elevation
   real, allocatable :: yrev_elev(:,:)      ! Y-reversed mosaicked-tile elevation
   real, allocatable :: subset_elev(:,:)    ! Read input parameter
   real, allocatable :: subset_slope(:,:)   ! Derived from input parameter

   character(len=LDT_CONST_PATH_LEN) :: tempfile
!________________________________________________________________________

  fgrd = 0.
  slopeave = 0.

!- Set parameter grid fgrd inputs:
   LDT_rc%topo_proj = 'latlon'
   param_gridDesc = 0.
   param_gridDesc(1)  = 0.             ! Latlon
   param_gridDesc(2)  = input_cols
   param_gridDesc(3)  = input_rows
   param_gridDesc(4)  = -90.0  + (input_yres/2) ! LL lat
   param_gridDesc(5)  = -180.0 + (input_xres/2) ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (input_yres/2)  ! UR lat
   param_gridDesc(8)  = 180.0 - (input_xres/2)  ! UR lon
   param_gridDesc(9)  = input_yres     ! dy: 0.0083333
   param_gridDesc(10) = input_xres     ! dx: 0.0083333
   param_gridDesc(20) = 64
  
! MERIT DEM naming convention:
   ns(1)="n90"
   ns(2)="n60"
   ns(3)="n30"
   ns(4)="n00"
   ns(5)="s30"
   we(1)="w180"
   we(2)="w150"
   we(3)="w120"
   we(4)="w090"
   we(5)="w060"
   we(6)="w030"
   we(7)="e000"
   we(8)="e030"
   we(9)="e060"
   we(10)="e090"
   we(11)="e120"
   we(12)="e150"

!- Upated GTOPO30 - Antarctica only files:
!  https://lta.cr.usgs.gov/GTOPO30
   ns(6)="s60"          ! Added to account for Antarctica
   we_antarc(1)="w180"
   we_antarc(2)="w120"
   we_antarc(3)="w060"
   we_antarc(4)="w000"
   we_antarc(5)="e060"
   we_antarc(6)="e120"
! __________________________

  if( LDT_rc%lis_map_proj(n) == "latlon"   .or. &
      LDT_rc%lis_map_proj(n) == "mercator" .or. &
      LDT_rc%lis_map_proj(n) == "lambert" ) then
     if( param_gridDesc(10) .ne. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) .and.&
         LDT_rc%topo_gridtransform(n) .eq. "none" ) then
        write(LDT_logunit,*) "[ERR] MERIT_1K has been selected which has a resolution"
        write(LDT_logunit,*) "    (0.00833deg), but the LIS run domain resolution"
        write(LDT_logunit,*) "    selected is not equal to that. So please select a spatial"
        write(LDT_logunit,*) "    transform other than 'none', if you want to read in"
        write(LDT_logunit,*) "    the resolution of MERIT_1K."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif
  endif

! Check if MERIT file directory exists:
  tempfile = trim(LDT_rc%elevfile(n))//"/merit_"//we(2)//ns(2)//".dem"
  inquire(file=tempfile, exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] MERIT elevation map file, ",&
                           trim(tempfile),", not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(LDT_logunit,*) "[INFO] Reading MERIT elevation files",&
                       " found in directory: ",trim(LDT_rc%elevfile(n))
! --
!  Open and read in tile files:
! --
   allocate( read_elevtile(tile_nc, tile_nr, tile_col) )
   allocate( read_elevtile_antarc(tile_nc_antarc, tile_nr_antarc, tile_col_antarc) )
   allocate( mosaic_elev(input_cols,input_rows) )
   read_elevtile = LDT_rc%udef
   read_elevtile_antarc = LDT_rc%udef
   mosaic_elev = LDT_rc%udef
   row_cnt1 = 0;  col_cnt1 = 0   ! tile_nc=3600, tile_nr=3600
   row_cnt2 = 0;  col_cnt2 = 0

 ! Read in tiled files:
   do l = 1, tile_row      ! = 6 now to account for Antarctica
      if( l < 6 ) then
        row_cnt1 = (l-1)*tile_nr+1
        row_cnt2 = row_cnt1+tile_nr-1
        tilecol_temp = tile_col
      elseif( l == 6 ) then   ! Antarctica
        row_cnt1 = (l-1)*tile_nr+1
        row_cnt2 = row_cnt1+tile_nr_antarc-1
        tilecol_temp = tile_col_antarc
      endif

      ! MERIT - Globe:
      do k = 1, tilecol_temp   ! read one row of tiles (total=12)
         ftn = LDT_getNextUnitNumber()
         if( l < 6 ) then
         open(ftn, file=trim(LDT_rc%elevfile(n))//"/merit_"//we(k)//ns(l)//".dem", &
              form="unformatted", access="direct", status="old", &
              convert='little_endian',recl=tile_nc*tile_nr*4)
         read(ftn, rec=1) read_elevtile(:, :, k)

      !- Mosaic all elevation tiles together:
         col_cnt1 = (k-1)*tile_nc+1
         col_cnt2 = col_cnt1+tile_nc-1
         mosaic_elev(col_cnt1:col_cnt2, row_cnt1:row_cnt2) = real(read_elevtile(:,:,k))

         ! GTOPO30 - Antarctica:
         elseif( l == 6 ) then
           tempfile = trim(LDT_rc%elevfile(n))//"/gt30"//we_antarc(2)//ns(6)//".dem"
           inquire(file=tempfile, exist=file_exists)
           if(.not.file_exists) then
             write(LDT_logunit,*) &
              "[WARN] MISSING GTOPO30 ANTARCTICA ELEVATION FILE: ",&
                trim(tempfile)
             write(LDT_logunit,*) "[WARN] GRIDCELLS OVER ANTARCTICA WILL BE FILLED."
           else
             open(ftn, file=trim(LDT_rc%elevfile(n))//"/gt30"//we_antarc(k)//ns(l)//".dem", &
                form="unformatted", access="direct", status="old", &
                convert='big_endian',recl=tile_nc_antarc*tile_nr_antarc*2)
             read(ftn, rec=1) read_elevtile_antarc(:, :, k)
           endif
           ! Mosaic all elevation tiles together:
           col_cnt1 = (k-1)*tile_nc_antarc+1
           col_cnt2 = col_cnt1+tile_nc_antarc-1
           mosaic_elev(col_cnt1:col_cnt2, row_cnt1:row_cnt2) = real(read_elevtile_antarc(:,:,k))
         endif
         call LDT_releaseUnitNumber(ftn)

      end do ! k
    end do   ! l
    deallocate( read_elevtile )
    deallocate( read_elevtile_antarc )

  ! Reverse-Y rows of read-in input file:
    allocate( yrev_elev(input_cols,input_rows) )
    yrev_elev = LDT_rc%udef
    i = 0
    do r = input_rows, 1, -1
       i = i + 1
       do c = 1, input_cols
          yrev_elev(c,i) = mosaic_elev(c,r)
       end do
    end do
    deallocate( mosaic_elev )

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%topo_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )


!- Subset parameter read-in array:
   allocate( subset_elev(subpnc, subpnr) )
   subset_elev = LDT_rc%udef
   do r = 1, subpnr
      do c = 1, subpnc
         subset_elev(c,r) = yrev_elev(lon_line(c,r),lat_line(c,r))
       
         ! for coastal areas
         if (subset_elev(c,r) .eq. LDT_rc%udef) then
            subset_elev(c,r) = 0.
         endif
      enddo
   enddo
   deallocate( yrev_elev )

!- Estimate slope field from elevation field:
   allocate( subset_slope(subpnc, subpnr) )
   subset_slope = LDT_rc%udef

   call calc_Slope_fromElev( subset_elev, subset_slope, input_xres, &
             subparam_gridDesc(5), subparam_gridDesc(4), subpnc, subpnr )

   deallocate( subset_elev )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   allocate( gi1(mi), li1(mi) )
   gi1 = LDT_rc%udef; li1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi1(i) = subset_slope(c,r)
         if( gi1(i) .ne. LDT_rc%udef ) li1(i) = .true.
      enddo
   enddo
   deallocate( subset_slope )

!- Transform parameter grid to LIS run domain:
   select case ( LDT_rc%topo_gridtransform(n) )

 !- Transforming 2-D elevation field: 
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )

   !- Transform parameter from original grid to LIS output grid:
      lo1 = .false.
      call LDT_transform_paramgrid(n, LDT_rc%topo_gridtransform(n), &
               subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n); i = i + 1
            if( go1(i) < -1. ) go1(i) = LDT_rc%udef
         !- Single elevation layer, write to first bin:
            fgrd(c,r,1) = 1.0
            slopeave(c,r,1) = go1(i)
         enddo
      enddo

 !- Transforming 3-D elevation field: 
    case( "tile" )

        allocate( n11(mi) )
     !- Create mapping between parameter domain and LIS grid domain:
        call upscaleByAveraging_input( subparam_gridDesc, LDT_rc%gridDesc(n,:), &
                                       mi, mo, n11 )

        call param_1dbin_areacalc( n, num_bins, mi, mo, n11, LDT_rc%udef, & 
                                   gi1, fgrd, slopeave )
        deallocate( n11 )

    case default
       write(LDT_logunit,*) "[ERR] The spatial transform, ",trim(LDT_rc%topo_gridtransform(n)),&
                            ", for MERIT elevation is not available at this time ... "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select
   deallocate( gi1, li1 )

   write(LDT_logunit, *) "[INFO] Done reading MERIT elevation files."

end subroutine read_MERIT1K_slope

