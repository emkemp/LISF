!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module jules54_lsmMod
!BOP
!
! !MODULE: jules54_lsmMod
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!  
! \begin{description}
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[numout]
!    number of output times 
!  \item[outInterval]
!    output writing interval
!  \item[jules54open]
!    variable to keep track of opened files
!  \item[jules54]
!   Template LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
! 28 Nov 2018; Shugong Wang; updated for JULES 5.2 
! 12 Dec 2018; Shugong Wang; updated for JULES 5.4

! !USES:        

  use jules54_module
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: jules54_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: jules54_struc
!EOP
  type, public :: jules54_type_dec
     integer                    :: jules54open
     integer                    :: numout
     integer                    :: nsurft     ! nsurf
     integer                    :: nsoilt
     integer                    :: sm_levels    ! sm_levels, default 4 
     integer                    :: nsmax     ! maximum number of snow layers, default 0, read from jules_snow.nml 
     integer                    :: temp_tiles 
     integer                    :: temp_layers
     integer                    :: ns_deep   ! number of bedrock layers, default 100
     integer                    :: dim_cs1   ! Soil Organic Nitrogen
     integer                    :: dim_cslayer   ! Soil Organic Nitrogen
     integer                    :: nice      ! Number of sea ice categories available
     integer                    :: nice_use  ! Number of sea ice categories in use
     integer                    :: npft_trif 
     integer                    :: npft, nnvg, ntype
     integer                    :: forc_count 
     real                       :: ts
     real                       :: z1_tq
     real                       :: z1_uv
     integer                    :: count
     real                       :: rstInterval
     integer                    :: outInterval
     character(len=LIS_CONST_PATH_LEN) :: namelist_dir
     character(len=LIS_CONST_PATH_LEN) :: rfile 
     character(len=32)          :: rformat
     type(jules54dec), allocatable :: jules54(:)
  end type jules54_type_dec
  type(jules54_type_dec), allocatable :: jules54_struc(:)
  
  SAVE
contains
!BOP
! 
! !ROUTINE: jules54_ini
! \label{jules54_ini}
! 
! !INTERFACE:
  subroutine jules54_ini()
! !USES:
   use ESMF
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_coreMod, only : LIS_rc
   use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
        LIS_update_timestep, LIS_registerAlarm
   use LIS_logMod,      only : LIS_logunit,LIS_endrun
   USE jules_surface_types_mod,  ONLY: npft, nnvg, ntype
   USE jules_soil_mod, ONLY : dzsoil, sm_levels, ns_deep
   use ancil_info, only : dim_cs1, nsurft, npft_trif, dim_cslayer, nsoilt, z1_tq_ij, z1_uv_ij
   use jules_snow_mod, only : nsmax, cansnowtile
   use jules_sea_seaice_mod, only : nice, nice_use
   use init_mod, only : init 
   use jules_surface_mod,      only: l_aggregate
   use model_time_mod, only: timestep_len
! !DESCRIPTION:        
!
!EOP
   implicit none
   integer :: n, t, k
   integer :: yr, mo, da, hr, mn, ss
   integer :: status
   integer :: counter  
   real, parameter :: tolerance=1e-6
   character*3 :: fnest ! EMK

   allocate(jules54_struc(LIS_rc%nnest))
   call jules54_readcrd()
   ! call jules initialization 
   do n=1, lis_rc%nnest
      call init(jules54_struc(n)%namelist_dir)
      ! In init subroutine, CALL init_time(nml_dir) 
      ! init_time reads in time steps from JULES timestep.nml
      ! double check here to make sure LIS and JULES have the same timestep settings
      ! Shugong Wang 10/03/2018 
!      timestep_len = int(LIS_rc%ts) 
!      timestep_real = real(timestep_len)
      if(timestep_len .ne. int(jules54_struc(n)%ts)) then 
         write(LIS_logunit,*)'Fatal JULES error:' 
         write(LIS_logunit,*)'  The timestep in the LIS configuration file is not consistent with'
         write(LIS_logunit,*)'  the timestep (timestep_len) in the JULES configuration namelist'
         write(LIS_logunit,*)'  file '//trim(jules54_struc(n)%namelist_dir)//'/timesteps.nml'
         write(LIS_logunit,*)'  LIS   timestep (seconds):', int(jules54_struc(n)%ts) 
         write(LIS_logunit,*)'  JULES timestep (seconds):', timestep_len 
         write(LIS_logunit,*)'  Stopping run....'
         call LIS_endrun 
      endif

      ! double check the reference heights for wind and temperature/humidity 
      if(z1_uv_ij(1,1) .ne. jules54_struc(n)%z1_uv) then 
         write(LIS_logunit,*)'Fatal JULES error:' 
         write(LIS_logunit,*)'  The "reference height for forcing u and v" in the LIS configuration file'
         write(LIS_logunit,*)'  is not consistent with "z1_uv_in" in the JULES configuration namelist'
         write(LIS_logunit,*)'  file '//trim(jules54_struc(n)%namelist_dir)//'/drive.nml'
         write(LIS_logunit,*)'  LIS   reference height for forcing u and v:', jules54_struc(n)%z1_uv
         write(LIS_logunit,*)'  JULES reference height for forcing u and v:', z1_uv_ij(1,1)
         write(LIS_logunit,*)'  Stopping run....'
         call LIS_endrun 
      endif 
      
      if(z1_tq_ij(1,1) .ne. jules54_struc(n)%z1_tq) then 
         write(LIS_logunit,*)'Fatal JULES error:' 
         write(LIS_logunit,*)'  The "reference height for forcing T and q" in the LIS configuration file'
         write(LIS_logunit,*)'  is not consistent with "z1_tq_in" in the JULES configuration namelist'
         write(LIS_logunit,*)'  file '//trim(jules54_struc(n)%namelist_dir)//'/drive.nml'
         write(LIS_logunit,*)'  LIS   reference height for forcing T and q:', jules54_struc(n)%z1_tq
         write(LIS_logunit,*)'  JULES reference height for forcing T and q:', z1_tq_ij(1,1)
         write(LIS_logunit,*)'  Stopping run....'
         call LIS_endrun 
      endif 
      
      
      !! from jules namelist
      jules54_struc(n)%sm_levels = sm_levels 
      jules54_struc(n)%ns_deep = ns_deep
      jules54_struc(n)%dim_cs1 = dim_cs1
      jules54_struc(n)%dim_cslayer = dim_cslayer
      
      jules54_struc(n)%nsurft = nsurft

      jules54_struc(n)%nsoilt = nsoilt 
      jules54_struc(n)%npft_trif = npft_trif 
      
      if (nsmax >0) then
        jules54_struc(n)%nsmax = nsmax 
        jules54_struc(n)%temp_tiles = nsurft
        jules54_struc(n)%temp_layers = nsmax
      else
        jules54_struc(n)%nsmax = 1
        jules54_struc(n)%temp_tiles = 1
        jules54_struc(n)%temp_layers = 1
      endif

      jules54_struc(n)%nice = nice 
      jules54_struc(n)%nice_use = nice_use 
      jules54_struc(n)%npft = npft
      jules54_struc(n)%nnvg = nnvg 
      jules54_struc(n)%ntype = ntype
      allocate(jules54_struc(n)%jules54(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      jules54_struc(n)%numout = 0

       
      ! allocate memory for state variables
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
      ! prognostics 
         allocate(jules54_struc(n)%jules54(t)%nsnow(jules54_struc(n)%nsurft))                               ! nsurf  number of snow layers on ground on tiles
         allocate(jules54_struc(n)%jules54(t)%tsoil_deep(jules54_struc(n)%ns_deep))                        ! ns_deep  deep soil temperatures (k)
         allocate(jules54_struc(n)%jules54(t)%sice(jules54_struc(n)%temp_tiles, jules54_struc(n)%temp_layers))      ! temp_tiles, temp_layers, Snow layer ice mass on tiles (Kg/m2)
         allocate(jules54_struc(n)%jules54(t)%sliq(jules54_struc(n)%temp_tiles, jules54_struc(n)%temp_layers))        ! temp_tile, temp_layers, nsmax, Snow layer liquid mass on tiles (Kg/m2)
         allocate(jules54_struc(n)%jules54(t)%snowdepth(jules54_struc(n)%nsurft))                           ! nsurf, Snow depth on ground on tiles (m)
         allocate(jules54_struc(n)%jules54(t)%tsnow(jules54_struc(n)%temp_tiles, jules54_struc(n)%temp_layers))      ! nsurf, nsmax,Snow layer temperature (K)
         allocate(jules54_struc(n)%jules54(t)%rgrainl(jules54_struc(n)%temp_tiles, jules54_struc(n)%temp_layers))    ! nsurf, nsmax, Snow layer grain size on tiles (microns)
         allocate(jules54_struc(n)%jules54(t)%rho_snow_grnd(jules54_struc(n)%nsurft))                       ! nsurf, Snowpack bulk density (kg/m3)
         allocate(jules54_struc(n)%jules54(t)%rho_snow(jules54_struc(n)%nsurft, jules54_struc(n)%nsmax))    ! nsurf, nsmax, Snow layer densities (m)
         allocate(jules54_struc(n)%jules54(t)%ds(jules54_struc(n)%temp_tiles, jules54_struc(n)%temp_layers))  ! temp_tiles, temp_layers, Snow layer thickness (m)
         allocate(jules54_struc(n)%jules54(t)%ns(jules54_struc(n)%dim_cslayer, jules54_struc(n)%dim_cs1))                                ! dim_cs1, Soil Organic Nitrogen (kg N m-2)
         allocate(jules54_struc(n)%jules54(t)%canht_ft(jules54_struc(n)%npft))                             ! npft, Canopy height (m)
         allocate(jules54_struc(n)%jules54(t)%canopy(jules54_struc(n)%nsurft))                              ! nsurf, Surface/canopy water for snow-free land tiles (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%cs_pool_soilt( jules54_struc(n)%dim_cslayer, jules54_struc(n)%dim_cs1))                                ! dim_cs1, Soil carbon (kg C/m2)
         allocate(jules54_struc(n)%jules54(t)%di_ncat(jules54_struc(n)%nice))                              ! nice, "Equivalent thickness" of sea-ice catagories (m)
         allocate(jules54_struc(n)%jules54(t)%k_sice(jules54_struc(n)%nice))                               ! nice, Sea ice effective conductivity (2*kappai/de)
         allocate(jules54_struc(n)%jules54(t)%gc_surft(jules54_struc(n)%nsurft))                                  ! nsurf, Stomatal" conductance to evaporation for land tiles(m s-1)
         allocate(jules54_struc(n)%jules54(t)%lai(jules54_struc(n)%npft))                                  ! npft LAI of plant functional types
         allocate(jules54_struc(n)%jules54(t)%rgrain(jules54_struc(n)%nsurft))                              ! nsurf, Snow surface grain size on tiles (microns)
         !allocate(jules54_struc(n)%jules54(t)%smc_soilt(jules54_struc(n)%nsoilt))                                ! sm_levels, Soil moisture content of layers (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%smcl_soilt(jules54_struc(n)%sm_levels))                                ! sm_levels, Soil moisture content of layers (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%snow_tile(jules54_struc(n)%nsurft))                           ! nsurf, Lying snow on tiles (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%snow_grnd(jules54_struc(n)%nsurft))                           ! nsurf, Snow on the ground (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%snow_mass_sea_sicat(jules54_struc(n)%nice_use))               !nice_use,  Snow on category sea-ice (Kg/m2)
         allocate(jules54_struc(n)%jules54(t)%t_soil(jules54_struc(n)%sm_levels))                              ! sm_levels, Sub-surface temperatures (K)
         allocate(jules54_struc(n)%jules54(t)%tstar_tile(jules54_struc(n)%nsurft))                          ! nsurf, Tile surface temperatures (K)
         jules54_struc(n)%jules54(t)%nsnow(:) = 0
         jules54_struc(n)%jules54(t)%tsoil_deep(:) = 0
         jules54_struc(n)%jules54(t)%sice(:,:) = 0
         jules54_struc(n)%jules54(t)%sliq(:,:) = 0
         jules54_struc(n)%jules54(t)%snowdepth(:) = 0
         jules54_struc(n)%jules54(t)%tsnow(:,:) = 0
         jules54_struc(n)%jules54(t)%rgrainl(:,:) = 0
         jules54_struc(n)%jules54(t)%rho_snow_grnd(:) = 0
         jules54_struc(n)%jules54(t)%rho_snow(:,:) = 0
         jules54_struc(n)%jules54(t)%ds(:,:) = 0
         jules54_struc(n)%jules54(t)%ns(:,:) = 0
         jules54_struc(n)%jules54(t)%canht_ft(:) = 0
         jules54_struc(n)%jules54(t)%canopy(:) = 0
         jules54_struc(n)%jules54(t)%cs_pool_soilt(:,:) = 0
         jules54_struc(n)%jules54(t)%di_ncat(:) = 0
         jules54_struc(n)%jules54(t)%k_sice(:) = 0
         jules54_struc(n)%jules54(t)%gc_surft(:) = 0
         jules54_struc(n)%jules54(t)%lai(:) = 0
         jules54_struc(n)%jules54(t)%rgrain(:) = 0
         !jules54_struc(n)%jules54(t)%smc_soilt(:) = 0
         jules54_struc(n)%jules54(t)%smcl_soilt(:) = 0
         jules54_struc(n)%jules54(t)%snow_tile(:) = 0
         jules54_struc(n)%jules54(t)%snow_grnd(:) = 0
         jules54_struc(n)%jules54(t)%snow_mass_sea_sicat = 0
         jules54_struc(n)%jules54(t)%t_soil(:) = 0
         jules54_struc(n)%jules54(t)%tstar_tile(:) = 0
         
         !p_s_parm
         allocate(jules54_struc(n)%jules54(t)%p_s_b( jules54_struc(n)%sm_levels ))                 !  Exponent for soil moisture characteristic function Clapp-Hornberger model: b is the Clapp-Hornberger exponent,  van Genuchten model: b=1/(n-1)  (metres)
         allocate(jules54_struc(n)%jules54(t)%p_s_sathh( jules54_struc(n)%sm_levels ))             !  Parameter for soil moisture characteristic functions Clapp-Hornberger model: sathh is the saturated soil water pressure (m), van Genuchten model: sathh=1/alpha
         allocate(jules54_struc(n)%jules54(t)%p_s_hcap( jules54_struc(n)%sm_levels ))              !  Soil heat capacity (J/K/m3)
         allocate(jules54_struc(n)%jules54(t)%p_s_hcon( 0:jules54_struc(n)%sm_levels ))              !  Soil thermal conductivity (W/m/K)
         allocate(jules54_struc(n)%jules54(t)%p_s_satcon( 0:jules54_struc(n)%sm_levels ))            !  Saturated hydraulic conductivity (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%p_s_smvccl( jules54_struc(n)%sm_levels ))            !  Critical volumetric SMC (cubic m per cubic m of soil)    
         allocate(jules54_struc(n)%jules54(t)%p_s_smvcst( jules54_struc(n)%sm_levels ))            !  Volumetric saturation point (m^3 m-3 of soil)
         allocate(jules54_struc(n)%jules54(t)%p_s_smvcwt( jules54_struc(n)%sm_levels ))            !  Volumetric wilting point (cubic m per cubic m of soil)
         allocate(jules54_struc(n)%jules54(t)%p_s_catch(jules54_struc(n)%nsurft ))             !  Surface/canopy water capacity of snow-free land tiles (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%p_s_catch_snow(jules54_struc(n)%nsurft ))        !  Snow interception capacity (kg m-2)
         allocate(jules54_struc(n)%jules54(t)%p_s_infil_tile(jules54_struc(n)%nsurft ))        !  Maximum possible surface infiltration for tiles (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%p_s_z0_tile(jules54_struc(n)%nsurft ))           !  Surface roughness on tiles (m).
         allocate(jules54_struc(n)%jules54(t)%p_s_z0h_tile_bare(jules54_struc(n)%nsurft ))     !  Surface thermal roughness on tiles before allowance for snow cover (m).
         allocate(jules54_struc(n)%jules54(t)%p_s_sthu( jules54_struc(n)%sm_levels ))              !  Unfrozen soil moisture content of the layers as a fraction of saturation.
         allocate(jules54_struc(n)%jules54(t)%p_s_sthf( jules54_struc(n)%sm_levels ))              !  Frozen soil moisture content of the layers as a fraction of saturation.
         allocate(jules54_struc(n)%jules54(t)%p_s_sthu_min( jules54_struc(n)%sm_levels ))          ! Minimum unfrozen water content for each layer. Used to normalise thaw depth calculation based on unfrozen water content fraction.
         allocate(jules54_struc(n)%jules54(t)%p_s_clay_soilt(jules54_struc(n)%dim_cslayer))                                ! dim_cs1, Soil Organic Nitrogen (kg N m-2)
         allocate(jules54_struc(n)%jules54(t)%p_s_v_close_pft(jules54_struc(n)%sm_levels, jules54_struc(n)%npft))
         allocate(jules54_struc(n)%jules54(t)%p_s_v_open_pft(jules54_struc(n)%sm_levels, jules54_struc(n)%npft))
         allocate(jules54_struc(n)%jules54(t)%p_s_soil_ph_soilt(jules54_struc(n)%sm_levels))
         allocate(jules54_struc(n)%jules54(t)%n_inorg_soilt_lyrs(jules54_struc(n)%dim_cslayer))                                !Gridbox Inorganic N pool on soil levels (kg N/m2) 
         allocate(jules54_struc(n)%jules54(t)%n_inorg_avail_pft(jules54_struc(n)%npft, jules54_struc(n)%dim_cslayer))
         allocate(jules54_struc(n)%jules54(t)%sice_pts_ncat(jules54_struc(n)%nice ))                                      ! number of points for each sea-ice category
         allocate(jules54_struc(n)%jules54(t)%sice_index_ncat(jules54_struc(n)%nice ))                                    ! index of points for each sea-ice category
         allocate(jules54_struc(n)%jules54(t)%sice_frac_ncat(jules54_struc(n)%nice ))                                     !  Fraction of gridbox covered by each sea-ice category (converted to single vector array)
         allocate(jules54_struc(n)%jules54(t)%tile_index(jules54_struc(n)%ntype ))                                        !  indices of land points which include the nth surface type
         allocate(jules54_struc(n)%jules54(t)%tile_pts(jules54_struc(n)%ntype ))                                          !  Number of land points which include the nth surface type
         allocate(jules54_struc(n)%jules54(t)%frac(jules54_struc(n)%ntype ))                                              !  fractional cover of each surface type
         allocate(jules54_struc(n)%jules54(t)%ice_fract_ncat(jules54_struc(n)%nice ))                                     !  fraction of gridbox covered by sea-ice on catagories
         allocate(jules54_struc(n)%jules54(t)%ti_cat(jules54_struc(n)%nice ))                                             ! sea ice surface temperature on categories
         allocate(jules54_struc(n)%jules54(t)%ti_sicat(jules54_struc(n)%nice ))                                             ! sea ice surface temperature on categories
         allocate(jules54_struc(n)%jules54(t)%tsurf_elev_surft(jules54_struc(n)%nsurft ))                                             ! sea ice surface temperature on categories
         allocate(jules54_struc(n)%jules54(t)%pond_frac_cat(jules54_struc(n)%nice ))                                      ! Meltpond fraction on sea ice categories 
         allocate(jules54_struc(n)%jules54(t)%pond_depth_cat(jules54_struc(n)%nice ))                                     ! Meltpond depth on sea ice categories (m)  
         allocate(jules54_struc(n)%jules54(t)%l_lice_surft(jules54_struc(n)%ntype ))                                     ! Meltpond depth on sea ice categories (m)  
         jules54_struc(n)%jules54(t)%p_s_b(:) = 0
         jules54_struc(n)%jules54(t)%p_s_sathh(:) = 0
         jules54_struc(n)%jules54(t)%p_s_hcap(:) = 0
         jules54_struc(n)%jules54(t)%p_s_hcon(:) = 0
         jules54_struc(n)%jules54(t)%p_s_satcon(:) = 0 
         jules54_struc(n)%jules54(t)%p_s_smvccl(:) = 0
         jules54_struc(n)%jules54(t)%p_s_smvcst(:) = 0
         jules54_struc(n)%jules54(t)%p_s_smvcwt(:) = 0
         jules54_struc(n)%jules54(t)%p_s_catch(:) = 0
         jules54_struc(n)%jules54(t)%p_s_catch_snow(:) = 0
         jules54_struc(n)%jules54(t)%p_s_infil_tile(:) = 0
         jules54_struc(n)%jules54(t)%p_s_z0_tile(:) = 0
         jules54_struc(n)%jules54(t)%p_s_z0h_tile_bare(:) = 0
         jules54_struc(n)%jules54(t)%p_s_sthu(:) = 0
         jules54_struc(n)%jules54(t)%p_s_sthf(:) = 0
         jules54_struc(n)%jules54(t)%p_s_sthu_min(:) = 0
         jules54_struc(n)%jules54(t)%sice_pts_ncat(:) = 0
         jules54_struc(n)%jules54(t)%sice_index_ncat(:) = 0
         jules54_struc(n)%jules54(t)%sice_frac_ncat(:) = 0
         jules54_struc(n)%jules54(t)%tile_index(:) = 0
         jules54_struc(n)%jules54(t)%tile_pts(:) = 0
         jules54_struc(n)%jules54(t)%frac(:) = 0
         jules54_struc(n)%jules54(t)%ice_fract_ncat(:) = 0
         jules54_struc(n)%jules54(t)%ti_cat(:) = 0
         jules54_struc(n)%jules54(t)%ti_sicat(:) = 0
         jules54_struc(n)%jules54(t)%pond_frac_cat(:) = 0
         jules54_struc(n)%jules54(t)%pond_depth_cat(:) = 0

         allocate(jules54_struc(n)%jules54(t)%isoprene_ft(jules54_struc(n)%npft))          ! Isoprene emission flux on PFTs (kgC/m2/s) 
         allocate(jules54_struc(n)%jules54(t)%terpene_ft(jules54_struc(n)%npft))           ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s) 
         allocate(jules54_struc(n)%jules54(t)%methanol_ft(jules54_struc(n)%npft))          ! Methanol emission flux on PFTs (kgC/m2/s) 
         allocate(jules54_struc(n)%jules54(t)%acetone_ft(jules54_struc(n)%npft))           ! Acetone emission flux on PFTs (kgC/m2/s) 
         allocate(jules54_struc(n)%jules54(t)%g_leaf_acc(jules54_struc(n)%npft))           ! Accumulated leaf turnover rate
         allocate(jules54_struc(n)%jules54(t)%npp_ft_acc(jules54_struc(n)%npft_trif))      ! Accumulated NPP_FT
         allocate(jules54_struc(n)%jules54(t)%g_leaf_phen_acc(jules54_struc(n)%npft))      ! Accumulated leaf turnover rate including phenology
         allocate(jules54_struc(n)%jules54(t)%resp_w_ft_acc(jules54_struc(n)%npft_trif))   ! Accum RESP_W_FT
         allocate(jules54_struc(n)%jules54(t)%resp_s_acc_soilt(jules54_struc(n)%dim_cslayer,jules54_struc(n)%dim_cs1))        ! Accumulated RESP_S
         allocate(jules54_struc(n)%jules54(t)%g_leaf(jules54_struc(n)%npft))               ! Leaf turnover rate (/360days)
         allocate(jules54_struc(n)%jules54(t)%g_leaf_phen(jules54_struc(n)%npft))          ! Mean leaf turnover rate over phenology period(/360days)
         allocate(jules54_struc(n)%jules54(t)%gpp_ft(jules54_struc(n)%npft))               ! Gross primary productivity on PFTs (kg C/m2/s)
         allocate(jules54_struc(n)%jules54(t)%npp_ft(jules54_struc(n)%npft))               ! Net primary productivity on PFTs (kg C/m2/s)
         allocate(jules54_struc(n)%jules54(t)%resp_p_ft(jules54_struc(n)%npft))            ! Plant respiration on PFTs (kg C/m2/s)
         allocate(jules54_struc(n)%jules54(t)%resp_s_soilt(jules54_struc(n)%dim_cslayer, jules54_struc(n)%dim_cs1))            ! Soil respiration (kg C/m2/s)
         allocate(jules54_struc(n)%jules54(t)%resp_w_ft(jules54_struc(n)%npft))            ! Wood maintenance respiration (kg C/m2/s)
         allocate(jules54_struc(n)%jules54(t)%lai_phen(jules54_struc(n)%npft))             ! LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A
         allocate(jules54_struc(n)%jules54(t)%c_veg(jules54_struc(n)%npft))                ! Total carbon content of the vegetation (kg C/m2)
         allocate(jules54_struc(n)%jules54(t)%g_leaf_day(jules54_struc(n)%npft))           ! Mean leaf turnover rate for input to PHENOL (/360days)
         allocate(jules54_struc(n)%jules54(t)%g_leaf_dr_out(jules54_struc(n)%npft))        ! Mean leaf turnover rate for driving TRIFFID (/360days)
         allocate(jules54_struc(n)%jules54(t)%lit_c(jules54_struc(n)%npft))                ! Carbon Litter (kg C/m2/360days)
         allocate(jules54_struc(n)%jules54(t)%npp_dr_out(jules54_struc(n)%npft))           ! Mean NPP for driving TRIFFID (kg C/m2/360days)
         allocate(jules54_struc(n)%jules54(t)%resp_w_dr_out(jules54_struc(n)%npft))        ! Mean wood respiration for driving TRIFFID (kg C/m2/360days)
         allocate(jules54_struc(n)%jules54(t)%resp_s_dr_out_gb(jules54_struc(n)%dim_cslayer, 5))                            ! Mean soil respiration for driving TRIFFID (kg C/m2/360days)
         jules54_struc(n)%jules54(t)%isoprene_ft(:) = 0
         jules54_struc(n)%jules54(t)%terpene_ft(:) = 0
         jules54_struc(n)%jules54(t)%methanol_ft(:) = 0
         jules54_struc(n)%jules54(t)%acetone_ft(:) = 0
         jules54_struc(n)%jules54(t)%g_leaf_acc(:) = 0
         jules54_struc(n)%jules54(t)%npp_ft_acc(:) = 0
         jules54_struc(n)%jules54(t)%g_leaf_phen_acc(:) = 0
         jules54_struc(n)%jules54(t)%resp_w_ft_acc(:) = 0
         jules54_struc(n)%jules54(t)%resp_s_acc_soilt(:,:) = 0
         jules54_struc(n)%jules54(t)%g_leaf(:) = 0
         jules54_struc(n)%jules54(t)%g_leaf_phen(:) = 0
         jules54_struc(n)%jules54(t)%gpp_ft(:) = 0
         jules54_struc(n)%jules54(t)%npp_ft(:) = 0
         jules54_struc(n)%jules54(t)%resp_p_ft(:) = 0
         jules54_struc(n)%jules54(t)%resp_s_soilt(:,:) = 0
         jules54_struc(n)%jules54(t)%resp_w_ft(:) = 0
         jules54_struc(n)%jules54(t)%lai_phen(:) = 0
         jules54_struc(n)%jules54(t)%c_veg(:) = 0
         jules54_struc(n)%jules54(t)%g_leaf_day(:) = 0
         jules54_struc(n)%jules54(t)%g_leaf_dr_out(:) = 0
         jules54_struc(n)%jules54(t)%lit_c(:) = 0
         jules54_struc(n)%jules54(t)%npp_dr_out(:) = 0
         jules54_struc(n)%jules54(t)%resp_w_dr_out(:) = 0
         jules54_struc(n)%jules54(t)%resp_s_dr_out_gb(:,:) = 0

         if (any(cansnowtile(1:npft)) .eqv. .true.) then
            allocate( jules54_struc(n)%jules54(t)%unload_backgrnd(npft))
         else
            allocate( jules54_struc(n)%jules54(t)%unload_backgrnd(1) )
         endif

         ! fluxes 

         allocate(jules54_struc(n)%jules54(t)%anthrop_heat_surft(jules54_struc(n)%nsurft))! 
         allocate(jules54_struc(n)%jules54(t)%surf_ht_store_surft(jules54_struc(n)%nsurft))! 
         allocate(jules54_struc(n)%jules54(t)%snow_soil_htf(jules54_struc(n)%nsurft))! 
         allocate(jules54_struc(n)%jules54(t)%alb_tile(jules54_struc(n)%nsurft, 4)      )!  Albedo for surface tiles, 1: direct beam visible, 2: diffuse visible, 3: direct beam near-IR, 4: diffuse near-IR   
         allocate(jules54_struc(n)%jules54(t)%fsmc(jules54_struc(n)%npft )             )!  Moisture availability factor
         allocate(jules54_struc(n)%jules54(t)%ftl_tile(jules54_struc(n)%nsurft)         )!  Surface FTL for land tiles (W m-2)
         allocate(jules54_struc(n)%jules54(t)%le_tile(jules54_struc(n)%nsurft)          )!  Surface latent heat flux for land tiles (W m-2) 
         allocate(jules54_struc(n)%jules54(t)%fqw_tile(jules54_struc(n)%nsurft)         )!  Surface FQW (Moisture flux between layers) for land tiles, (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%fqw_ice(jules54_struc(n)%nice_use)       )!  Surface FQW (Moisture flux between layers) for sea-ice, (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%ftl_ice(jules54_struc(n)%nice_use)       )!  Surface FTL for sea-ice (W m-2) 
         allocate(jules54_struc(n)%jules54(t)%esoil_tile(jules54_struc(n)%nsurft)       )!  ESOIL for snow-free land tiles
         allocate(jules54_struc(n)%jules54(t)%sea_ice_htf(jules54_struc(n)%nice)       )!  Heat flux through sea-ice (W m-2, positive downwards)
         allocate(jules54_struc(n)%jules54(t)%surf_htf_tile(jules54_struc(n)%nsurft)    )!  Surface heat flux on land tiles (W m-2)
         allocate(jules54_struc(n)%jules54(t)%land_albedo(4)                           )!  GBM albedo, 1-direct beam visible, 2-diffuse visible, 3-direct beam near-IR, 4-diffuse near-IR
         allocate(jules54_struc(n)%jules54(t)%ei_tile(jules54_struc(n)%nsurft)          )!  EI for land tiles (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%ecan_tile(jules54_struc(n)%nsurft )       )!  Canopy evaporation from for snow-free land tiles (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%ext(jules54_struc(n)%sm_levels)              )!  Extraction of water from each soil layer (kg m-2/s)
         allocate(jules54_struc(n)%jules54(t)%radnet_tile(jules54_struc(n)%nsurft)      )!  Surface net radiation on tiles ( W m-2)
         allocate(jules54_struc(n)%jules54(t)%sw_tile(jules54_struc(n)%nsurft)          )!  Surface net shortwave on tiles (W m-2)
         allocate(jules54_struc(n)%jules54(t)%emis_tile(jules54_struc(n)%nsurft)        )!  Tile emissivity
         allocate(jules54_struc(n)%jules54(t)%melt_tile(jules54_struc(n)%nsurft)        )!  Snowmelt on land tiles (kg m-2/s)

         if(l_aggregate) then
            if(LIS_rc%surface_maxt /= 1) then 
              write(LIS_logunit,*)'Fatal JULES error:' 
              write(LIS_logunit,*)'  Maximum number of surface type tiles per grid should be 1'
              write(LIS_logunit,*)'  when l_aggreate is set to .true. in jules_surface.nml'
              call LIS_endrun 
            else
              allocate(jules54_struc(n)%jules54(t)%surft_frac(jules54_struc(n)%ntype)) ! surface type fractions 
            endif 
         endif
      enddo
      
      ! initialize forcing variables to zeros
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         jules54_struc(n)%jules54(t)%lwdown = 0.0
         jules54_struc(n)%jules54(t)%swdown = 0.0
         jules54_struc(n)%jules54(t)%psurf = 0.0
         jules54_struc(n)%jules54(t)%rainf = 0.0
         jules54_struc(n)%jules54(t)%snowf = 0.0
         jules54_struc(n)%jules54(t)%tair = 0.0
         jules54_struc(n)%jules54(t)%qair = 0.0
         jules54_struc(n)%jules54(t)%wind_e = 0.0
         jules54_struc(n)%jules54(t)%wind_n = 0.0
      enddo ! end of tile (t) loop
      !------------------------------------------------------------------------
      ! Model timestep Alarm
      !------------------------------------------------------------------------
      jules54_struc(n)%forc_count = 0
      call LIS_update_timestep(LIS_rc, n, jules54_struc(n)%ts)
      LIS_sfmodel_struc(n)%ts = jules54_struc(n)%ts
      call LIS_update_timestep(LIS_rc, n, jules54_struc(n)%ts)
      call LIS_registerAlarm("JULES.5.4 model alarm",&
                             jules54_struc(n)%ts, &
                             jules54_struc(n)%ts)
      call LIS_registerAlarm("JULES.5.4 restart alarm", &
                             jules54_struc(n)%ts,&
                             jules54_struc(n)%rstInterval)

      ! EMK Add alarm to reset tair_agl_min for RHMin.  This should 
      ! match the output interval, since that is used for calculating 
      ! Tair_F_min.
      write(fnest,'(i3.3)') n
      call LIS_registerAlarm("JULES.5.4 RHMin alarm "//trim(fnest),&
           jules54_struc(n)%ts,&
           LIS_sfmodel_struc(n)%outInterval)
      if (LIS_sfmodel_struc(n)%outInterval .gt. 86400 .or. &
           trim(LIS_sfmodel_struc(n)%outIntervalType) .eq. "dekad") then
         write(LIS_logunit,*) &
              '[WARN] If RHMin is selected for output, please reset ', &
              'surface model output interval to no more than 24 hours.'
      end if
      
      !
      !
      !
      LIS_sfmodel_struc(n)%nsm_layers = jules54_struc(n)%sm_levels
      LIS_sfmodel_struc(n)%nst_layers = jules54_struc(n)%sm_levels
      allocate(LIS_sfmodel_struc(n)%lyrthk(jules54_struc(n)%sm_levels))
      ! EMK...Output soil layer thicknesses in centimeters.
!      LIS_sfmodel_struc(n)%lyrthk(1:jules54_struc(n)%sm_levels) = dzsoil(1:jules54_struc(n)%sm_levels) 
      LIS_sfmodel_struc(n)%lyrthk(1:jules54_struc(n)%sm_levels) = &
           100*dzsoil(1:jules54_struc(n)%sm_levels) 

      ! EMK...Additional variables for Air Force
      jules54_struc(n)%jules54(:)%tair_agl_min = 999.0
      jules54_struc(n)%jules54(:)%rhmin = 999.0

    enddo


    call jules54_search_end_k 

  end subroutine jules54_ini

  subroutine jules54_search_end_k()
     ! !USES:
     use LIS_coreMod,            only: LIS_rc, LIS_domain, LIS_surface
  !
  ! !DESCRIPTION:
  !
  !  This routine search the start and end tile indices for a grid 
  !
  !EOP

     implicit none
     integer :: t, n, j, k
     integer :: c, r
     real    :: lat, lon
     integer :: pft
     integer :: gid, cur_grid, start_k, end_k

     do n=1, LIS_rc%nnest
       k = 1
       do
          if ( k > LIS_rc%npatch(n,LIS_rc%lsm_index) ) then
             exit
          endif

          ! Find all tiles in current grid
          cur_grid = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%index
          start_k = k
          do
             k = k + 1
             if ( k > LIS_rc%npatch(n,LIS_rc%lsm_index) ) then
                gid = -9999 ! invalid grid number
             else
                gid = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%index
             endif
             if ( gid /= cur_grid ) then
                end_k = k - 1
                jules54_struc(n)%jules54(start_k)%end_k = end_k
                exit
             endif
          enddo
          
        enddo
     enddo
  end subroutine jules54_search_end_k

end module jules54_lsmMod

