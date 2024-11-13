!                                                                       
!  =================================================================    
!  DECIPHeR-GW VERSION 1
!
!  Bristol University 2018 (Gemma Coxon, Jim Freer, Rosie Lane and Toby Dunne)
!  Based on fortran 77 version of dynamic TOPMODEL produced in
!  Lancaster University 12/01/00 (Keith Beven & Jim Freer)              
!  Migrated to std=F2003 Toby Dunne 25/01/2015
!  Modified extensively by Gemma Coxon 2016-2018
!  Integrated with multiscale parameter regionalization - Rosie Lane 2017-2019
!  add with GW model - Yanchen Zheng 2023-2024
!  =================================================================
!                                                                       
program main

    use dyna_common_types
    use dyna_project
    use dyna_random
    use dyna_tread_dyna
    use dyna_inputs
    use dyna_mc_setup
    use dyna_modelstruct_setup
    use dyna_genpar
    use dyna_main_loop
    use dyna_file_open
    use dyna_read_routingdata

    use dyna_utility
    use dyna_route_processing
    use dyna_route_riv_tree_node

    use mpr_main
    use mpr_extract_BasPred
    use mpr_control
    use mpr_errors
    use mpr_upscaling
    use mpr_random

    !MPR settings and functions for each DynaTOP parameter
    use mpr_SZM
    use mpr_LnTo
    use mpr_SRmaxGW
    use mpr_SRinit
    use mpr_CHV
    use mpr_Td
    use mpr_Smax

    !YZ GW modules
    use GW_preprocess
    use GW_read
    use GW_modelNEW
    use mpr_B_KsGW
    use GWpr_T_Sy

    implicit none

    !  Local variable declares
    integer :: seed_1, seed_2
    integer :: i,ii,jj
    integer :: nac
    integer :: nstep
    integer :: num_rivers

    ! allocated from 'Input'
    double precision, dimension(:,:), allocatable :: pe_step
    double precision, dimension(:), allocatable :: qobs_riv_step_start
    double precision, dimension(:,:), allocatable :: r_gau_step

    ! allocated after 'Input'
    type(dyna_hru_type), dimension(:), allocatable :: dyna_hru
    type(dyna_riv_type), dimension(:), allocatable :: dyna_riv

    ! MULTI POINT RIVER ROUTING
    ! Toby Dunne April 2016 + GC June 2016
    ! these type are defined in modules as part of the dta files
    type(route_river_info_type) :: route_riv
    type(route_time_delay_hist_type) :: route_tdh
    integer :: cat_route_vmode
    integer :: routing_mode
    integer, dimension(:), allocatable :: node_to_flow_mapping
    character(900) :: out_dir_path

    ! Monte Carlo Simulations
    integer :: numsim, i_mc, num_par_types
    integer, allocatable, dimension(:) :: num_mcpar
    character(len=20), dimension(7) :: all_pm_names !names of all possible parameters - must be added to if new params are developed!
    doubleprecision, allocatable, dimension(:,:) :: mcpar_ll
    doubleprecision, allocatable, dimension(:,:) :: mcpar_ul
    doubleprecision, allocatable, dimension(:,:) :: mcpar

    ! Kinematic solution parameters and time step
    doubleprecision :: dt
    doubleprecision :: acc
    doubleprecision :: wt
    integer :: ntt
    doubleprecision :: dtt

    double precision, dimension(:,:), allocatable :: rivers !(nac, n_riv)
    double precision, dimension(:), allocatable :: sum_ac_riv  !

    character(1024):: arg
    character(1024):: auto_start_file

    !! ADDED VARIABLES FOR MPR

        !key filepaths and folder locations
    character(len=1024)     :: MPRfolder        !location of MPR folder containing INPUT and OUTPUT folders, and MPR_control.dat file
    character(len=1024)     :: fpath_HRUmap     !filepath pointing to HRU .asc file created in DTA
    character(len=1024)     :: folder_input   !folder containing Basin Predictor files - INPUT folder within MPR main folder
    character(len=1024)     :: folder_output   !output folder - OUTPUT folder within MPR main folder
    character(len=900)      :: folder_outfull          !full path to output folder
    character(len=1024)     :: fname_control    !full path to control file
    character(len=1024)     :: fname_filemgr    !full path to mpr filemanger file
    character(len=1024)     :: fpath_gp_file  !full path to the global parameter output text file.
    character(len=1024), dimension(27)   :: fnames_baspred !filenames of all basin predictors - read from filepath manager.
    
    !YZ GW new declare
    character(len=1024)     :: GW_alldata_fn,catch_fn,GWmaskfile
    character(len=6)     :: cats_name,tem_name
    double precision, allocatable, dimension(:,:) :: dta_data,GW_data
    double precision, allocatable, dimension(:,:) :: GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff,GEO_index
    double precision, allocatable, dimension(:,:) :: GW_Sy_clip,GW_T_clip,GW_Sytemp,GW_Ttemp
    double precision, allocatable, dimension(:,:) :: GW_hi_ini,GW_hNew_ini,GW_runoff_ini
    integer :: ncols, nrows, GW_ncols, GW_nrows,NR,NC
    double precision :: xllcorner, yllcorner, cellsize, nodata
    double precision :: GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata
    double precision, allocatable, dimension(:,:) ::GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w
    double precision, allocatable, dimension(:,:) ::GWgauge_meta
    integer :: ia
    integer :: j
    character(900) :: filename
    character(len=1024) :: line_format



    !increments for do statements etc.
    !    integer                 :: i_HRU
    !    integer                 :: i_n
    !    double precision        :: zero=0

    !variables in the control file
    integer                 ::  n_pm_maps       !number of parameter files to create
    integer                 ::  start_seed      !starting seed for generation of global params
    integer                 ::  save_pm_maps    !set to 1 if user wants to save parameter maps
    integer                 ::  save_bp_maps    !set to 1 if user wants to save basin predictor maps
    integer                 ::  save_GW_maps    !set to 1 if user wants to save GW maps
    integer                 ::  set_GW_buffer
    integer                 ::  GW_maxiteration
    integer                 ::  GW_tolerance
    double precision        ::  tolerance
    !    integer                 ::  pedo_tf_SZM
    !    integer                 ::  pedo_tf_LnTo
    !    integer                 ::  pedo_tf_SRmax   !pedo-transfer equation numbers for all params
    !    integer                 ::  pedo_tf_SRinit
    !    integer                 ::  pedo_tf_CHV
    !    integer                 ::  pedo_tf_Td
    !    integer                 ::  pedo_tf_Smax

    !parameter min, default and max values from parameter file
    !    double precision, dimension(3)  :: SZM_range
    !    double precision, dimension(3)  :: LnTo_range
    !    double precision, dimension(3)  :: SRmax_range
    !    double precision, dimension(3)  :: SRinit_range
    !    double precision, dimension(3)  :: CHV_range
    !    double precision, dimension(3)  :: Td_range
    !    double precision, dimension(3)  :: Smax_range
    double precision, dimension(9,3):: Params_range

    !variables combining results from the control file
    integer, dimension(11)   ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax

    !variables defining HRU map and parameter maps
    double precision, allocatable, dimension(:,:) :: HRU_map
    double precision, allocatable, dimension(:,:) :: pm_map_Drz
    !    integer                         :: ncols_HRU
    !    integer                         :: nrows_HRU        !all variables relating to size and location
    !    double precision                :: xllcorner_HRU    !of HRU map created during DTA
    !    double precision                :: yllcorner_HRU
    !    double precision                :: cellsize_HRU
    !    double precision                :: nodata_HRU
    !    integer                         :: ncol
    !    integer                         :: nrow             !all variables relating to size and location
    !    double precision                :: xll              !of param maps - assumed same as HRU map
    !    double precision                :: yll
    !    double precision                :: cellsize
    integer                         :: nHRUS            !number of HRUs
    double precision, dimension(4)  :: stats_HRU        !variable comprising HRU map xll, yll, cellsize & nodata value

    !basin predictor maps needed - logical arrays with true values if maps required
    logical, dimension(20) :: req_bp    !Sand, silt, clay and OM rows, split by depth columns

    !basin predictor maps
    real, allocatable, dimension(:,:,:)    :: bp_maps
    real, allocatable, dimension(:,:)   :: sand_0_10
    real, allocatable, dimension(:,:)   :: silt_0_10
    real, allocatable, dimension(:,:)   :: clay_0_10
    real, allocatable, dimension(:,:)   :: orgm_0_10
    double precision, allocatable, dimension(:,:,:) :: soilmusiddata    !SZM:tables containing depth decline data for each soil series and musid.
    real, allocatable, dimension(:,:)   :: bp_root_depths   !SRmax: landcover class rooting depth table
    double precision, allocatable, dimension(:,:)   :: B_Ks_table !YZ 2023
    real, allocatable, dimension(:,:)   :: soil_class !YZ 2023
    double precision, allocatable, dimension(:,:)   :: T_Sy_table !YZ 2024

    !global parameters
    integer, dimension(11)                           :: n_gp_all     !vector of number of global params for all params
    double precision, allocatable, dimension(:,:)   :: glob_pms_SZM     !dim(glob pm number, 1:n_pm_maps)
    double precision, allocatable, dimension(:,:)   :: glob_pms_LnTo
    double precision, allocatable, dimension(:,:)   :: glob_pms_SRmax
    double precision, allocatable, dimension(:,:)   :: glob_pms_SRinit
    double precision, allocatable, dimension(:,:)   :: glob_pms_CHV
    double precision, allocatable, dimension(:,:)   :: glob_pms_Td
    double precision, allocatable, dimension(:,:)   :: glob_pms_Smax

    !YZ 2023
    double precision, allocatable, dimension(:,:)   :: glob_pms_B
    double precision, allocatable, dimension(:,:)   :: glob_pms_Ks
    double precision, allocatable, dimension(:,:)   :: glob_pms_T
    double precision, allocatable, dimension(:,:)   :: glob_pms_Sy
    integer:: GW_print_control

    !parameter maps at the finest scale, created using transfer functions
    !    double precision, allocatable, dimension(:,:)   :: pm_map_SZM
    !    double precision, allocatable, dimension(:,:)   :: pm_map_LnTo
    !    double precision, allocatable, dimension(:,:)   :: pm_map_SRmax
    !    double precision, allocatable, dimension(:,:)   :: pm_map_SRinit
    !    double precision, allocatable, dimension(:,:)   :: pm_map_CHV
    !    double precision, allocatable, dimension(:,:)   :: pm_map_Td
    !    double precision, allocatable, dimension(:,:)   :: pm_map_Smax
    double precision, allocatable, dimension(:,:,:)   :: pm_map_all !mp_map_all(x,y,[SZM, LnTo, SRmax, SRinit, CHV, Td, Smax])

    !variables required for parameter upscaling
    !    double precision, allocatable, dimension(:)     :: pms_SZM       !vector of all param values in a HRU
    !    double precision, allocatable, dimension(:)     :: pms_LnTo
    !    double precision, allocatable, dimension(:)     :: pms_SRmax
    !    double precision, allocatable, dimension(:)     :: pms_SRinit
    !    double precision, allocatable, dimension(:)     :: pms_CHV
    !    double precision, allocatable, dimension(:)     :: pms_Td
    !    double precision, allocatable, dimension(:)     :: pms_Smax
    double precision, allocatable, dimension(:,:)   :: up_pms       !vector of upscaled params
    double precision, allocatable, dimension(:,:)   :: up_pms_Drz 

    !variables for saving Dynamic TOPMODEL parameter file
    !    character(len=1024)                             :: temp_fn      !stores dynamic filename
    !    character(len = 5)                              :: i_n_str      !string storing current number of pm files made
    !    character(len=64), dimension(8)                 :: headers      !parameter file headers
    !    character(len = 6)                              :: fmt = '(I5.5)' !format descriptor - integer of width 5 with zeros on left
    !    integer  :: x,y
    !    integer, dimension(2) :: test_lhs

    !variables for doing repeat runs altering some parameter maps
    integer, dimension(7)   ::  param_repeatruns !number of parameter map modifications needed for each parameter.
    double precision, dimension(7,4)    :: parammap_pert_ranges !user-defined ranges for additions and multiplications to parameter map.
    double precision, allocatable,dimension(:,:,:) :: parammap_add_mult !param map perturbations (nparams,nperturbations,add/mult)
    ! initialise pf_data
    type(pf_data_type) :: pf_data

    pf_data%pf_enabled = .false.


    auto_start_file = ''
    i = 0

    print *, '--- Starting DECIPHeR rainfall-runoff modelling ---'

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-auto')) then
            CALL get_command_argument(i+1, auto_start_file)
        endif
        i = i + 1
    end do

    !  Call project to read in project files
    call project (auto_start_file, &
        out_dir_path, &
        fname_control, &
        fname_filemgr)

    write(999,*) ''
    write(999,*) 'Reading in HRU files'

    ! Call tread_dyna to read in the HRU file
    call Tread_dyna (nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv, &
        rivers, &
        sum_ac_riv)

    write(999,*) ''
    write(999,*) 'Reading in Input files'

    ! Call inputs to read in the input data
    call Inputs (nstep, &
        num_rivers, &
        dt, &
        pe_step, &
        qobs_riv_step_start, &
        r_gau_step)

    write(999,*) ''
    write(999,*) 'Reading in Parameter files'

        ! Read in parameter data
    ! call mc_setup (ACC, &
    !     dyna_hru, &
    !     mcpar_ll, &
    !     mcpar_ul, &
    !     NTT, &
    !     num_mcpar, &
    !     num_par_types, &
    !     numsim, &
    !     seed_1, &
    !     seed_2, &
    !     WT, &
    !     all_pm_names)

    !  Initialise randowm number generator
    !
    ! call ranin (seed_1, seed_2)

    write(999,*) ''
    write(999,*) 'Reading in Model Structure files'

    ! Read in model structure info
    call modelstruct_setup (dyna_hru, &
        nac)

    write(999,*) ''
    write(999,*) 'Reading in Routing Data'

    ! call read_routingdata to read in the routing data
    call read_routingdata (cat_route_vmode, &
        dyna_hru, &
        nac, &
        node_to_flow_mapping, &
        num_rivers, &
        route_riv, &
        route_tdh, &
        routing_mode)

    write(999,*) ''
    write(999,*) 'Reading in MPR files'

    !call mpr_read_controls to read key variables from the MPR control file
    call mpr_read_controls(MPRfolder, &    !key filepaths and folder locations
        fpath_HRUmap, &
        folder_input, &
        folder_output, &
        fname_control, &
        fname_filemgr, &
        fnames_baspred, &
        Params_range, &      !parameter min, default and max values from parameter file
        n_pm_maps, &         !variables in the control file
        start_seed, &
        save_pm_maps, &
        save_bp_maps, &
        pedo_tf_all, &
        param_repeatruns,&
        parammap_pert_ranges, &
        fpath_gp_file, &
        save_GW_maps, &
        set_GW_buffer,GW_maxiteration,GW_tolerance)
        if (GW_tolerance==0) then
           !tolerance in GW model (default)
           tolerance=1e-6
        endif
    
    
    !call mpr_read_hru to read in the .asc file defining HRU map (or zones requiring different parameters)
    call mpr_read_hru(fpath_HRUmap, &
        HRU_map, &
        nHRUs, &
        stats_HRU)

    !call mpr_read_BasPred to read all required basin predictors, and clip to HRU map size
    call mpr_read_BasPred(folder_input, &
        MPRfolder, &
        fnames_baspred, &
        req_bp, &
        pedo_tf_all, &
        bp_maps, &
        sand_0_10, &
        silt_0_10, &
        clay_0_10, &
        orgm_0_10, &
        soilmusiddata, &
        stats_hru, &
        hru_map, &
        save_bp_maps, &
        bp_root_depths,&
        B_Ks_table, &
        T_Sy_table)
    
    !YZ 2023 soil classifcation
    call get_soilclass(sand_0_10, silt_0_10,clay_0_10,soil_class)
    write(999,*) 'Soil classification for B&Ks finished!'
    
    write(999,*) ''
    write(999,*) 'Generating MPR global parameters'

    !call mpr_init_gparams to generate n global parameter sets
    call mpr_init_gparams(n_pm_maps, &
        start_seed, &
        pedo_tf_all, &
        Params_range, &
        n_gp_all, &
        glob_pms_SZM, &
        glob_pms_LnTo, &
        glob_pms_SRmax, &
        glob_pms_SRinit, &
        glob_pms_CHV, &
        glob_pms_Td, &
        glob_pms_Smax, &
        glob_pms_B, &
        glob_pms_Ks, &
        bp_root_depths,&
        B_Ks_table)

!    !call mpr_generate_parammap_pert to generate n modifications to parameter maps.
!    call mpr_generate_parammap_pert(param_repeatruns, &
!        parammap_pert_ranges, &
!        parammap_add_mult)

    !write global parameters in output folder - and random parameters if repeatruns.
    folder_outfull = (trim(MPRfolder) // trim(folder_output))
    call mpr_write_gparams(n_pm_maps, &
        pedo_tf_all, &
        folder_output, &
        n_gp_all, &
        glob_pms_SZM, &
        glob_pms_LnTo, &
        glob_pms_SRmax, &
        glob_pms_SRinit, &
        glob_pms_CHV, &
        glob_pms_Td, &
        glob_pms_Smax, &
        glob_pms_B, &
        glob_pms_Ks, &
        param_repeatruns, &
        parammap_add_mult, &
        fpath_gp_file)
    
    !YZ June 2023
    !Clip&prepare GW model data
    write(999,*) ''
    write(999,*) 'Clip/prepare or read GW model data'
    
    
    !If it is the first time run
    !calculate GW input/output weight
    write(999,*) ''
    write(999,*) 'calculate HRU-GW input/output weight'
    
  
    tem_name=fname_filemgr(len_trim(fname_filemgr)-10+1:len_trim(fname_filemgr)-10+1)
    
    !if the catchment name is 6-digit number
    if (ichar(tem_name(1:1)) >= ichar('0') .and. ichar(tem_name(1:1)) <= ichar('9')) then
       cats_name=trim(fname_filemgr(len_trim(fname_filemgr)-10+1:len_trim(fname_filemgr)-5+1))
    else
        !The catchment name is 5-digit number
        cats_name=trim(fname_filemgr(len_trim(fname_filemgr)-9+1:len_trim(fname_filemgr)-5+1))
    end if

    print*,'GW clipped data catchment name: ',cats_name
    call GW_read_data(cats_name,hru_map,stats_HRU,folder_input,&
                save_GW_maps, fnames_baspred,&
                GW_data,GW_Sy_clip,GW_T_clip,GW_hi_ini,GW_hNew_ini,GW_runoff_ini,GEO_index,&
                GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w,&
                GWgauge_meta,GW_print_control,set_GW_buffer)

    
    !YZ Feb 2024
    !Sample GW grid T/Sy 
    write(999,*) ''
    write(999,*) 'Sampling GW grid T/Sy parameters'
    call init_gp_T(pedo_tf_all(10),glob_pms_T, n_pm_maps, start_seed, n_gp_all,&
        T_Sy_table,GW_T_clip,GW_T)
    call init_gp_Sy(pedo_tf_all(11),glob_pms_Sy, n_pm_maps, start_seed, n_gp_all,&
        T_Sy_table,GW_Sy_clip,GW_Sy)
    
    !transmissvity T
    if (pedo_tf_all(10)==2) then
    !if (size(glob_pms_T,2)>1) then !for sensitivity test
        filename=trim(out_dir_path)//'Tpar.txt'
        open(unit=101,file=filename,status='unknown')
        write(101, '(A)') '! T transmissvity sampled for 5000 sims'
        write (line_format, '(A,I0, A)') '(',size(glob_pms_T, 1),'(1X, f0.3), A)'
        do i = 1, size(glob_pms_T, 2)
            write(101, line_format) glob_pms_T(:, i)
        end do
        close(101)
    
    end if

    !Sy
    if (pedo_tf_all(11)==2) then
    !if (size(glob_pms_Sy,2)>1) then !for sensitivity test
        filename=trim(out_dir_path)//'Sypar.txt'
        open(unit=102,file=filename,status='unknown')
        write(102, '(A)') '! Sy sampled for 5000 sims'
        write (line_format, '(A,I0, A)') '(',size(glob_pms_Sy, 1),'(1X, f0.5), A)'
        do i = 1, size(glob_pms_Sy, 2)
            write(102, line_format) glob_pms_Sy(:, i)
        end do
        close(102)

    end if

    ! Loop through for each simulation
    write(999,*) ''
    write(999,*) 'Looping through simulations...'


    do i_mc = 1, n_pm_maps  !previously numsim

        !print *,
        print *, '**************************************************************************'
        print *, 'Running sim ', i_mc, ' out of ',n_pm_maps
        print *, '**************************************************************************'

        !Apply pedotransfer functions to get parameter values at 50m resolution
        call mpr_pedotransfer(i_mc, &
            n_pm_maps, &
            pedo_tf_all, &
            save_pm_maps, &
            start_seed, &
            folder_output, &
            glob_pms_SZM, &
            glob_pms_LnTo, &
            glob_pms_SRmax, &
            glob_pms_SRinit, &
            glob_pms_CHV, &
            glob_pms_Td, &
            glob_pms_Smax,&
            glob_pms_B,&
            glob_pms_Ks,&
            soil_class,&
            soilmusiddata, &
            bp_maps , &
            pm_map_all, &
            stats_HRU, &
            HRU_map,&
            pm_map_Drz)


        !upscale parameters to get parameter sets per parameter map (HRU map)
        call mpr_upscale(HRU_map, &
            nHRUs, &
            pm_map_all, &
            up_pms, &
            pedo_tf_all, &
            save_pm_maps, &
            stats_HRU, &
            i_mc, &
            start_seed, &
            folder_output,&
            pm_map_Drz,&
            up_pms_Drz)

        

        !write upscaled parameters into DYMOND expected format - mcpar
        !YZ 2023 mcpar Drz
        allocate(mcpar(nHRUs,10))
        mcpar(1:nHRUs,1:9) = up_pms(1:nHRUs,2:10)
        mcpar(1:nHRUs,10) = up_pms_Drz(1:nHRUs,1)
        NTT=1
        WT=0.8              !kinematic solution parameters
        ACC=0.000000001

        !YZ 2024
        !Generate T/Sy geo spatial map GW_T/Sy
        !according to glob_pms_t/Sy
        !if T/Sy is fixed estimated value, then glob_pms_T=0
        if (pedo_tf_all(10)==2) then
        !if (size(glob_pms_T,2)>1) then
            print*,'Updating T'
            allocate(GW_Ttemp(size(GEO_index,1),size(GEO_index,2)))
            GW_Ttemp = -1 !easy to find error if it's negative value
            call Tglob_pmsTomap(glob_pms_T,i_mc,Geo_index,&
                GW_Ttemp,T_Sy_table)
            
            GW_T = GW_Ttemp
            deallocate(GW_Ttemp)

        end if 
          
          if (pedo_tf_all(11)==2) then
         ! if (size(glob_pms_Sy,2)>1) then
            print*,'Updating Sy'
            allocate(GW_Sytemp(size(GEO_index,1),size(GEO_index,2)))
            !easy to find error if it's negative value
            GW_Sytemp = -1
            call Syglob_pmsTomap(glob_pms_Sy,i_mc,Geo_index,&
                GW_Sytemp,T_Sy_table)
            
            GW_Sy = GW_Sytemp
            deallocate(GW_Sytemp)

        end if 



    ! Open up output files
        call file_open (i_mc+start_seed-1, &
            out_dir_path, pf_data,GW_print_control)
        print*,'Actual simulation No.:',i_mc+start_seed-1

        ! ADD AS A NON-MPR VERSION ONLY
        !        print *, 'Running sim ', i_mc

        ! ADD AS A NON-MPR VERSION ONLY
        ! Get a new set of random parameters
        !        call genpar (mcpar, &
        !            mcpar_ll, &
        !            mcpar_ul, &
        !            num_mcpar, &
        !            num_par_types)
        
        !15 Nov 2023 YZ
        !Initialization GW intial table, GW new table, GW runoff matrix
        GW_hi=GW_hi_ini
        GW_hNew=GW_hNew_ini
        GW_runoff= GW_runoff_ini                                    
        call mainloop (all_pm_names, &
            nac, &
            nstep, &
            num_rivers, &
            acc, &
            dt, &
            dtt, &
            dyna_hru, &
            mcpar, &
            ntt, &
            node_to_flow_mapping, &
            nHRUs, & !ACTUALLY NUMBER OF PARAMETER GRIDS
            pe_step, &
            qobs_riv_step_start, &
            r_gau_step, &
            rivers, &
            route_riv, &
            route_tdh, &
            sum_ac_riv, &
            wt,&
            GW_data,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff,&
            stats_hru,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w,&
            GW_print_control,GW_maxiteration,tolerance)
            
        
            

        close (40)
        close (41)
        if (GW_print_control==1) then
            close (43)
        end if

    end do
    
    close(999)
    !
    !============================================================
    !  END MAIN MC LOOP
    !============================================================
    !
    stop

end program main
