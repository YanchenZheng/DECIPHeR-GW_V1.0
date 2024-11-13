!! Module containing main subroutine for MPR

!! Rosie Lane - 3rd November 2017

module mpr_main
contains




        ! **********************************************************************************************
        ! subroutine: 1.mpr_read_controls: reads control file and param file
        ! **********************************************************************************************
        ! This function checks that all required MPR files are in correct paths, prints error messages
        ! otherwise. Then it reads the parameter and control file - reading each variable in from the
        ! control file seperately so that the order they are entered is irrelevant.
        ! Returns values of user_inputted requirements.
        ! 01/06/2018 - also reads in filepath from MPR_filenames.dat
    subroutine mpr_read_controls(MPRfolder, &    !key filepaths and folder locations
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
        pedo_tf_all, &       !list of pedo-transfer functions selected for each parameter
        param_repeatruns, &
        parammap_pert_ranges, &
        fpath_gp_file, &
        save_GW_maps, &
        set_GW_buffer,GW_maxiteration,GW_tolerance)    !numbver of runs to do per parameter map created - modifying parameter map.

        use mpr_control
        use mpr_errors

        implicit none

        !DECLARE DUMMY VARIABLES
        !key filepaths and folder locations
        character(len=1024),intent(out)                 :: MPRfolder        !location of MPR folder containing INPUT and OUTPUT folders, and MPR_control.dat file
        character(len=1024),intent(out)                 :: fpath_HRUmap     !filepath pointing to HRU .asc file created in DTA
        character(len=1024),intent(out)                 :: folder_input     !folder containing Basin Predictor files - INPUT folder within MPR main folder
        character(len=1024),intent(out)                 :: folder_output    !output folder to store parameter input files - OUTPUT folder within MPR main folder
        character(len=1024),intent(in )                 :: fname_control    !full path to control file
        character(len=1024),intent(in )                 :: fname_filemgr    !full path to file manager file
        character(len=1024),dimension(27),intent(out)   :: fnames_baspred   !list of names of basin predictor files.

        !parameter min, default and max values from parameter file
        double precision, dimension(9,3),intent(out)  :: Params_range

        !variables in the control file
        integer,intent(out)                 ::  n_pm_maps       !number of parameter files to create
        integer,intent(out)                 ::  start_seed      !starting seed for generation of global params
        integer,intent(out)                 ::  save_pm_maps    !set to 1 if user wants to save parameter maps
        integer,intent(out)                 ::  save_bp_maps    !set to 1 if user wants to save basin predictor maps
        integer,intent(out)                 ::  save_GW_maps    !set to 1 if user wants to save GW data maps
        integer,intent(out)                 ::  set_GW_buffer
        integer,intent(out)                 ::  GW_maxiteration
        integer,intent(out)        ::  GW_tolerance
        !variables combining results from the control file
        integer, dimension(11),intent(out)   ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax

        !parameter map perturbation information
        integer, dimension(7),intent(out)               :: param_repeatruns !number of parameter map modifications and repeat runs for SZM, LnTo, SRmax,SRinit, CHV, Td, Smax
        double precision, dimension(7,4), intent(out)   :: parammap_pert_ranges !returned ranges for additions and multiplications to parameter map.

        !DECLARE LOCAL VARIABLES
        character(len=1024)     ::  fpath_PmSettings   !file containing parameter details
        character(len=1024)     ::  fpath_PmPertSettings  !file containing parameter map perturbation details.
        character(len=1024)     ::  fpath_gp_file
        character(len=1024)     ::  temp_name
        character(len=13)       ::  comment
        character(len=1024)     ::  folder_settings
        character(len=1024)     ::  def_msg = ""
        integer                 ::  pedo_tf_SZM
        integer                 ::  pedo_tf_LnTo
        integer                 ::  pedo_tf_SRmax   !pedo-transfer equation numbers for all params
        integer                 ::  pedo_tf_SRinit
        integer                 ::  pedo_tf_CHV
        integer                 ::  pedo_tf_Td
        integer                 ::  pedo_tf_Smax
        integer                 ::  i
        
        !YZ 2023
        integer                 ::  pedo_tf_B
        integer                 ::  pedo_tf_Ks
        integer                 ::  pedo_tf_T
        integer                 ::  pedo_tf_Sy

        !print *,
        print *, '**************************************************************************'
        print *, '0. Get filepaths from filemanager'
        print *, '**************************************************************************'

        !print *,
        print *, '!Read filenames from the filemanager control file:'
        CALL read_control_file_char(fname_filemgr, "dir_root....", MPRfolder,def_msg)
        CALL read_control_file_char(fname_filemgr, "dir_input...", folder_input,def_msg)
        CALL read_control_file_char(fname_filemgr, "dir_output..", folder_output,def_msg)
        CALL read_control_file_char(fname_filemgr, "dir_settings", folder_settings,def_msg)
        CALL read_control_file_char(fname_filemgr, "name_hru_map", temp_name,def_msg)
        fpath_HRUmap = trim(folder_input) // trim(temp_name)
        CALL read_control_file_char(fname_filemgr, "name_pm_stng", temp_name,def_msg)
        fpath_PmSettings = trim(folder_settings) // trim(temp_name)
        CALL read_control_file_char(fname_filemgr, "name_pm_pert", temp_name,def_msg)
        fpath_PmPertSettings = trim(folder_settings) // trim(temp_name)
        CALL read_control_file_char(fname_filemgr, "name_gp_file", temp_name,def_msg)
        fpath_gp_file = trim(folder_output) // trim(temp_name)

        !print *,
        print *, '!Read basin predictor filenames from filemanager control file:'
        CALL read_control_file_char(fname_filemgr, "lnto_sand_d1", fnames_baspred(1),def_msg)
        CALL read_control_file_char(fname_filemgr, "lnto_silt_d1", fnames_baspred(2),def_msg)
        CALL read_control_file_char(fname_filemgr, "lnto_clay_d1", fnames_baspred(3),def_msg)
        CALL read_control_file_char(fname_filemgr, "lnto_orgm_d1", fnames_baspred(4),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_musids..", fnames_baspred(5),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_mid_nde.", fnames_baspred(6),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_table_d1", fnames_baspred(7),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_table_d2", fnames_baspred(8),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_table_d3", fnames_baspred(9),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_table_d4", fnames_baspred(10),def_msg)
        CALL read_control_file_char(fname_filemgr, "szm_table_d5", fnames_baspred(11),def_msg)
        CALL read_control_file_char(fname_filemgr, "srmax_lcm...", fnames_baspred(12),def_msg)
        CALL read_control_file_char(fname_filemgr, "srmax_table.", fnames_baspred(13),def_msg)
        CALL read_control_file_char(fname_filemgr, "smax_d2b_map", fnames_baspred(14),def_msg)
        CALL read_control_file_char(fname_filemgr, "lnto_bulk_d1", fnames_baspred(15),def_msg)
        CALL read_control_file_char(fname_filemgr, "lnto_isorg..", fnames_baspred(16),def_msg)
        CALL read_control_file_char(fname_filemgr, "sr_smax_poro", fnames_baspred(17),def_msg)
        CALL read_control_file_char(fname_filemgr, "hydrogeology", fnames_baspred(18),def_msg)
        CALL read_control_file_char(fname_filemgr, "B_Ks_table..", fnames_baspred(19),def_msg) !YZ 2023
        CALL read_control_file_char(fname_filemgr, "T_Sy_table..", fnames_baspred(20),def_msg) !YZ 2024
        !YZ May 2024
        !GW related data inputs
        CALL read_control_file_char(fname_filemgr, "GW_topo_clip", fnames_baspred(21),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_mask_clip", fnames_baspred(22),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_Sy_clip..", fnames_baspred(23),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_T_clip...", fnames_baspred(24),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_Hini_clip", fnames_baspred(25),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_GEO_clip.", fnames_baspred(26),def_msg)
        CALL read_control_file_char(fname_filemgr, "GW_obsgauge.", fnames_baspred(27),def_msg)

!        !OLD WAY OF READING FILENAMES FROM THE FILEMANAGER CONTROL FILE
!        !KEEP THIS CODE, AS POSSIBLY MORE ROBUST THAN ABOVE.
!        print *, 'fname_filemanager = ',trim(fname_filemgr)
!        CALL file_open_err(fname_filemgr,2)
!        READ(2,*) !ignore header line
!        READ(2,525) comment, MPRfolder
!        READ(2,525) comment, folder_input
!        READ(2,525) comment, folder_output
!        READ(2,525) comment, folder_settings
!        READ(2,525) comment, temp_name
!        fpath_HRUmap = trim(folder_input) // trim(temp_name)
!        READ(2,525) comment, temp_name
!        fpath_PmSettings = trim(folder_settings) // trim(temp_name)
!        READ(2,525) comment, temp_name
!        fpath_PmPertSettings = trim(folder_settings) // trim(temp_name)
!        !! Read filenames of all basin predictors.
!        READ(2, *)
!        READ(2, *)
!        DO i = 1,11
!            READ(2,525) comment, fnames_baspred(i)
!        END DO
!        CLOSE(2)

        !check all files exist and print useful info to screen
        !print *,
        !print *, 'Double check these locations: '
        !print *, 'Param map : ',trim(fpath_HRUmap)
        !print *, 'Param settings : ',trim(fpath_PmSettings)
        !print *, 'fname_control = ',trim(fname_control)
        !print *, 'Param map settings = ',trim(fpath_PmPertSettings)
        CALL check_files_exist(MPRfolder, fpath_HRUmap, folder_input, folder_output, fname_control,fpath_PmSettings)

        !print *,
        print *, '**************************************************************************'
        print *, '1. read parameter and control file'
        print *, '**************************************************************************'

        !parameter file
        CALL read_param_file(fpath_PmSettings,Params_range(1,:),Params_range(2,:),Params_range(3,:),&
            Params_range(4,:),Params_range(5,:),Params_range(6,:),Params_range(7,:),&
            Params_range(8,:),Params_range(9,:))

        !parameter map perturbation settings file
        CALL read_parammap_pert_file(fpath_PmPertSettings,parammap_pert_ranges)

        !mpr control file.
        CALL read_control_file_int(fname_control, "n_pm_maps", n_pm_maps,1)
        CALL read_control_file_int(fname_control, "start_seed", start_seed,1)
        CALL read_control_file_int(fname_control, "save_pm_maps", save_pm_maps,0)
        CALL read_control_file_int(fname_control, "save_bp_maps", save_bp_maps,0)
        CALL read_control_file_int(fname_control, "pedo_tf_SZM", pedo_tf_SZM,0)
        CALL read_control_file_int(fname_control, "pedo_tf_LnTo", pedo_tf_LnTo,0)
        CALL read_control_file_int(fname_control, "pedo_tf_SRmax", pedo_tf_SRmax,0)
        CALL read_control_file_int(fname_control, "pedo_tf_SRinit", pedo_tf_SRinit,0)
        CALL read_control_file_int(fname_control, "pedo_tf_CHV", pedo_tf_CHV,0)
        CALL read_control_file_int(fname_control, "pedo_tf_Td", pedo_tf_Td,0)
        CALL read_control_file_int(fname_control, "pedo_tf_Smax", pedo_tf_Smax,0)
        CALL read_control_file_int(fname_control, "n_runs_SZM",param_repeatruns(1),0)
        CALL read_control_file_int(fname_control, "n_runs_LnTo",param_repeatruns(2),0)
        CALL read_control_file_int(fname_control, "n_runs_SRmax",param_repeatruns(3),0)
        CALL read_control_file_int(fname_control, "n_runs_SRinit",param_repeatruns(4),0)
        CALL read_control_file_int(fname_control, "n_runs_CHV",param_repeatruns(5),0)
        CALL read_control_file_int(fname_control, "n_runs_Td",param_repeatruns(6),0)
        CALL read_control_file_int(fname_control, "n_runs_Smax",param_repeatruns(7),0)
        !YZ May 2024
        CALL read_control_file_int(fname_control, "save_GW_maps",save_GW_maps,0)
        CALL read_control_file_int(fname_control, "set_GW_buffer",set_GW_buffer,3)
        CALL read_control_file_int(fname_control, "GW_maxiteration",GW_maxiteration,1500)
        CALL read_control_file_int(fname_control, "GW_tolerance",GW_tolerance,1)
        !save_GW_maps=1

        pedo_tf_all(1:4) = (/pedo_tf_SZM, pedo_tf_LnTo, pedo_tf_SRmax, pedo_tf_SRinit /)
        pedo_tf_all(5:7) = (/pedo_tf_CHV, pedo_tf_Td, pedo_tf_Smax /)

        !YZ 2023
        CALL read_control_file_int(fname_control, "pedo_tf_B", pedo_tf_B,0)
        CALL read_control_file_int(fname_control, "pedo_tf_Ks", pedo_tf_Ks,0)
        pedo_tf_all(8:9) = (/pedo_tf_B, pedo_tf_Ks /)

        !YZ 2024
        CALL read_control_file_int(fname_control, "pedo_tf_T", pedo_tf_T,0)
        CALL read_control_file_int(fname_control, "pedo_tf_Sy", pedo_tf_Sy,0)
        pedo_tf_all(10:11) = (/pedo_tf_T, pedo_tf_Sy /)

        !only allow repeat runs for parameters where MPR is being applied.
        DO i=1,7
            IF (param_repeatruns(i) > 0) THEN
                IF (pedo_tf_all(i) < 2) THEN
                    print *, 'IMPORTANT: parameter ', i,' cannot have repeat runs as no spatial field applied'
                    param_repeatruns(i) = 0
                END IF
            END IF
        END DO

    525     format(a13,1X,a700)

    end subroutine mpr_read_controls




        ! **********************************************************************************************
        ! subroutine: 2. mpr_read_hru: Reads the HRU/ spatial parameter map
        ! **********************************************************************************************
        ! This function reads the HRU ascii grid and header, from file specified in fpath_HRUmap
        ! It saves key information from the header into stats_HRU
        ! It makes sure all nodata values are set to -9999 for consistency.
    subroutine mpr_read_hru(fpath_HRUmap, &
        HRU_map, &
        nHRUs, &
        stats_HRU)

        use dyna_utility

        implicit none

        !DECLARE DUMMY VARIABLES
        character(len=1024), intent(in) :: fpath_HRUmap     !filepath pointing to HRU .asc file created in DTA

        !variables defining HRU map and parameter maps
        double precision, allocatable, dimension(:,:),intent(out) :: HRU_map
        integer,intent(out)                                       :: nHRUs            !number of HRUs
        double precision, dimension(4),intent(out)                :: stats_HRU        !variable comprising HRU map xll, yll, cellsize & nodata value

        !DECLARE LOCAL VARIABLES
        integer                         :: ncols_HRU
        integer                         :: nrows_HRU        !all variables relating to size and location
        double precision                :: xllcorner_HRU    !of HRU map created during DTA
        double precision                :: yllcorner_HRU
        double precision                :: cellsize_HRU
        double precision                :: nodata_HRU
        integer                         :: ncol
        integer                         :: nrow             !all variables relating to size and location
        double precision                :: xll              !of param maps - assumed same as HRU map
        double precision                :: yll
        double precision                :: cellsize

        !print *,
        print *, '**************************************************************************'
        print *, '2. Read HRU map'
        print *, '**************************************************************************'

        print *, 'Reading map from: ',trim(fpath_HRUmap)

        call read_ascii_grid(fpath_HRUmap, HRU_map, ncols_HRU, nrows_HRU, xllcorner_HRU, yllcorner_HRU, cellsize_HRU, nodata_HRU)

        nHRUs = maxval(int(HRU_map))
        stats_HRU(1:4) = (/ xllcorner_HRU, yllcorner_HRU, cellsize_HRU, nodata_HRU/)

        !stats values of hru map - all other maps will be converted to these
        nrow = nrows_HRU
        ncol = ncols_HRU
        xll = xllcorner_HRU
        yll = yllcorner_HRU
        cellsize = cellsize_HRU

        print *, 'Number of HRUs found: ',nHRUs
        print *, 'Rows in HRU file: ',nrows_HRU
        print *, 'Cols in HRU file: ',ncols_HRU

        !convert nodata value to -9999 to avoid confusion
        IF (stats_HRU(4) /= -9999) THEN
            stats_HRU(4) = -9999
            print *,'old nodata value: ',int(nodata_HRU),' being converted to ',int(stats_HRU(4))
            call set_nodata_value( HRU_map, nodata_HRU, stats_HRU(4))
            nodata_HRU = -9999
        ELSE
            print *, 'nodata value: ', nodata_HRU
        END IF

    end subroutine mpr_read_hru





        ! **********************************************************************************************
        ! subroutine: 3. mpr_read_BasPred: Reads required Basin Predictor maps and trims to HRU extent
        ! **********************************************************************************************
        ! a. checks which basin predictors are needed for selected pedo-transfer functions
        ! b. reads in only those maps
        !   - maps are clipped to same extent as HRU map when being read-in
        !   - maps can be saved if specified in MPR_control file
    subroutine mpr_read_BasPred(folder_input, &
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
        bp_root_depths, &
        B_Ks_table, &
        T_Sy_table)

        use dyna_utility
        use mpr_extract_BasPred
        use mpr_SZM
        use mpr_LnTo
        use mpr_SRmaxGW
        use mpr_SRinit
        use mpr_CHV
        use mpr_Td
        use mpr_Smax
        

        implicit none

        !DECLARE DUMMY VARIABLES
        !key filepaths
        character(len=1024), intent(in)                 :: MPRfolder        !location of MPR folder containing INPUT, SETTINGS and OUTPUT folders
        character(len=1024), intent(in)                 :: folder_input   !folder containing Basin Predictor files
        character(len=1024), dimension(27),intent(in)   :: fnames_baspred

        !specific requirements based on user input to control file
        integer, intent(in)                 ::  save_bp_maps    !set to 1 if user wants to save parameter maps
        !logical, dimension(6,5), intent(out)::  req_soils       !basin predictor maps needed - logical arrays with true values if maps required -Sand, silt, clay and OM rows
        logical, dimension(20), intent(out)::  req_bp          !basin predictor maps needed - logical arrays with true values if maps required -same dimensions as fnames_baspred.
        integer, dimension(11),intent(inout) ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax

        !basin predictor maps
        real, allocatable, dimension(:,:,:),intent(out)  :: bp_maps     !big array of all basin predictor maps
        real, allocatable, dimension(:,:),intent(out)    :: sand_0_10
        real, allocatable, dimension(:,:),intent(out)    :: silt_0_10
        real, allocatable, dimension(:,:),intent(out)    :: clay_0_10
        real, allocatable, dimension(:,:),intent(out)    :: orgm_0_10
        double precision, allocatable, dimension(:,:,:), intent(out) :: soilmusiddata  !SZM: tables of musid data
        real, allocatable, dimension(:,:), intent(out)   :: bp_root_depths    !SRmax: landcover class rooting depth table
        double precision, allocatable, dimension(:,:), intent(out)   :: B_Ks_table !YZ 2023
        double precision, allocatable, dimension(:,:), intent(out)   :: T_Sy_table !YZ 2024

        !HRU map and properties
        double precision, allocatable, dimension(:,:),intent(in) :: HRU_map
        double precision, dimension(4), intent(in)               :: stats_HRU        !variable comprising HRU map xll, yll, cellsize & nodata value

        !DECLARE LOCAL VARIABLES
        integer :: i

        !print *,
        print *, '**************************************************************************'
        print *, '3. Read required basin predictor maps'
        print *, '**************************************************************************'
        print *, 'Reading maps from: ',trim(folder_input)
        !print *,

!        ! assume no basin predictors are required unless specified otherwise
!        req_soils(1, 1:5) = (/ .false. , .false., .false., .false., .false. /) !sand (5 depth profiles)
!        req_soils(2, 1:5) = (/ .false. , .false., .false., .false., .false. /) !silt
!        req_soils(3, 1:5) = (/ .false. , .false., .false., .false., .false. /) !clay
!        req_soils(4, 1:5) = (/ .false. , .false., .false., .false., .false. /) !orgm
!        req_soils(5, 1:5) = (/ .false. , .false., .false., .false., .false. /) !Ksat form of decline
!        req_soils(6, 1:5) = (/ .false. , .false., .false., .false., .false. /) !Soils depth data for SZM - table not ascii. Multiple per MUSID.
!        req_soils(6, 1:5) = (/ .false. , .false., .false., .false., .false. /) !Soils depth data for SZM - table not ascii. Multiple per MUSID.

        DO i = 1,14
            req_bp(i) = .false.
        END DO

        ! Check which basin predictors are required
        ! get_BasPred_req returns a logical list of the required soils maps req_soils(type,depth),
        ! which contains .true. where soils maps are required
        CALL get_BasPred_SZM(pedo_tf_all(1), req_bp)
        CALL get_BasPred_LnTo(pedo_tf_all(2), req_bp)
        CALL get_BasPred_SRmax(pedo_tf_all(3), req_bp)
        CALL get_BasPred_SRinit(pedo_tf_all(4), req_bp)
        CALL get_BasPred_CHV(pedo_tf_all(5), req_bp)
        CALL get_BasPred_Td(pedo_tf_all(6), req_bp)
        CALL get_BasPred_Smax(pedo_tf_all(7), req_bp)

        !YZ Oct 2023 read silt data for soil texture classification
        req_bp(2) = .true. !read silt data
        req_bp(19) = .true. !read B_Ks_lookup table
        !YZ Feb 2024 read T_Sy lookup table
        req_bp(20) = .true.

        ! read_req_soils reads in any basin predictor maps required for user-selected pedo transfer functions
        CALL read_req_soils(folder_input,&
            MPRfolder,&
            fnames_baspred,&
            req_bp,&
            stats_HRU,&
            HRU_map,&
            save_bp_maps,&
            bp_maps, &
            soilmusiddata, &
            bp_root_depths,&
            B_Ks_table,&
            T_Sy_table)

        sand_0_10 = bp_maps(:,:,1)
        silt_0_10 = bp_maps(:,:,2)
        clay_0_10 = bp_maps(:,:,3)
        orgm_0_10 = bp_maps(:,:,4)


    end subroutine mpr_read_BasPred



        ! **********************************************************************************************
        ! subroutine: 4. mpr_init_gparams: initialises global parameters using random num generator
        ! **********************************************************************************************
        ! a. checks how many global params are needed, given pedo-transfer functions selected
        ! b. initialises global parameters within given bounds
        ! c. generates n sets of parameters, starting at given start seed
    subroutine mpr_init_gparams(n_pm_maps, &
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
        rd,&
        B_Ks_table)

        use dyna_utility
        use mpr_extract_BasPred
        use mpr_SZM
        use mpr_LnTo
        use mpr_SRmaxGW
        use mpr_SRinit
        use mpr_CHV
        use mpr_Td
        use mpr_Smax
        !YZ 2023
        use mpr_B_KsGW

        implicit none

        !DECLARE DUMMY VARIABLES
        !variables from the control file
        integer,intent(in)                ::  n_pm_maps       !number of parameter files to create
        integer,intent(in)                ::  start_seed      !starting seed for generation of global params
        integer, dimension(11),intent(in)   ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax


        !parameter min, default and max values from parameter file
        double precision, dimension(9,3),intent(in)  :: Params_range

        !global parameters
        integer, dimension(11),intent(out)                          :: n_gp_all     !vector of number of global params for all params
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_SZM     !dim(glob pm number, 1:n_pm_maps)
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_LnTo
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_SRmax
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_SRinit
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_CHV
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_Td
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_Smax
        
        !YZ 2023
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_B
        double precision, allocatable, dimension(:,:),intent(out)   :: glob_pms_Ks
        !special parameter information files
        real, allocatable, dimension(:,:),intent(in)    :: rd !rooting depth table bp_root_depths
        double precision, allocatable, dimension(:,:),intent(in)    :: B_Ks_table ! soil class lookup table

        !print *,
        print *, '**************************************************************************'
        print *, '4. Initialise global parameters within given bounds'
        print *, '**************************************************************************'

        n_gp_all = (/ 0,0,0,0,0,0,0,0,0,0,0 /)  !initialise vector storing the number of global parameters needed for each parameter.
        !        !call init_gp_ functions to fill in the glob_pms_ arrays, glob_pms(pm number, 1:n)

        if (pedo_tf_all(1)==4) THEN
            CALL init_gp_SZM(2,glob_pms_SZM, n_pm_maps, start_seed, n_gp_all, Params_range(1,:))
        ELSE
            CALL init_gp_SZM(pedo_tf_all(1),glob_pms_SZM, n_pm_maps, start_seed, n_gp_all, Params_range(1,:))
        END IF
        CALL init_gp_LnTo(pedo_tf_all(2),glob_pms_LnTo, n_pm_maps, start_seed, n_gp_all, Params_range(2,:))
        CALL init_gp_SRmax(pedo_tf_all(3),glob_pms_SRmax, n_pm_maps, start_seed, n_gp_all, Params_range(3,:),rd)
        CALL init_gp_SRinit(pedo_tf_all(4),glob_pms_SRinit, n_pm_maps, start_seed, n_gp_all, Params_range(4,:))
        CALL init_gp_CHV(pedo_tf_all(5),glob_pms_CHV, n_pm_maps, start_seed, n_gp_all,Params_range(5,:))
        CALL init_gp_Td(pedo_tf_all(6),glob_pms_Td, n_pm_maps, start_seed, n_gp_all,Params_range(6,:))
        CALL init_gp_Smax(pedo_tf_all(7),glob_pms_Smax, n_pm_maps, start_seed, n_gp_all,Params_range(7,:))
        
        !YZ 2023
        CALL init_gp_B(pedo_tf_all(8),glob_pms_B, n_pm_maps, start_seed, n_gp_all,Params_range(8,:),B_Ks_table)
        CALL init_gp_Ks(pedo_tf_all(9),glob_pms_Ks, n_pm_maps, start_seed, n_gp_all,Params_range(9,:),B_Ks_table)

    end subroutine mpr_init_gparams




        ! **********************************************************************************************
        ! subroutine: 4b. mpr_write_gparams: writes initialised global parameters in output file
        ! **********************************************************************************************
        !
    subroutine mpr_write_gparams(n_pm_maps, &
        pedo_tf_all, &
        out_dir_path, &
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

        use dyna_utility
        !use mpr_extract_BasPred
        !use mpr_SZM
        !use mpr_LnTo
        !use mpr_SRmax
        !use mpr_SRinit
        !use mpr_CHV
        !use mpr_Td
        !use mpr_Smax

        implicit none

        !DECLARE DUMMY VARIABLES
        !variables from the control file
        integer,intent(in)                  ::  n_pm_maps       !number of parameter files
        integer, dimension(11),intent(in)    ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax
        character(900), intent(in)          ::  out_dir_path
        character(len=1024),intent(in)      :: fpath_gp_file

        !global parameters
        integer, dimension(11),intent(in)                           :: n_gp_all     !vector of number of global params for all params
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_SZM     !dim(glob pm number, 1:n_pm_maps)
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_LnTo
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_SRmax
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_SRinit
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_CHV
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_Td
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_Smax
        
        !YZ 2023
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_B
        double precision, allocatable, dimension(:,:),intent(in)   :: glob_pms_Ks

        !variables for parameter map perturbation runs.
        integer, dimension(7), intent(in) :: param_repeatruns   !number of perturbed repeats for each parameter.
        double precision, allocatable, dimension(:,:,:), intent(in) :: parammap_add_mult

        !local parameters
        character(900) :: filename
        character(1000000)  :: numformat
        double precision, allocatable, dimension(:,:,:) :: glob_pms_all
        integer :: i, j
        integer :: max_ngp
        character(1) :: c

        max_ngp = maxval(n_gp_all(:))

        !put all global parameters into one large array
        allocate(glob_pms_all(max_ngp,n_pm_maps,9))

        glob_pms_all(1:n_gp_all(1),:,1) = glob_pms_SZM
        glob_pms_all(1:n_gp_all(2),:,2) = glob_pms_LnTo
        glob_pms_all(1:n_gp_all(3),:,3) = glob_pms_SRmax
        glob_pms_all(1:n_gp_all(4),:,4) = glob_pms_SRinit
        glob_pms_all(1:n_gp_all(5),:,5) = glob_pms_CHV
        glob_pms_all(1:n_gp_all(6),:,6) = glob_pms_Td
        glob_pms_all(1:n_gp_all(7),:,7) = glob_pms_Smax
        glob_pms_all(1:n_gp_all(8),:,8) = glob_pms_B
        glob_pms_all(1:n_gp_all(9),:,9) = glob_pms_Ks

        ! Then open the output file
        open(unit=24,file=fpath_gp_file,status='unknown')

        !write title followed by formatted output for each parameter
        write(24,*) 'Param (szm,lnto,srmax,srinit,chv,td,smax,B,Ks), GP_num, values'
        !numformat = "(I5,I5,F15.10 " // repeat("F15.10 ",n_pm_maps) // ")"
        !numformat = "(I5,I5," // repeat("F15.10 ",n_pm_maps) // ")"
        numformat = "(I2,',',I2" // repeat("',',F15.5",n_pm_maps) // ")"
        do i = 1,9
            do j = 1,n_gp_all(i)
                write(24,numformat) i,j,glob_pms_all(j,:,i)
             end do
        end do
        close(unit=24)

        !numformat

        !write pedotransfer function chosen, global parameter values (-9999 default?)

        ! SAVE PARAMETER PERTURBATIONS IF ALSO NEEDED!
                   !parammap_add_mult(:,i,1) = norm_rand_add
            !parammap_add_mult(:,i,2) = norm_rand_mult
        IF (maxval(param_repeatruns(:))>=1) THEN
            filename = trim(out_dir_path) // 'repeat_perturbations.txt'
            open(unit=24, file=filename, status='unknown')
            write(24,*) 'Param perturbation values:'
            write(24,*) 'Param maps have been perturbed by the following additions/multiplications and re-run'
            write(24,*) 'Parameter, add(1)/multiply(2),n_perturbations, perturbation values'
            numformat = "(I5,I5,I5,F15.10 " // repeat("F15.10 ",maxval(param_repeatruns)) // ")"

            !only include entries for parameters where repeat perturbed runs are going
            do i=1,7
                IF (param_repeatruns(i)>=1) THEN
                    write(24,numformat) i,1,param_repeatruns(i),parammap_add_mult(i,:,1)
                    write(24,numformat) i,2,param_repeatruns(i),parammap_add_mult(i,:,2)
                END IF
            end do
        END IF
        close(unit=24)

        !write(40,FMT0) [pm_class,all_pm_names(1:7)]


   end subroutine mpr_write_gparams





        ! **********************************************************************************************
        ! subroutine: 4c. mpr_generate_parammap_pert
        ! **********************************************************************************************
        ! Aim is to get the random perturbations to parameter maps.
        ! Enables running of one parameter map multiple times.
        ! To reduce computational time.
    subroutine mpr_generate_parammap_pert(param_repeatruns, &
        parammap_pert_ranges, &
        parammap_add_mult)

        use dyna_utility

        implicit none

        !DECLARE DUMMY VARIABLES
        integer,dimension(7),intent(in)                 ::  param_repeatruns    !number of repeat runs for each parameter
        double precision, dimension(7,4), intent(in)    ::  parammap_pert_ranges    !user-defined upper and lower bounds
        double precision, dimension(:,:,:), allocatable,intent(out)   :: parammap_add_mult !param map perturbations - randomly generated between bounds.(nparams,nperturbations,add/mult)

        !DECLARE LOCAL VARIABLES
        integer                 :: i
        integer, dimension(12)  :: seed
        double precision, dimension(:,:), allocatable :: rand_nums
        double precision, dimension(7,3)    :: add_ranges           !add lower bounds, upper bounds and range
        double precision, dimension(7,3)    :: mult_ranges          !multiply lower bounds, upper bounds and range
        double precision, dimension(7)      :: rand_add
        double precision, dimension(7)      :: rand_mult
        double precision, dimension(7)      :: norm_rand_add
        double precision, dimension(7)      :: norm_rand_mult

        !allocate output grid to max size needed - from maximum number of repeat runs needed
        allocate(parammap_add_mult(7,maxval(param_repeatruns),2))
        allocate(rand_nums(14,maxval(param_repeatruns)))

        !get random arrays - one for add and one for multiply.
        Call RANDOM_SEED( GET = seed )
        seed(1:12) = (/ 3,3,3,3,3,3,3,3,3,3,3,3/)
        Call RANDOM_SEED( PUT = seed )
        call RANDOM_NUMBER ( rand_nums )

        !get add and multiply ranges
        add_ranges(:,1) = parammap_pert_ranges(:,1)
        add_ranges(:,2) = parammap_pert_ranges(:,2)
        mult_ranges(:,1) = parammap_pert_ranges(:,3)
        mult_ranges(:,2) = parammap_pert_ranges(:,4)

        add_ranges(:,3) = add_ranges(:,2) - add_ranges(:,1)
        mult_ranges(:,3) = mult_ranges(:,2) - mult_ranges(:,1)

        !normalise random numbers and reformat for output
        DO i = 1,maxval(param_repeatruns)
            rand_add = rand_nums(1:7,i)
            rand_mult = rand_nums(8:14,i)
            norm_rand_add = (rand_add * add_ranges(:,3))+add_ranges(:,1)
            norm_rand_mult = (rand_mult * mult_ranges(:,3))+mult_ranges(:,1)
            parammap_add_mult(:,i,1) = norm_rand_add
            parammap_add_mult(:,i,2) = norm_rand_mult
!            print *,
!            print *, 'parammap_add_mult(1:3,i,1) = ',parammap_add_mult(1:3,i,1)
!            print *, 'parammap_add_mult(1:3,i,2) = ',parammap_add_mult(1:3,i,2)
        END DO


    end subroutine



        ! **********************************************************************************************
        ! subroutine: 5. mpr_pedotransfer: applies tf functions to create parameter maps
        ! **********************************************************************************************
        ! For run n in N parameter maps
        ! a. Apply pedotransfer functions to map n to get parameter maps
        ! b. Save these parameter maps if save_pm_maps specified in the mpr control file
        ! c. output parameter maps
    subroutine mpr_pedotransfer(i_n, &
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

        use dyna_utility
        use mpr_SZM
        use mpr_LnTo
        use mpr_SRmaxGW
        use mpr_SRinit
        use mpr_CHV
        use mpr_Td
        use mpr_Smax
        !YZ 2023
        use mpr_B_KsGW

        implicit none

        !!DECLARE DUMMY VARIABLES
        integer                 :: i_n
        integer                 :: n_pm_maps
        integer, dimension(11)   :: pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax
        integer                 :: save_pm_maps    !set to 1 if user wants to save parameter maps
        integer                 :: start_seed      !starting seed for generation of global params
        character(len=1024)     :: folder_output   !output folder to store parameter input files - OUTPUT folder within MPR main folder

        !global parameters
        double precision, allocatable, dimension(:,:)   :: glob_pms_SZM     !dim(glob pm number, 1:n_pm_maps)
        double precision, allocatable, dimension(:,:)   :: glob_pms_LnTo
        double precision, allocatable, dimension(:,:)   :: glob_pms_SRmax
        double precision, allocatable, dimension(:,:)   :: glob_pms_SRinit
        double precision, allocatable, dimension(:,:)   :: glob_pms_CHV
        double precision, allocatable, dimension(:,:)   :: glob_pms_Td
        double precision, allocatable, dimension(:,:)   :: glob_pms_Smax

        !YZ2023
        double precision, allocatable, dimension(:,:)   :: glob_pms_B
        double precision, allocatable, dimension(:,:)   :: glob_pms_Ks
        real, allocatable, dimension(:,:)   :: soil_class

        !basin predictor maps
        real, allocatable, dimension(:,:,:)    :: bp_maps !bp maps (x,y,sand/silt/clay/orgm/ksat decline,depth_class)
        double precision, allocatable, dimension(:,:,:), intent(in) :: soilmusiddata !musid data tables.

        !parameter maps at the finest scale, created using transfer functions
        double precision, allocatable, dimension(:,:,:)   :: pm_map_all !mp_map_all(x,y,[SZM, LnTo, SRmax, SRinit, CHV, Td, Smax])

        !variables defining HRU map and parameter maps
        double precision, allocatable, dimension(:,:) :: HRU_map
        double precision, dimension(4)  :: stats_HRU !variable comprising HRU map xll, yll, cellsize & nodata value


        !DECLARE LOCAL VARIABLES
        !stats of HRU map
        integer                         :: ncols_HRU
        integer                         :: nrows_HRU

        !parameter maps at the finest scale, created using transfer functions
        double precision, allocatable, dimension(:,:)   :: pm_map_SZM
        double precision, allocatable, dimension(:,:)   :: pm_map_LnTo
        double precision, allocatable, dimension(:,:)   :: pm_map_SRmax
        double precision, allocatable, dimension(:,:)   :: pm_map_SRinit
        double precision, allocatable, dimension(:,:)   :: pm_map_CHV
        double precision, allocatable, dimension(:,:)   :: pm_map_Td
        double precision, allocatable, dimension(:,:)   :: pm_map_Smax
        double precision, allocatable, dimension(:,:,:) :: pm_map_all_temp
        !YZ 2023
        double precision, allocatable, dimension(:,:)   :: pm_map_Drz
        double precision, allocatable, dimension(:,:)   :: pm_map_B
        double precision, allocatable, dimension(:,:)   :: pm_map_Ks

        !basin predictors
        real, allocatable, dimension(:,:)   :: sand_0_10
        real, allocatable, dimension(:,:)   :: silt_0_10
        real, allocatable, dimension(:,:)   :: clay_0_10
        real, allocatable, dimension(:,:)   :: orgm_0_10
        double precision, allocatable, dimension(:,:)   :: ksat_decline
        real, allocatable, dimension(:,:)   :: srmax_lcm
        real, allocatable, dimension(:,:)   :: smax_d2r


        !variables used when writing out a file
        character(len=1024)                             :: temp_fn      !stores dynamic filename
        character(len = 5)                              :: i_n_str      !string storing current number of pm files made
        character(len = 6)                              :: fmt = '(I5.5)' !format descriptor - integer of width 5 with zeros on left
        integer                                         :: x,y,i


        !print *,
        !print *,'          ****************************************************************'
        print *, '             Apply pedo-transfer functions to create parameter map'
        print *, '          ****************************************************************'

        !print *,
        print *, 'n = ',i_n,' out of ',n_pm_maps

        !fill in local variables
        ncols_HRU = size(HRU_map(1,:))
        nrows_HRU = size(HRU_map(:,1))
        sand_0_10 = bp_maps(:,:,1)
        silt_0_10 = bp_maps(:,:,2)
        clay_0_10 = bp_maps(:,:,3)
        orgm_0_10 = bp_maps(:,:,4)
        !ksat_decline=bp_maps(:,:,5)
        srmax_lcm=bp_maps(:,:,12)
        smax_d2r=bp_maps(:,:,14)

        if (i_n == 1) then
            allocate(pm_map_all(nrows_HRU, ncols_HRU, 9))
            allocate(pm_map_all_temp(nrows_HRU, ncols_HRU, 9))
        end if

        !allocate(pm_map_Drz(nrows_HRU, ncols_HRU))

        !apply pedotransfer functions to map n to get parameter maps
        CALL pedotf_LnTo(i_n,pedo_tf_all(2), glob_pms_LnTo, sand_0_10, clay_0_10, pm_map_LnTo,bp_maps)

        IF (pedo_tf_all(1)==3) THEN !apply SZM with lnT0 parameters
            CALL pedotf_SZM(i_n, 2, glob_pms_LnTo,bp_maps, pm_map_SZM, soilmusiddata,folder_output,pm_map_LnTo)
        ELSE  !apply SZM ptf with SZM global parameters
            CALL pedotf_SZM(i_n, pedo_tf_all(1), glob_pms_SZM,bp_maps, pm_map_SZM, soilmusiddata,folder_output,pm_map_LnTo)
        END IF

        CALL pedotf_SRmax(i_n, pedo_tf_all(3), glob_pms_SRmax,srmax_lcm,bp_maps(:,:,17),pm_map_SRmax, pm_map_Drz)
        CALL pedotf_SRinit(i_n, pedo_tf_all(4), glob_pms_SRinit,sand_0_10, pm_map_SRinit)
        CALL pedotf_CHV(i_n, pedo_tf_all(5), glob_pms_CHV,sand_0_10,pm_map_CHV)
        CALL pedotf_Td(i_n, pedo_tf_all(6), glob_pms_Td,sand_0_10, pm_map_Td,pm_map_LnTo)
        CALL pedotf_Smax(i_n, pedo_tf_all(7), glob_pms_Smax,smax_d2r,bp_maps(:,:,17),bp_maps(:,:,18),pm_map_Smax)
        CALL pedotf_B(i_n, pedo_tf_all(8), glob_pms_B,soil_class, pm_map_B)
        CALL pedotf_Ks(i_n, pedo_tf_all(9), glob_pms_Ks,soil_class, pm_map_Ks,pm_map_Drz,pm_map_szm)

        pm_map_all(:,:,1) = pm_map_SZM
        pm_map_all(:,:,2) = pm_map_LnTo
        pm_map_all(:,:,3) = pm_map_SRmax
        pm_map_all(:,:,4) = pm_map_SRinit
        pm_map_all(:,:,5) = pm_map_CHV
        pm_map_all(:,:,6) = pm_map_Td
        pm_map_all(:,:,7) = pm_map_Smax
        pm_map_all(:,:,8) = pm_map_B
        pm_map_all(:,:,9) = pm_map_Ks
        pm_map_all_temp=pm_map_all

        !save parameter maps if required
        IF (save_pm_maps >= 1) THEN
        !print *,
        print *, 'Saving the following parameter maps:'
            write(i_n_str,fmt) (i_n+start_seed-1)

                        !write nodata values so that maps match HRU maps
            DO x = 1,ncols_HRU
                DO y = 1,nrows_HRU
                    IF (HRU_map(y,x) == stats_HRU(4)) THEN !if nodata in HRU map
                        pm_map_all_temp(y,x,:) = stats_HRU(4)
                    END IF
                    DO i = 1,7
                        IF (pm_map_all_temp(y,x,i)>=100000) THEN
                            pm_map_all_temp(y,x,i) = 100000
                        END IF
                    END DO
                END DO
            END DO

            DO i = 1,9

            !get filename for each different parameter
            IF (i==1) THEN
                temp_fn = trim(folder_output) // "param_map_szm" // trim(i_n_str)//".asc"
            ELSE IF(i==2) THEN
                temp_fn = trim(folder_output) // "param_map_lnto" // trim(i_n_str)//".asc"
            ELSE IF(i==3) THEN
                temp_fn = trim(folder_output) // "param_map_srmax" // trim(i_n_str)//".asc"
            ELSE IF(i==4) THEN
                temp_fn = trim(folder_output) // "param_map_srinit" // trim(i_n_str)//".asc"
            ELSE IF(i==5) THEN
                temp_fn = trim(folder_output) // "param_map_chv" // trim(i_n_str)//".asc"
            ELSE IF(i==6) THEN
                temp_fn = trim(folder_output) // "param_map_td" // trim(i_n_str)//".asc"
            ELSE IF(i==7) THEN
                temp_fn = trim(folder_output) // "param_map_smax" // trim(i_n_str)//".asc"
            ELSE IF(i==8) THEN
                temp_fn = trim(folder_output) // "param_map_B" // trim(i_n_str)//".asc"
            ELSE IF(i==9) THEN
                temp_fn = trim(folder_output) // "param_map_Ks" // trim(i_n_str)//".asc"
            END IF

            !only save parameter map if this parameter has had MPR applied
            IF (pedo_tf_all(i)>1) THEN
                call write_ascii_grid(temp_fn,pm_map_all_temp(:,:,i),ncols_HRU,nrows_HRU,&
                    stats_HRU(1),stats_HRU(2),stats_HRU(3),stats_HRU(4),6)
                print *,'   Saving parameter map: ',trim(temp_fn)
            END IF
            END DO
        END IF


    end subroutine mpr_pedotransfer




        ! **********************************************************************************************
        ! subroutine: 6. mpr_upscale: Upscale param map to get param value for each HRU
        ! **********************************************************************************************
        ! For run n in N parameter maps
        ! a. Loop through all HRUs (or zones with different parameter values)
        ! b. Extract all parameter values from parameter map covering selected HRU
        ! c. Upscale parameters by taking an average
    subroutine mpr_upscale(HRU_map, &
        nHRUs, &
        pm_map_all, &
        up_pms, &
        pedo_tf_all, &
        save_pm_maps, &
        stats_HRU, &
        i_n, &
        start_seed, &
        folder_output,&
        pm_map_Drz,&
        up_pms_Drz)

        use dyna_utility
        use mpr_SZM
        use mpr_LnTo
        use mpr_SRmaxGW
        use mpr_SRinit
        use mpr_CHV
        use mpr_Td
        use mpr_Smax
        use mpr_upscaling

        implicit none

        !!DECLARE DUMMY VARIABLES
        integer, intent(in)                 :: i_n
        double precision, allocatable, dimension(:,:),intent(in)     :: HRU_map
        integer,intent(in)                                           :: nHRUS            !number of HRUs
        double precision, allocatable, dimension(:,:,:),intent(in)   :: pm_map_all !mp_map_all(x,y,[SZM, LnTo, SRmax, SRinit, CHV, Td, Smax])
        double precision, allocatable, dimension(:,:),intent(out)    :: up_pms       !vector of upscaled params (nHRUs,8)
        integer, dimension(11) , intent(in)                           :: pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax
        integer                 :: save_pm_maps    !set to 1 if user wants to save high-res parameter maps
        double precision, dimension(4)  :: stats_HRU !variable comprising HRU map xll, yll, cellsize & nodata value
        integer                 :: start_seed      !starting seed for generation of global params
        character(len=1024)     :: folder_output   !output folder to store parameter input files - OUTPUT folder within MPR main folder
        
        double precision, allocatable, dimension(:,:),intent(out)    :: up_pms_Drz 
        double precision, allocatable, dimension(:,:)   :: pm_map_Drz


        !!DECLARE LOCAL VARIABLES
        integer                 :: i_HRU

        !variables required for parameter upscaling
        double precision, allocatable, dimension(:)     :: pms_SZM       !vector of all param values in a HRU
        double precision, allocatable, dimension(:)     :: pms_LnTo
        double precision, allocatable, dimension(:)     :: pms_SRmax
        double precision, allocatable, dimension(:)     :: pms_SRinit
        double precision, allocatable, dimension(:)     :: pms_CHV
        double precision, allocatable, dimension(:)     :: pms_Td
        double precision, allocatable, dimension(:)     :: pms_Smax
        
        !YZ Drz
        double precision, allocatable, dimension(:)     :: pms_Drz_temp
        double precision, allocatable, dimension(:)     :: pms_Drz
        !YZ 2023 B&Ks
        double precision, allocatable, dimension(:)     :: pms_B
        double precision, allocatable, dimension(:)     :: pms_Ks

        !real :: start, finish

        !allocate plotting variables
        double precision, allocatable, dimension(:,:,:) :: pm_map_all_temp
        character(len=1024)                             :: temp_fn      !stores dynamic filename
        character(len = 5)                              :: i_n_str      !string storing current number of pm files made
        character(len = 6)                              :: fmt = '(I5.5)' !format descriptor - integer of width 5 with zeros on left
        integer                                         :: x,y,i,nrows,ncols
        integer                                         :: positive_count

        !allocate variable to store upscaled parameters
        allocate(up_pms(nHRUs,10))
        allocate(up_pms_Drz(nHRUs,1))
        

                !print *,
        !print *, '          ****************************************************************'
        print *, '             Upscale parameters '
        print *, '          ****************************************************************'


        DO i_HRU = 1,nHRUs !loop through for each HRU

            IF (mod(i_HRU,250)==0) THEN
                print *, '     upscaling i_HRU = ',i_HRU,' out of ',nHRUs
            END IF
            !call cpu_time(start)

            !!  ----- E) FIND AND EXTRACT PARAMETER VALUES FOR EACH HRU ------------------------------
            !pms_LnTo is a vector of all LnTo parameter values in the HRU number i_HRU
             !!  ----- F) UPSCALE PARAMETERS BY TAKING AVERAGE FOR EACH HRU ---------------------------
            !up_pms will determine the formatting of the parameter file
            up_pms(i_HRU,1) = i_HRU                                 !keep log of HRU numbers

            if (pedo_tf_all(1) > 1) then !skip this if set as fixed or global parameter as very time consuming!
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,1), pms_SZM)
                CALL upscale_harmonic(pms_SZM, up_pms(i_HRU,2))   !upscaled parameter value for this HRU
            else
                up_pms(i_HRU,2) = pm_map_all(1,1,1)
            end if

            if (pedo_tf_all(2) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,2), pms_LnTo)
                !CALL upscale_harmonic(pms_LnTo, up_pms(i_HRU,3))
                CALL upscale_arithmetic(pms_LnTo, up_pms(i_HRU,3)) !normally harmonic, but values being both positive and negative is resulting in huge values which should be near 0.
            else
                up_pms(i_HRU,3) = pm_map_all(1,1,2)
            end if
            
            !SRmax 
            if (pedo_tf_all(3) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,3), pms_SRmax)
                CALL upscale_arithmetic(pms_SRmax, up_pms(i_HRU,4))
            else
                up_pms(i_HRU,4) = pm_map_all(1,1,3)
            end if


            !YZ 2024 Drz fix values also need to upscale
            !Drz
            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_Drz, pms_Drz_temp)
            !remove all missing value before upscaling
            pms_Drz=pms_Drz_temp
            positive_count=0
            do i = 1, size(pms_Drz_temp)
                if (pms_Drz_temp(i) > 0.0) then
                    positive_count = positive_count + 1
                    pms_Drz(positive_count) = pms_Drz_temp(i)
                end if
            end do
            pms_Drz=pms_Drz(1:positive_count)

            CALL upscale_arithmetic(pms_Drz, up_pms_Drz(i_HRU,1))


            if (pedo_tf_all(4) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,4), pms_SRinit)
                CALL upscale_arithmetic(pms_SRinit, up_pms(i_HRU,5))
            else
                up_pms(i_HRU,5) = pm_map_all(1,1,4)
            end if

            if (pedo_tf_all(5) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,5),  pms_CHV)
                CALL upscale_geometric(pms_CHV, up_pms(i_HRU,6))
            else
                up_pms(i_HRU,6) = pm_map_all(1,1,5)
            end if

            if (pedo_tf_all(6) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,6), pms_Td)
                CALL upscale_harmonic(pms_Td, up_pms(i_HRU,7))
            else
                up_pms(i_HRU,7) = pm_map_all(1,1,6)
            end if

            if (pedo_tf_all(7) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,7),  pms_Smax)
                CALL upscale_arithmetic(pms_Smax, up_pms(i_HRU,8))
            else
                up_pms(i_HRU,8) = pm_map_all(1,1,7)
            end if

            !YZ 2023
            !B
             if (pedo_tf_all(8) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,8),  pms_B)
                CALL upscale_arithmetic(pms_B, up_pms(i_HRU,9))
            else
                up_pms(i_HRU,9) = pm_map_all(1,1,8)
            end if

            !Ks
             if (pedo_tf_all(9) > 1) then
                CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_all(:,:,9),  pms_Ks)
                CALL upscale_arithmetic(pms_Ks, up_pms(i_HRU,10))
            else
                up_pms(i_HRU,10) = pm_map_all(1,1,9)
            end if


            !call cpu_time(finish)
            !print '("Time = ",f6.3," seconds.")',finish-start
        END DO


             !save parameter maps if required
        IF (save_pm_maps > 1) THEN
        !print *,
        print *, 'Saving the following parameter maps:'
            write(i_n_str,fmt) (i_n+start_seed-1)

            !match maps to HRU map
            nrows=size(HRU_map(:,1))
            ncols=size(HRU_map(1,:))

            allocate(pm_map_all_temp(nrows,ncols,10))

            DO i = 1,9
                DO x = 1,ncols
                    DO y = 1,nrows
                        IF (HRU_map(y,x) == -9999) THEN !if nodata in HRU map
                            pm_map_all_temp(y,x,:) = -9999
                        ELSEIF (HRU_map(y,x) == 0) THEN !nodata in HRU map
                            pm_map_all_temp(y,x,:) = 0
                        ELSE
                            pm_map_all_temp(y,x,i) = up_pms(int(HRU_map(y,x)),i+1)
                        END IF
                        IF (pm_map_all_temp(y,x,i)>=100000) THEN
                            pm_map_all_temp(y,x,i) = 100000
                        END IF
                    END DO
                END DO

            !get filename for each different parameter
            IF (i==1) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_szm" // trim(i_n_str)//".asc"
            ELSE IF(i==2) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_lnto" // trim(i_n_str)//".asc"
            ELSE IF(i==3) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_srmax" // trim(i_n_str)//".asc"
            ELSE IF(i==4) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_srinit" // trim(i_n_str)//".asc"
            ELSE IF(i==5) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_chv" // trim(i_n_str)//".asc"
            ELSE IF(i==6) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_td" // trim(i_n_str)//".asc"
            ELSE IF(i==7) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_smax" // trim(i_n_str)//".asc"
            ELSE IF(i==8) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_B" // trim(i_n_str)//".asc"
            ELSE IF(i==9) THEN
                temp_fn = trim(folder_output) // "upscaled_param_map_Ks" // trim(i_n_str)//".asc"    
            
            END IF

            !only save parameter map if this parameter has had MPR applied
            IF (pedo_tf_all(i)>1) THEN
                call write_ascii_grid(temp_fn,pm_map_all_temp(:,:,i),size(HRU_map(1,:)),size(HRU_map(:,1)),&
                    stats_HRU(1),stats_HRU(2),stats_HRU(3),stats_HRU(4),6)
                print *,'   Saving parameter map: ',trim(temp_fn)
            END IF
            END DO
        END IF


    print *, "upscaling done."
    end subroutine mpr_upscale





        ! **********************************************************************************************
        ! subroutine: 6. mpr_write_paramfile
        ! **********************************************************************************************
        ! Writes a DYMOND parameter file - using parameters created within MPR.
    subroutine mpr_write_paramfile(i_n, &
        start_seed, &
        folder_output, &
        nHRUs, &
        up_pms)

        use dyna_utility

        implicit none

        !DECLARE DUMMY VARIABLES
        integer, intent(in)                                         :: i_n
        integer, intent(in)                                         ::  start_seed
        character(len=1024), intent(in)                             :: folder_output
        integer, intent(in)                                         :: nHRUS
        double precision, allocatable, dimension(:,:), intent(in)   :: up_pms

        !DECLARE LOCAL VARIABLES

        !variables for saving Dynamic TOPMODEL parameter file
        character(len=1024)                             :: temp_fn      !stores dynamic filename
        character(len = 5)                              :: i_n_str      !string storing current number of pm files made
        character(len=64), dimension(8)                 :: headers      !parameter file headers
        character(len = 6)                              :: fmt = '(I5.5)' !format descriptor - integer of width 5 with zeros on left

        write(i_n_str,fmt) (i_n+start_seed-1)

        headers(1) = 'param_layer'
        headers(2) = 'szm_def'
        headers(3) = 'lnto_def'
        headers(4) = 'srmax_def'
        headers(5) = 'srinit_def'
        headers(6) = 'chv_def'
        headers(7) = 'td_def'
        headers(8) = 'smax_def'

        temp_fn = trim(folder_output) // '/param_file_' // trim(i_n_str) // '.txt'
        print *, 'Writing param file: ',trim(temp_fn)
        open(unit=24,file=temp_fn,status='REPLACE',action = 'WRITE')
        !write(24, FMT=*) 'hello'
        write(24, FMT=*) nHRUs,size(up_pms(1,:)),' !number of parameter layers(rows/hrus), number of param names (columns-1)'

        CALL write_numeric_list_append(headers , up_pms, 6,24)
        write(24, FMT=*) '1 0.8 0.000000001 KINEMATIC   SOLUTION    PARS'
        close(unit=24)

    end subroutine mpr_write_paramfile





        ! **********************************************************************************************
        ! subroutine: 7. mpr_perturb_parammap
        ! **********************************************************************************************
        ! Perturbs the parameter map - to produce a different parameter map for one parameter.
    subroutine mpr_perturb_parammap(param_num, &
                it, &
                mcpar_orig, &
                mcpar_new, &
                p_add_mult)

        use dyna_utility

        implicit none

        !DECLARE DUMMY VARIABLES
        integer, intent(in)             :: param_num        !1=SZM, 2=lnto, etc
        integer, intent(in)             :: it       !iteration - Perterbation number for this parameter.
        double precision, allocatable, dimension(:,:), intent(in)      :: mcpar_orig
        double precision, allocatable, dimension(:,:), intent(inout)   :: mcpar_new
        double precision, allocatable, dimension(:,:,:), intent(in)      :: p_add_mult !parammap_add_mult

        !DECLARE LOCAL VARIABLES
        integer :: i
        integer :: npms

!        print *, 'param_num = ',param_num
!        print *, 'it = ',it
!!        print *, 'mcpar_orig = ',mcpar_orig
!!        print *, 'mcpar_new = ',mcpar_new


        !if iteration is 0 then make no changes
        if (it == 0) THEN
            mcpar_new(:,param_num) = mcpar_orig(:,param_num)
        else
            npms = size(mcpar_orig(:,param_num))
            print *, 'npms = ',npms
            DO i = 1,npms
!                print *,
!                print *, 'i = ',i
!                print *, mcpar_orig(i,param_num)
!                print *, mcpar_new
!                print *, mcpar_new(i,param_num)
!                print *, p_add_mult(param_num,it,1)
                mcpar_new(i,param_num)=mcpar_orig(i,param_num)+p_add_mult(param_num,it,1)
                mcpar_new(i,param_num)=mcpar_new(i,param_num)*p_add_mult(param_num,it,2)
            END DO
        end if


    end subroutine


!
!    subroutine MPR
!
!    ! Rosie Lane - 22nd August 2017
!
!    ! To run: ./MPR.e -MPRfolder /home/rl1023/2017_07_20_MPR_setup/ -fpath_HRUmap /projects/The_Env_Virtual_observatory/DynaTOP_data/Severn_Rosie2/preprocess_55025_10km/55025_classarray.asc
!
!    !fpath_HRUmap = "/projects/The_Env_Virtual_observatory/DynaTOP_data/Severn_Rosie2/preprocess_54057_uniform/54057_classarray.asc"
!    !fpath_HRUmap = "/projects/The_Env_Virtual_observatory/DynaTOP_data/Severn_Rosie2/preprocess_55025_10km/55025_classarray.asc"
!
!    !using classarray to set parameters per catchment rather than per hru
!    !fpath_HRUmap = "/projects/The_Env_Virtual_observatory/DynaTOP_data/Severn_Rosie2/preprocess_54057_uniform/54057_catch.asc"
!    !./MPR.e -MPRfolder /home/rl1023/2017_07_20_MPR_setup/ -fpath_HRUmap /projects/The_Env_Virtual_observatory/DynaTOP_data/Severn_Rosie2/preprocess_54057_uniform/54057_catch.asc
!
!    !use dta_utility
!    !use dta_param_upscaling
!    !use dta_MPR
!    use mpr_extract_BasPred
!    use mpr_control
!    use mpr_errors
!    use mpr_upscaling
!    use mpr_utility
!
!    !MPR settings and functions for each DynaTOP parameter
!    use mpr_SZM
!    use mpr_LnTo
!    use mpr_SRmax
!    use mpr_SRinit
!    use mpr_CHV
!    use mpr_Td
!    use mpr_Smax
!
!    implicit none
!
!    !key filepaths and folder locations
!    character(len=1024)     :: MPRfolder        !location of MPR folder containing INPUT, SETTINGS and OUTPUT folders
!    character(len=1024)     :: fpath_HRUmap     !filepath pointing to HRU .asc file created in DTA
!    character(len=1024)     :: folder_input   !folder containing Basin Predictor files - INPUT folder within MPR main folder
!    character(len=1024)     :: folder_output   !output folder to store parameter input files - OUTPUT folder within MPR main folder
!    character(len=1024)     :: fname_control    !full path to control file
!
!    !increments for do statements etc.
!    integer                 :: i_HRU
!    integer                 :: i_n
!    double precision        :: zero=0
!
!    !variables in the control file
!    integer                 ::  n_pm_maps       !number of parameter files to create
!    integer                 ::  start_seed      !starting seed for generation of global params
!    integer                 ::  save_pm_maps    !set to 1 if user wants to save parameter maps
!    integer                 ::  save_bp_maps    !set to 1 if user wants to save basin predictor maps
!    integer                 ::  pedo_tf_SZM
!    integer                 ::  pedo_tf_LnTo
!    integer                 ::  pedo_tf_SRmax   !pedo-transfer equation numbers for all params
!    integer                 ::  pedo_tf_SRinit
!    integer                 ::  pedo_tf_CHV
!    integer                 ::  pedo_tf_Td
!    integer                 ::  pedo_tf_Smax
!
!    !parameter min, default and max values from parameter file
!    double precision, dimension(3)  :: SZM_range
!    double precision, dimension(3)  :: LnTo_range
!    double precision, dimension(3)  :: SRmax_range
!    double precision, dimension(3)  :: SRinit_range
!    double precision, dimension(3)  :: CHV_range
!    double precision, dimension(3)  :: Td_range
!    double precision, dimension(3)  :: Smax_range
!
!    !variables combining results from the control file
!    integer, dimension(7)   ::  pedo_tf_all     !pedo-transfer eq. nums for SZM, LnTo, SRmax, SRinit, CHV, Td, Smax
!
!    !variables defining HRU map and parameter maps
!    double precision, allocatable, dimension(:,:) :: HRU_map
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
!    integer                         :: nHRUS            !number of HRUs
!    double precision, dimension(4)  :: stats_HRU        !variable comprising HRU map xll, yll, cellsize & nodata value
!
!    !basin predictor maps needed - logical arrays with true values if maps required
!    logical, dimension(4,5) :: req_soils    !Sand, silt, clay and OM rows, split by depth columns
!
!    !basin predictor maps
!    double precision, allocatable, dimension(:,:,:,:)    :: bp_maps !bp maps (x,y,sand/silt/clay/orgm,depth_class)
!    double precision, allocatable, dimension(:,:)   :: sand_0_10
!    double precision, allocatable, dimension(:,:)   :: silt_0_10
!    double precision, allocatable, dimension(:,:)   :: clay_0_10
!    double precision, allocatable, dimension(:,:)   :: orgm_0_10
!
!    !global parameters
!    integer, dimension(7)                           :: n_gp_all     !vector of number of global params for all params
!    double precision, allocatable, dimension(:,:)   :: glob_pms_SZM     !dim(glob pm number, 1:n_pm_maps)
!    double precision, allocatable, dimension(:,:)   :: glob_pms_LnTo
!    double precision, allocatable, dimension(:,:)   :: glob_pms_SRmax
!    double precision, allocatable, dimension(:,:)   :: glob_pms_SRinit
!    double precision, allocatable, dimension(:,:)   :: glob_pms_CHV
!    double precision, allocatable, dimension(:,:)   :: glob_pms_Td
!    double precision, allocatable, dimension(:,:)   :: glob_pms_Smax
!
!    !parameter maps at the finest scale, created using transfer functions
!    double precision, allocatable, dimension(:,:)   :: pm_map_SZM
!    double precision, allocatable, dimension(:,:)   :: pm_map_LnTo
!    double precision, allocatable, dimension(:,:)   :: pm_map_SRmax
!    double precision, allocatable, dimension(:,:)   :: pm_map_SRinit
!    double precision, allocatable, dimension(:,:)   :: pm_map_CHV
!    double precision, allocatable, dimension(:,:)   :: pm_map_Td
!    double precision, allocatable, dimension(:,:)   :: pm_map_Smax
!
!    !variables required for parameter upscaling
!    double precision, allocatable, dimension(:)     :: pms_SZM       !vector of all param values in a HRU
!    double precision, allocatable, dimension(:)     :: pms_LnTo
!    double precision, allocatable, dimension(:)     :: pms_SRmax
!    double precision, allocatable, dimension(:)     :: pms_SRinit
!    double precision, allocatable, dimension(:)     :: pms_CHV
!    double precision, allocatable, dimension(:)     :: pms_Td
!    double precision, allocatable, dimension(:)     :: pms_Smax
!    double precision, allocatable, dimension(:,:)   :: up_pms       !vector of upscaled params
!
!    !variables for saving Dynamic TOPMODEL parameter file
!    character(len=1024)                             :: temp_fn      !stores dynamic filename
!    character(len = 5)                              :: i_n_str      !string storing current number of pm files made
!    character(len=64), dimension(8)                 :: headers      !parameter file headers
!    character(len = 6)                              :: fmt = '(I5.5)' !format descriptor - integer of width 5 with zeros on left
!    integer                                         :: x,y
!
!
!
!    !   ******************************************************************************************
!    !       0. GET LOCATION OF CONTROL FILE AND OTHER FILEPATHS FROM COMMAND LINE - ERROR CHECK
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '0. Get filepaths from command line'
!    print *, '**************************************************************************'
!
!    call read_mpr_commands(MPRfolder, fpath_HRUmap, folder_input, folder_output, fname_control)
!    call check_files_exist(MPRfolder, fpath_HRUmap, folder_input, folder_output, fname_control)
!
!
!    !   ******************************************************************************************
!    !       1. READ PARAMETER AND CONTROL FILE
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '1. read parameter and control file'
!    print *, '**************************************************************************'
!
!    CALL read_param_file(MPRfolder,SZM_range,LnTo_range,SRmax_range,SRinit_range,CHV_range,Td_range,Smax_range)
!
!    CALL read_control_file_int(fname_control, "n_pm_maps", n_pm_maps)
!    CALL read_control_file_int(fname_control, "start_seed", start_seed)
!    CALL read_control_file_int(fname_control, "save_pm_maps", save_pm_maps)
!    CALL read_control_file_int(fname_control, "save_bp_maps", save_bp_maps)
!    CALL read_control_file_int(fname_control, "pedo_tf_SZM", pedo_tf_SZM)
!    CALL read_control_file_int(fname_control, "pedo_tf_LnTo", pedo_tf_LnTo)
!    CALL read_control_file_int(fname_control, "pedo_tf_SRmax", pedo_tf_SRmax)
!    CALL read_control_file_int(fname_control, "pedo_tf_SRinit", pedo_tf_SRinit)
!    CALL read_control_file_int(fname_control, "pedo_tf_CHV", pedo_tf_CHV)
!    CALL read_control_file_int(fname_control, "pedo_tf_Td", pedo_tf_Td)
!    CALL read_control_file_int(fname_control, "pedo_tf_Smax", pedo_tf_Smax)
!
!    pedo_tf_all(1:4) = (/pedo_tf_SZM, pedo_tf_LnTo, pedo_tf_SRmax, pedo_tf_SRinit /)
!    pedo_tf_all(5:7) = (/pedo_tf_CHV, pedo_tf_Td, pedo_tf_Smax /)
!
!    !   ******************************************************************************************
!    !       2. READ HRU MAP
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '2. Read HRU map'
!    print *, '**************************************************************************'
!
!    print *, 'Reading map from: ',trim(fpath_HRUmap)
!
!    call read_ascii_grid(fpath_HRUmap, HRU_map, ncols_HRU, nrows_HRU, xllcorner_HRU, yllcorner_HRU, cellsize_HRU, nodata_HRU)
!
!    nHRUs = maxval(HRU_map)
!    stats_HRU(1:4) = (/ xllcorner_HRU, yllcorner_HRU, cellsize_HRU, nodata_HRU/)
!
!    !assume stats of param maps are the same as HRU maps
!    nrow = nrows_HRU
!    ncol = ncols_HRU
!    xll = xllcorner_HRU
!    yll = yllcorner_HRU
!    cellsize = cellsize_HRU
!
!    print *, 'Number of HRUs found: ',nHRUs
!    print *, 'Rows in HRU file: ',nrows_HRU
!    print *, 'Cols in HRU file: ',ncols_HRU
!
!    !convert nodata value to -9999 to avoid confusion
!    IF (stats_HRU(4) /= -9999) THEN
!        stats_HRU(4) = -9999
!        print *,'old nodata value: ',nodata_HRU,' being converted to ',stats_HRU(4)
!        call set_nodata_value( HRU_map, nodata_HRU, stats_HRU(4))
!        nodata_HRU = -9999
!    ELSE
!       print *, 'nodata value: ', nodata_HRU
!    END IF
!
!    !convert all 0 values in HRU map to nodata as well
!    !call set_nodata_value(HRU_map,zero,nodata_HRU)
!
!
!    !   ******************************************************************************************
!    !       3. READ REQUIRED BASIN PREDICTOR MAPS AND TRIM TO HRU MAP EXTENT
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '3. Read required basin predictor maps'
!    print *, '**************************************************************************'
!    print *, 'Reading maps from: ',trim(folder_input)
!
!    ! assume no basin predictors are required unless specified otherwise
!    req_soils(1, 1:5) = (/ .false. , .false., .false., .false., .false. /)
!    req_soils(2, 1:5) = (/ .false. , .false., .false., .false., .false. /)
!    req_soils(3, 1:5) = (/ .false. , .false., .false., .false., .false. /)
!    req_soils(4, 1:5) = (/ .false. , .false., .false., .false., .false. /)
!
!    ! Check which basin predictors are required
!    ! get_BasPred_req returns a logical list of the required soils maps req_soils(type,depth),
!    ! which contains .true. where soils maps are required
!    CALL get_BasPred_SZM(pedo_tf_SZM, req_soils)
!    CALL get_BasPred_LnTo(pedo_tf_LnTo, req_soils)
!    CALL get_BasPred_SRmax(pedo_tf_SRmax, req_soils)
!    CALL get_BasPred_SRinit(pedo_tf_SRinit, req_soils)
!    CALL get_BasPred_CHV(pedo_tf_CHV, req_soils)
!    CALL get_BasPred_Td(pedo_tf_Td, req_soils)
!    CALL get_BasPred_Smax(pedo_tf_Smax, req_soils)
!
!
!    ! read_req_soils reads in any basin predictor maps required for user-selected pedo transfer functions
!    CALL read_req_soils(MPRfolder,req_soils,stats_HRU,HRU_map,save_bp_maps,bp_maps)
!    sand_0_10 = bp_maps(:,:,1)
!    silt_0_10 = bp_maps(:,:,2)
!    clay_0_10 = bp_maps(:,:,3)
!    orgm_0_10 = bp_maps(:,:,4)
!
!
!    !   ******************************************************************************************
!    !       4. DEFINE GLOBAL PARAMETERS WITHIN GIVEN BOUNDS AND GENERATE N SETS
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '4. Initialise global parameters within given bounds'
!    print *, '**************************************************************************'
!
!    n_gp_all = (/ 0,0,0,0,0,0,0 /)  !initialise vector storing the number of global parameters needed for each parameter.
!    !call init_gp_ functions to fill in the glob_pms_ arrays, glob_pms(pm number, 1:n)
!    CALL init_gp_SZM(pedo_tf_SZM,glob_pms_SZM, n_pm_maps, start_seed, n_gp_all, SZM_range)
!    CALL init_gp_LnTo(pedo_tf_LnTo,glob_pms_LnTo, n_pm_maps, start_seed, n_gp_all, LnTo_range)
!    CALL init_gp_SRmax(pedo_tf_SRmax,glob_pms_SRmax, n_pm_maps, start_seed, n_gp_all, SRmax_range)
!    CALL init_gp_SRinit(pedo_tf_SRinit,glob_pms_SRinit, n_pm_maps, start_seed, n_gp_all, SRinit_range)
!    CALL init_gp_CHV(pedo_tf_CHV,glob_pms_CHV, n_pm_maps, start_seed, n_gp_all,CHV_range)
!    CALL init_gp_Td(pedo_tf_Td,glob_pms_Td, n_pm_maps, start_seed, n_gp_all,Td_range)
!    CALL init_gp_Smax(pedo_tf_Smax,glob_pms_Smax, n_pm_maps, start_seed, n_gp_all,Smax_range)
!
!
!    !   ******************************************************************************************
!    !       5. APPLY PEDO-TRANSFER FUNCTIONS TO CREATE N PARAMETER MAPS
!    !   ******************************************************************************************
!    print *,
!    print *, '**************************************************************************'
!    print *, '5. Apply pedo-transfer functions to create n parameter maps'
!    print *, '**************************************************************************'
!
!    !allocate variable to store upscaled parameters
!    allocate(up_pms(nHRUs,8))
!
!    DO i_n = 1,n_pm_maps
!
!        print *,
!        print *, 'n = ',i_n,' out of ',n_pm_maps
!
!        !apply pedotransfer functions to map i to get parameter maps
!        CALL pedotf_SZM(i_n, pedo_tf_SZM, glob_pms_SZM,sand_0_10, pm_map_SZM)
!        CALL pedotf_LnTo(i_n,pedo_tf_LnTo, glob_pms_LnTo, sand_0_10, clay_0_10, pm_map_LnTo)
!        CALL pedotf_SRmax(i_n, pedo_tf_SRmax, glob_pms_SRmax,sand_0_10, pm_map_SRmax)
!        CALL pedotf_SRinit(i_n, pedo_tf_SRinit, glob_pms_SRinit,sand_0_10, pm_map_SRinit)
!        CALL pedotf_CHV(i_n, pedo_tf_CHV, glob_pms_CHV,sand_0_10, pm_map_CHV)
!        CALL pedotf_Td(i_n, pedo_tf_Td, glob_pms_Td,sand_0_10, pm_map_Td)
!        CALL pedotf_Smax(i_n, pedo_tf_Smax, glob_pms_Smax,sand_0_10, pm_map_Smax)
!
!        !save parameter maps if required
!        if (save_pm_maps == 1) then
!            write(i_n_str,fmt) (i_n+start_seed-1)
!            temp_fn = trim(folder_output) // "/param_map_LnTo" // trim(i_n_str)//".asc"
!            !write nodata values so that maps match HRU maps
!            DO x = 1,ncols_HRU
!                DO y = 1,nrows_HRU
!                    IF (HRU_map(y,x) == stats_HRU(4)) THEN !if nodata in HRU map
!                        pm_map_LnTo(y,x) = stats_HRU(4)
!                    END IF
!                END DO
!            END DO
!            call write_ascii_grid(temp_fn,pm_map_LnTo,ncol,nrow,xll,yll,cellsize,nodata_HRU,6)
!            print *, '     Saving parameter map in: ',trim(temp_fn)
!        end if
!
!
!        !   ******************************************************************************************
!        !       6. UPSCALE - 1 PARAMETER FOR EACH HRU
!        !   ******************************************************************************************
!
!    !allocate variable to store upscaled parameters
!    allocate(up_pms(nHRUs,8))
!
!        DO i_HRU = 1,nHRUs !loop through for each HRU
!
!            IF (mod(i_HRU,25)==0) THEN
!                print *, '     i_HRU = ',i_HRU,' out of ',nHRUs
!            END IF
!            !!  ----- E) FIND AND EXTRACT PARAMETER VALUES FOR EACH HRU ------------------------------
!            !pms_LnTo is a vector of all LnTo parameter values in the HRU number i_HRU
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_SZM, pms_SZM)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_LnTo, pms_LnTo)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_SRmax, pms_SRmax)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_SRinit, pms_SRinit)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_CHV,  pms_CHV)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_Td, pms_Td)
!            CALL extract_pm_HRU(i_HRU,HRU_map, pm_map_Smax,  pms_Smax)
!
!            !!  ----- F) UPSCALE PARAMETERS BY TAKING AVERAGE FOR EACH HRU ---------------------------
!            !up_pms will determine the formatting of the parameter file
!            up_pms(i_HRU,1) = i_HRU                                 !keep log of HRU numbers
!            CALL upscale_arithmetic(pms_SZM, up_pms(i_HRU,2))       !upscaled parameter values
!            CALL upscale_harmonic(pms_LnTo, up_pms(i_HRU,3))        !stored beside HRU number
!            CALL upscale_arithmetic(pms_SRmax, up_pms(i_HRU,4))
!            CALL upscale_arithmetic(pms_SRinit, up_pms(i_HRU,5))
!            CALL upscale_arithmetic(pms_CHV, up_pms(i_HRU,6))
!            CALL upscale_arithmetic(pms_Td, up_pms(i_HRU,7))
!            CALL upscale_arithmetic(pms_Smax, up_pms(i_HRU,8))
!
!
!        END DO
!
!        !   ******************************************************************************************
!        !       7. WRITE OUTPUT PARAMETER FILE
!        !   ******************************************************************************************
!
!        write(i_n_str,fmt) (i_n+start_seed-1)
!
!        headers(1) = 'param_layer'
!        headers(2) = 'szm_def'
!        headers(3) = 'lnto_def'
!        headers(4) = 'srmax_def'
!        headers(5) = 'srinit_def'
!        headers(6) = 'chv_def'
!        headers(7) = 'td_def'
!        headers(8) = 'smax_def'
!
!        temp_fn = trim(folder_output) // '/param_file_' // trim(i_n_str) // '.txt'
!        print *, 'Writing param file: ',trim(temp_fn)
!        open(unit=24,file=temp_fn,status='REPLACE',action = 'WRITE')
!        !write(24, FMT=*) 'hello'
!        write(24, FMT=*) nHRUs,size(up_pms(1,:)),' !number of parameter layers(rows/hrus), number of param names (columns-1)'
!
!        CALL write_numeric_list_append(headers , up_pms, 6,24)
!        write(24, FMT=*) '1 0.8 0.000000001 KINEMATIC   SOLUTION    PARS'
!        close(unit=24)
!    END DO
!
!
!
!
!
!end subroutine MPR




end module mpr_main


