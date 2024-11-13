!! Module containing subroutines for MPR.f90
!! This module contains all subroutines involving MPR for the SZM parameter:
!
!   get_BasPred_SZM    Complete logical list of basin predictors required for LnTo
!   init_gp_SZM        Initialise global parameters for LnTo
!   pedotf_SZM         Pedo-transfer function for LnTo

!! Rosie Lane - 31st August 2017

module mpr_SZM
contains


    ! **********************************************************************************************
    ! subroutine: get_BasPred_SZM : Complete logical list of basin predictors required
    ! **********************************************************************************************
    ! Given the user-specified pedo transfer function selected for SZM, this subroutine adds to a
    ! logical list of required basin predictors, by specifying which predictors it requires.

    subroutine get_BasPred_SZM(pedo_tf_SZM, req_soils)

        implicit none

        ! declare dummy arguments
        integer, intent(inout)                  :: pedo_tf_SZM !pedo_tf_equation to be used
        logical, dimension(*), intent(inout)  :: req_soils    !soils data required logical

        ! define required BasinPredictors for LnTo
        select case(pedo_tf_SZM)

            case (0)    !fixed parameter - no basin predictors required

            case (1)    !global parameter - no basin predictors required

            case (2)     !Use a pedo-transfer equation

                req_soils(5) = .true.   !musid map
                req_soils(7:11) = (/.true.,.true.,.true.,.true.,.true./) !Using tables of exponential decline profiles from soils data.

                !using form of exponential decline in ksat from soils data
                !req_soils(5,1) = (.true.)
                !print *, "ERROR: A transfer function has not yet been programmed for SZM"
                !print *, "Select pedo_tf = 0 to treat it as a fixed parameter"
                !print *, "Select pedo_tf = 1 to treat it as a global parameter"
                !stop
            case (3)    !Same ptf equation as case 2, but using same global parameters as lnto

                req_soils(5) = .true.   !musid map
                req_soils(7:11) = (/.true.,.true.,.true.,.true.,.true./) !Using tables of exponential decline profiles from soils data.

            case (4)

                req_soils(5) = .true.   !musid map
                req_soils(7:11) = (/.true.,.true.,.true.,.true.,.true./) !Using tables of exponential decline profiles from soils data.

           case (5)

                req_soils(5) = .true.   !musid map
                req_soils(7:11) = (/.true.,.true.,.true.,.true.,.true./) !Using tables of exponential decline profiles from soils data.


            case default
                print *, "WARNING: a valid pedo-tf equation for SZM must be specified in the control file"
                print *, "SZM will be set to the default of fixed parameter"
                print *, "The following options can be selected in the control file:"
                print *, "pedo_tf_SZM = 0, sets it as a fixed parameter"
                print *, "pedo_tf_SZM = 1, sets it as a global parameter"
                !print *,

                pedo_tf_SZM = 0

        end select

    end subroutine get_BasPred_SZM





    ! **********************************************************************************************
    ! subroutine: init_gp_SZM : Initialise global parameters for SZM
    ! **********************************************************************************************
    ! This subroutine does the following:
    !   1. Defines how many global parameters (n_glob_pms) are needed for the user-specified transfer function
    !       TF 0    : 1 global parameter
    !       TF 1    : 3 global parameters
    !   2. Defines min/max ranges for global parameters
    !   3. Generates list of global parameters using the rand function
    !       - user specified start_seed can be used to skip to any point in this list
    !       - global parameter list is of dimensions number of global params by number of param files needed
    !   4. Transforms global parameters to be within set min/max ranges

    subroutine init_gp_SZM(pedo_tf_SZM,glob_pms_SZM, n_pm_maps, start_seed, n_gp_all,pm_range)

        use dta_utility
        use mpr_LnTo

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_SZM
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_SZM      !global param list
        integer, intent(in)                                         :: n_pm_maps
        integer, intent(in)                                         :: start_seed   !to start rand num generator
        integer, dimension(11), intent(inout)                        :: n_gp_all     !vector of number of global params for all params
        double precision, dimension(3), intent(in)                  :: pm_range    !min, fixed, max values for this parameter

        !declare local variables
        integer                                        :: n_glob_pms   !number of global parameters
        double precision, allocatable, dimension(:)    :: min_gp    !min values for global parameters
        double precision, allocatable, dimension(:)    :: max_gp
        real, dimension(:,:),allocatable               :: rand_nums
        double precision                               :: num
        integer, dimension(12)                         :: seed
        integer                                        :: i,j, rn_i


        !1. Define how many global parameters should exist and min/max ranges, based on selected transfer function
        select case(pedo_tf_SZM)

            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "SZM is a fixed parameter, of value ",pm_range(2)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(2)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(2)

            case (1)    !Treat as global parameter

                print *, "SZM is a global parameter, in the range ",pm_range(1)," ",pm_range(3)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(1)    !min and max values set from parameter file
                max_gp(1) = pm_range(3)

            case (2) !pedo-transfer function based on soil ksat declines

                print *, "SZM will be parameterised using ksat declines from soils data"
                !n_glob_pms = 2
                !allocate(min_gp(2))
                !allocate(max_gp(2))
                !min_gp(1)=-0.01
                !max_gp(1)=0.05
                !min_gp(2)=0.5
                !max_gp(2)=1.5

                !pedotransfer function used will be the same as lnto transfer function 2
                n_glob_pms = 7         !a_const, a_sand, a_clay,
                allocate(min_gp(7))
                allocate(max_gp(7))
                min_gp(1) = -4.5 !-3.5 - changed following dotty plots of performance showing lower values do better
                min_gp(2)= 0.006
                min_gp(3) = -0.02
                max_gp(1) = -1!0.3 - changed following dotty plots of performance showing lower values do better
                max_gp(2)= 0.03
                max_gp(3) = -0.0032

            !HYPRES limits Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                min_gp(4) = 3
                max_gp(4) = 11

                min_gp(5) = -0.5
                max_gp(5) = -1.3

                min_gp(6) = -0.03
                max_gp(6) = -0.1

                min_gp(7) = -0.8
                max_gp(7) = -0.24

            case (3) !pedo-transfer function based on soils ksat declines - using same global parameters as SZM

                n_glob_pms = 1
                allocate(min_gp(1))
                allocate(max_gp(1))
                max_gp(1) = 1
                min_gp(1) = 1

                print *, "SZM will be parameterised using ksat declines from soils data and the same parameters as lnT0"

           !case(4) redirects to case(2)

           case (5) !pedo-transfer function based on soil ksat declines

                print *, "SZM will be parameterised using ksat declines from soils data"
                !n_glob_pms = 2
                !allocate(min_gp(2))
                !allocate(max_gp(2))
                !min_gp(1)=-0.01
                !max_gp(1)=0.05
                !min_gp(2)=0.5
                !max_gp(2)=1.5

                !pedotransfer function used will be the same as lnto transfer function 2
                n_glob_pms = 8         !a_const, a_sand, a_clay,
                allocate(min_gp(8))
                allocate(max_gp(8))
                min_gp(1) = -4.5 !-3.5 - changed following dotty plots of performance showing lower values do better
                min_gp(2)= 0.006
                min_gp(3) = -0.02
                max_gp(1) = -1!0.3 - changed following dotty plots of performance showing lower values do better
                max_gp(2)= 0.03
                max_gp(3) = -0.0032

            !HYPRES limits Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
            !changed to match ln(T0) values from testing theoretical min and max ranges.
                min_gp(4) = 1 !3
                max_gp(4) = 6!11

                min_gp(5) = -1.3!-0.5
                max_gp(5) = -0.5!-1.3

                min_gp(6) = -0.1
                max_gp(6) = -0.003

                min_gp(7) = -0.4!-0.8
                max_gp(7) = 0!-0.24

                min_gp(8) = 1
                max_gp(8) = 4



        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(1) = n_glob_pms
        !print *, 'n_gp_all(1) = ',n_gp_all

        allocate(glob_pms_SZM(n_glob_pms, n_pm_maps))
        allocate(rand_nums(sum(n_gp_all), n_pm_maps+start_seed-1))

        !generate lists of global parameters between min and max ranges, of length n
        !need to have list such that starting seed = 50 would give same result as 50th number from seed=1
        Call RANDOM_SEED( GET = seed )
        seed(1:12) = (/ 3,3,3,3,3,3,3,3,3,3,3,3/)
        Call RANDOM_SEED( PUT = seed )
        call RANDOM_NUMBER ( rand_nums )
        !print *, rand_nums

        print *, 'SZM pedo-transfer has ',n_glob_pms,' global parameters'
        !print *, 'max values for SZM global param(s) : ', max_gp
        !print *, 'min values for SZM global param(s) : ', min_gp
        !print *, 'sum n_gp_all SZM = ',sum(n_gp_all)


        !now normalise random numbers to within bounds
        DO i = 1,n_glob_pms
            DO j = start_seed, start_seed+n_pm_maps-1

                rn_i = (sum(n_gp_all) - n_glob_pms )+ i                 !rn_i is index for the random number
                num = rand_nums(rn_i,j)                                 !the appropriate random number, in range 0-1
                num = (num * (max_gp(i)-min_gp(i)))+min_gp(i)           !normalise the random number to given min/max
                glob_pms_SZM(i, j-start_seed+1) = num

            END DO
            IF (n_pm_maps >= 2) THEN
                print *, 'First 2 values for parameter ',i,' : ',glob_pms_SZM(i,1:2)
            ELSE
                print *, 'Value for parameter ',i,' : ',glob_pms_SZM(i,1)
            END IF
        END DO


        !print *,


    end subroutine init_gp_SZM





    ! **********************************************************************************************
    ! subroutine: pedotf_SZM: Pedo-transfer function for SZM
    ! **********************************************************************************************
    ! This routine takes n_i, signalling that we are now calculating parameter map n_i out of n,
    ! It takes the list of global parameters and required basin predictors, and from that
    ! applies the pedo-transfer functions to produce a parameter map.

    subroutine pedotf_SZM(n_i, pedo_tf, glob_pms,bp, pm_map,soilmusiddata,folder_output,lnto_map)

        use dyna_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                          :: n_i          !param map number to do
        integer, intent(in)                                          :: pedo_tf      !equation selection
        double precision, dimension(:,:), allocatable, intent(in)    :: glob_pms     !list of global params for SZM
        real, dimension(:,:,:), allocatable, intent(in):: bp           !basin predictor maps
        double precision, dimension(:,:), allocatable, intent(out)   :: pm_map       !output parameter map
        double precision, allocatable, dimension(:,:,:), intent(in) :: soilmusiddata !musiddatatables.
        character(len=1024), intent(in)                             :: folder_output !MPR output folder
        double precision, dimension(:,:), allocatable, intent(in)   :: lnto_map

        !declare local variables
        integer             :: n_glob_pms       !number of global parameters
        integer             :: n_pm_maps        !number of parameter maps - defined by ncols in glob_pms_LnTo
        integer             :: i,j,d,e          !incrementers.
        integer             :: nrows_bp         !number of rows in basin predictor map
        integer             :: ncols_bp
        integer             :: nmusids          !number of musids in table.
        integer             :: nentries         !number of entries (series'/landuses) for this musid.
        integer             :: mid_i
        integer             :: row_start
        integer             :: row_end
        integer             :: nrows
        integer             :: len_include_data
        integer,dimension(5):: d_s      !start and end depths of depth classes
        integer,dimension(5):: d_e
        integer,allocatable,dimension(:,:)             :: nodata_found  !keep track of where nodata values are.
        integer             :: id_s, id_e
        character(len=1024) :: temp_filename
        character(len=2)    :: depth_string
        character(len=5)    :: ni_string
        double precision    :: a,b              !exponential fit constants
        double precision, allocatable, dimension(:) :: a_all, b_all
        double precision, allocatable, dimension(:,:,:) :: dep_data   !depth decline data
        double precision, allocatable, dimension(:,:,:) :: ptf_data   !depth decline data - ksat values
        double precision, allocatable, dimension(:,:)   :: fit_data   !depth decline data - ready for fitting
        double precision, allocatable, dimension(:)     :: weights          !weights for each series/LU within a musid
        double precision, allocatable, dimension(:)     :: musid_szm
        double precision, allocatable, dimension(:)     :: musid_area
        double precision                                :: area_data  !areas under ksat curve (integreated)

        !find number of rows and columns in basin predictor files - pm map should be same size
        nrows_bp = size(bp(:,1,1))
        ncols_bp = size(bp(1,:,1))
        n_pm_maps = size(glob_pms(1,:))
        !find number of rows and columns in basin predictor files - pm map should be same size
        n_pm_maps = size(glob_pms(1,:))
        allocate(pm_map(nrows_bp, ncols_bp))

        !vary transfer function, depending on user-inputted equation selection.
        select case(pedo_tf)

            case (0)    !Fixed parameter

                n_glob_pms = 1
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        pm_map(i,j) = glob_pms(1, n_i)

                    END DO
                END DO

            case (1)    !Global parameter

                n_glob_pms = 1
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        pm_map(i,j) = glob_pms(1, n_i)

                    END DO
                END DO

            case (2)    !Pedotransfer function using sand and clay.
                n_glob_pms = 2

                !initialise log files

                write(ni_string,"(I5.5)") n_i
                temp_filename = trim(folder_output) // "SZM_pedotransfer_" // ni_string //".log"
                !write(temp_filename, "(A16,I2,A4)") 'SZM_pedotransfer', n_i, '.log'
                OPEN(unit=99,file=trim(temp_filename), status='unknown')
                write(99,*) 'This file contains SZM ksat values, calculated using pedotransfer case 2.'
                write(99,*) 'DEPTH, MUSID, ENTRY, WEIGHT, KSAT'
                CLOSE(99)

                write(ni_string,"(I5.5)") n_i
                temp_filename = trim(folder_output) // "SZM_fits_" // ni_string //".log"
                !write (temp_filename, "(A9,I2,A4)") "SZM_fits_", n_i,".log"
                OPEN(unit = 99,file=trim(temp_filename),status='unknown')
                WRITE(99,*) 'MUSID_NUM,SERIES,ENTRY_NUM, a(1)/b(2), FIT VALUES'
                CLOSE(99)

                !loop through every MUSID
                nmusids = maxval(soilmusiddata(:,1,1))
                allocate(musid_szm(nmusids))
                print *,maxval(soilmusiddata(:,1,1))
                row_start = 1

                DO i = 1,nmusids

                    !allocate arrays for number of entries in this musid
                    row_end = row_start + soilmusiddata(row_start,17,1) -1
                    nrows = (row_end - row_start) +1
                    allocate(dep_data(nrows,17,5)) !dep_data (n_entries, n_info_columns, n_depth_profiles)
                    allocate(ptf_data(nrows,4,5))  !ptf_data (n_entries, MUSID/ENTRY_NUM/WEIGHTING/KSAT, n_depths)
                    allocate(a_all(nrows))
                    allocate(b_all(nrows))
                    allocate(weights(nrows))
                    allocate(nodata_found(nrows,5))

                    !!Extract all depth declines for this musid - looping through depths
                    DO d = 1,5
                        dep_data(1:nrows,1:17,d) = soilmusiddata(row_start:row_end,1:17,d)
                    END DO

                    !apply pedotransfer function to get ksat values
                    ptf_data(1:nrows,1:3,1:5) = dep_data(1:nrows,1:3,1:5)
                    DO e = 1,nrows
                        DO d = 1,5
                            !IF (dep_data(e,16,d)==0) THEN !apply Cosby et al. (1984) for mineral soils
                            IF (dep_data(e,11,d)<35) THEN !apply Cosby et al. (1984) for minearal soils - organic content less than 35%
                                !ln_ksat = a + b(%sand) + c(%clay)
                                ptf_data(e,4,d)=glob_pms(1,n_i)+(glob_pms(2,n_i)*dep_data(e,8,d))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(3,n_i)*dep_data(e,10,d))
                                !ksat = exp(ln_ksat)
                                ptf_data(e,4,d)= 2.54*(10**(ptf_data(e,4,d)))

                                !keep record of nodata values
                                IF(int(dep_data(e,8,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF

                            ELSE    !apply HYPRES for peaty soils, Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                                ptf_data(e,4,d)=glob_pms(4,n_i)+(glob_pms(5,n_i)*dep_data(e,12,d)**2)
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(6,n_i)*dep_data(e,11,d)**(-1))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(7,n_i)*dep_data(e,12,d)*dep_data(e,11,d))
                                ptf_data(e,4,d) = EXP(ptf_data(e,4,d))

                                !keep record of nodata values
                                IF(int(dep_data(e,11,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF
                            END IF


                        END DO
                    END DO

                    !write ksat declines to a log file
                    write(ni_string,"(I5.5)") n_i
                    temp_filename = trim(folder_output) // "SZM_pedotransfer_" // ni_string //".log"
                    !write(temp_filename, "(A16,I2,A4)") 'SZM_pedotransfer_',n_i, '.log'
                    OPEN(unit=99,file=trim(temp_filename), position='append',status='old')
                    DO d = 1,5
                        DO e = 1, nrows
                            write(99,*) d,ptf_data(e,1:4,d)
                        END DO

                    END DO
                    CLOSE(unit=99)


                    !fit exponential to each ksat recession - to get m values.
                    DO e = 1,nrows
                        d_s = (/1,11,26,51,101/)
                        d_e = (/10,25,50,100,150/)
                        len_include_data = 0

                        !Work out which data to include, and how to allocate fit_data.
                        DO d = 1,5
                            IF (nodata_found(e,d)== 0) THEN
                                !len_include_data = len_include_data + (d_e(d)-d_s(d))+1 !if fitting depths 1 to 150
                                len_include_data = len_include_data + 1
                            END IF
                        END DO

                        !allocate fit data.
                        allocate(fit_data(len_include_data,3))

                        !extract fit data for each depth, and get integers of all the depths to include in fit calculation.
                        id_s = 1
                        DO d = 1,5
                            IF (nodata_found(e,d)== 0) THEN
                                fit_data(id_s,1) = ptf_data(e,4,d) !y = ksat values
                                fit_data(id_s,2) = -(d_e(id_s) + d_s(id_s))/2 !x = -depth
                                id_s = id_s + 1
                                !id_e = id_s + (d_e(d) - d_s(d))
                                !fit_data(id_s:id_e,1) =(/(j, j=d_s(d),d_e(d), 1)/)
                                !fit_data(id_s:id_e,2) = ptf_data(e,4,d)
                                !id_s = id_e + 1
                            END IF
                        END DO

                        !calculate the exponential fit! To the fit_data which excludes all nodata values.
                        CALL calc_exp_fit(fit_data, a, b)
                        a_all(e) = a
                        b_all(e) = b
                        deallocate(fit_data)

                    END DO

                    !save fit information in a log file - to check fits look reasonable in matlab
                    write(ni_string,"(I5.5)") n_i
                    temp_filename = trim(folder_output) // "SZM_fits_" // ni_string //".log"
                    !write (temp_filename, "(A9,I2,A4)") "SZM_fits_", n_i,".log"
                    OPEN(unit = 99,file=trim(temp_filename),position='append',status='old')
                    DO e = 1,nrows
                        WRITE(99,*) i,soilmusiddata(row_start+e-1,4,1),e,1,a_all(e)
                        WRITE(99,*) i,soilmusiddata(row_start+e-1,4,1),e,2,b_all(e)
                    END DO
                    CLOSE(unit=99)


                    !take a weighted average of all m values, to get a musid parameter value.
                    weights = soilmusiddata(row_start:row_end,3,1)!weight is column 3.
                    IF (sum(weights) > 1.001) THEN
                        print *, 'error - weights within musids exceed 1'
                    END IF
                    IF (sum(weights) < 0.999) THEN
                        print *, 'error - weights within musids dont add up to 1'
                    END IF
                    musid_szm(i) = sum(weights * b_all)

                    !deallocate/set variables for reading in next musid
                    deallocate(ptf_data)
                    deallocate(dep_data)
                    deallocate(a_all)
                    deallocate(b_all)
                    deallocate(weights)
                    deallocate(nodata_found)

                    row_start = row_end+1

                END DO

                !plot parameter values back in space - using musid ascii grid.
                !   where nodata, use average szm value across modelled region.
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp
                        e = INT(bp(i,j,6))
                        IF(e<-9000) THEN
                            pm_map(i,j) = (sum(musid_szm))/nmusids
                        ELSE
                            pm_map(i,j) = musid_szm(e)
                        END IF
                    END DO
                END DO

            case (4)    !Pedotransfer function using sand and clay -
                !n_glob_pms = 2

                !loop through every MUSID
                nmusids = maxval(soilmusiddata(:,1,1))
                allocate(musid_szm(nmusids))
                allocate(musid_area(nmusids))
                print *,maxval(soilmusiddata(:,1,1))
                row_start = 1

                DO i = 1,nmusids

                    !allocate arrays for number of entries in this musid
                    row_end = row_start + soilmusiddata(row_start,17,1) -1
                    nrows = (row_end - row_start) +1
                    allocate(dep_data(nrows,17,5)) !dep_data (n_entries, n_info_columns, n_depth_profiles)
                    allocate(ptf_data(nrows,4,5))  !ptf_data (n_entries, MUSID/ENTRY_NUM/WEIGHTING/KSAT, n_depths)
                    allocate(a_all(nrows))
                    allocate(b_all(nrows))
                    allocate(weights(nrows))
                    allocate(nodata_found(nrows,5))

                    !!Extract all depth declines for this musid - looping through depths
                    DO d = 1,5
                        dep_data(1:nrows,1:17,d) = soilmusiddata(row_start:row_end,1:17,d)
                    END DO

                    !apply pedotransfer function to get ksat values
                    ptf_data(1:nrows,1:3,1:5) = dep_data(1:nrows,1:3,1:5)
                    DO e = 1,nrows
                        DO d = 1,5
                            !IF (dep_data(e,16,d)==0) THEN !apply Cosby et al. (1984) for mineral soils
                            IF (dep_data(e,11,d)<35) THEN !apply Cosby et al. (1984) for minearal soils - organic content less than 35%
                                !ln_ksat = a + b(%sand) + c(%clay)
                                ptf_data(e,4,d)=glob_pms(1,n_i)+(glob_pms(2,n_i)*dep_data(e,8,d))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(3,n_i)*dep_data(e,10,d))
                                !ksat = exp(ln_ksat)
                                ptf_data(e,4,d)= EXP(ptf_data(e,4,d))

                                !keep record of nodata values
                                IF(int(dep_data(e,8,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF

                            ELSE    !apply HYPRES for peaty soils, Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                                ptf_data(e,4,d)=glob_pms(4,n_i)+(glob_pms(5,n_i)*dep_data(e,12,d)**2)
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(6,n_i)*dep_data(e,11,d)**(-1))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(7,n_i)*dep_data(e,12,d)*dep_data(e,11,d))
                                ptf_data(e,4,d) = EXP(ptf_data(e,4,d))

                                !keep record of nodata values
                                IF(int(dep_data(e,11,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF
                            END IF


                        END DO
                    END DO

                    !fit exponential to each ksat recession - to get m values.
                    area_data=0
                    weights = soilmusiddata(row_start:row_end,3,1)
                    DO e = 1,nrows
                        d_s = (/1,11,26,51,101/)
                        d_e = (/10,25,50,100,150/)
                        len_include_data = 0

                        !extract fit data for each depth, and get integers of all the depths to include in fit calculation.
                        id_s = 1

                        DO d = 1,5
                            IF (nodata_found(e,d)== 0) THEN
                                area_data = area_data + weights(e)*(ptf_data(e,4,d)*(d_e(id_s)-d_s(id_s)))
                                id_s = id_s + 1
                            END IF
                        END DO

                    END DO

                    !take a weighted average of all m values, to get a musid parameter value.
                    !weights = soilmusiddata(row_start:row_end,3,1)!weight is column 3.
                    IF (sum(weights) > 1.001) THEN
                        print *, 'error - weights within musids exceed 1'
                    END IF
                    IF (sum(weights) < 0.999) THEN
                        print *, 'error - weights within musids dont add up to 1'
                    END IF
                    musid_area(i) = area_data

                    !deallocate/set variables for reading in next musid
                    deallocate(ptf_data)
                    deallocate(dep_data)
                    deallocate(a_all)
                    deallocate(b_all)
                    deallocate(weights)
                    deallocate(nodata_found)

                    row_start = row_end+1

                END DO

                !plot parameter values back in space - using musid ascii grid.
                !   where nodata, use average szm value across modelled region.
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp
                        e = INT(bp(i,j,6))
                        IF(e<-9000) THEN
                            pm_map(i,j) = exp(lnto_map(i,j))/((sum(musid_area))/nmusids)
                        ELSE
                            pm_map(i,j) = exp(lnto_map(i,j))/musid_area(e)
                        END IF
                    END DO
                END DO



                case (5)    !Pedotransfer function using sand and clay -
                !n_glob_pms = 2

                !loop through every MUSID
                nmusids = maxval(soilmusiddata(:,1,1))
                allocate(musid_szm(nmusids))
                allocate(musid_area(nmusids))
                print *,maxval(soilmusiddata(:,1,1))
                row_start = 1

                DO i = 1,nmusids

                    !allocate arrays for number of entries in this musid
                    row_end = row_start + soilmusiddata(row_start,17,1) -1
                    nrows = (row_end - row_start) +1
                    allocate(dep_data(nrows,17,5)) !dep_data (n_entries, n_info_columns, n_depth_profiles)
                    allocate(ptf_data(nrows,4,5))  !ptf_data (n_entries, MUSID/ENTRY_NUM/WEIGHTING/KSAT, n_depths)
                    allocate(a_all(nrows))
                    allocate(b_all(nrows))
                    allocate(weights(nrows))
                    allocate(nodata_found(nrows,5))

                    !!Extract all depth declines for this musid - looping through depths
                    DO d = 1,5
                        dep_data(1:nrows,1:17,d) = soilmusiddata(row_start:row_end,1:17,d)
                    END DO

                    !apply pedotransfer function to get ksat values
                    ptf_data(1:nrows,1:3,1:5) = dep_data(1:nrows,1:3,1:5)
                    DO e = 1,nrows
                        DO d = 1,5
                            !IF (dep_data(e,16,d)==0) THEN !apply Cosby et al. (1984) for mineral soils
                            IF (dep_data(e,11,d)<35) THEN !apply Cosby et al. (1984) for minearal soils - organic content less than 35%
                                !ln_ksat = a + b(%sand) + c(%clay)
                                ptf_data(e,4,d)=glob_pms(1,n_i)+(glob_pms(2,n_i)*dep_data(e,8,d))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(3,n_i)*dep_data(e,10,d))
                                !ksat = exp(ln_ksat)
                                ptf_data(e,4,d)= EXP(ptf_data(e,4,d))

                                !keep record of nodata values
                                IF(int(dep_data(e,8,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF

                            ELSE    !apply HYPRES for peaty soils, Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                                ptf_data(e,4,d)=glob_pms(4,n_i)+(glob_pms(5,n_i)*dep_data(e,12,d)**2)
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(6,n_i)*dep_data(e,11,d)**(-1))
                                ptf_data(e,4,d)=ptf_data(e,4,d)+(glob_pms(7,n_i)*dep_data(e,12,d)*dep_data(e,11,d))
                                ptf_data(e,4,d) = EXP(ptf_data(e,4,d))

                                !keep record of nodata values
                                IF(int(dep_data(e,11,d))==-9999) THEN
                                    nodata_found(e,d) = 1
                                ELSE
                                    nodata_found(e,d) = 0
                                END IF
                            END IF


                        END DO
                    END DO

                    !fit exponential to each ksat recession - to get m values.
                    area_data=0
                    weights = soilmusiddata(row_start:row_end,3,1)
                    DO e = 1,nrows
                        d_s = (/1,11,26,51,101/)
                        d_e = (/10,25,50,100,150/)
                        len_include_data = 0

                        !extract fit data for each depth, and get integers of all the depths to include in fit calculation.
                        id_s = 1

                        DO d = 1,5
                            IF (nodata_found(e,d)== 0) THEN
                                area_data = area_data + weights(e)*(ptf_data(e,4,d)*(d_e(id_s)-d_s(id_s)))
                                id_s = id_s + 1
                            END IF
                        END DO

                    END DO

                    !take a weighted average of all m values, to get a musid parameter value.
                    !weights = soilmusiddata(row_start:row_end,3,1)!weight is column 3.
                    IF (sum(weights) > 1.001) THEN
                        print *, 'error - weights within musids exceed 1'
                    END IF
                    IF (sum(weights) < 0.999) THEN
                        print *, 'error - weights within musids dont add up to 1'
                    END IF
                    musid_area(i) = area_data

                    !deallocate/set variables for reading in next musid
                    deallocate(ptf_data)
                    deallocate(dep_data)
                    deallocate(a_all)
                    deallocate(b_all)
                    deallocate(weights)
                    deallocate(nodata_found)

                    row_start = row_end+1

                END DO

                !plot parameter values back in space - using musid ascii grid.
                !   where nodata, use average szm value across modelled region.
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp
                        e = INT(bp(i,j,6))
                        IF(e<-9000) THEN
                            pm_map(i,j) = exp(lnto_map(i,j))/((sum(musid_area))/nmusids)
                            pm_map(i,j) = pm_map(i,j) * glob_pms(8,n_i)
                        ELSE
                            pm_map(i,j) = exp(lnto_map(i,j))/musid_area(e)
                            pm_map(i,j) = pm_map(i,j) * glob_pms(8,n_i)
                        END IF
                    END DO
                END DO


            case default
                print *, "ERROR: invalid pedo-transfer setting."
                stop

        end select

        DO i = 1,nrows_bp
            DO j = 1,ncols_bp
                !add in upper and lower caps
                IF (pm_map(i,j)<0.0001) THEN
                    pm_map(i,j)=0.0001
                ELSEIF (pm_map(i,j)>15) THEN
                    pm_map(i,j)=15
                END IF
            END DO
        END DO

    end subroutine pedotf_SZM




end module mpr_SZM


