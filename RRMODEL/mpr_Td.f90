!! Module containing subroutines for MPR.f90
!! This module contains all subroutines involving MPR for the Td parameter:
!
!   get_BasPred_Td    Complete logical list of basin predictors required for LnTo
!   init_gp_Td        Initialise global parameters for LnTo
!   pedotf_Td         Pedo-transfer function for LnTo

!! Rosie Lane - 31st August 2017

module mpr_Td
contains


    ! **********************************************************************************************
    ! subroutine: get_BasPred_LnTo : Complete logical list of basin predictors required for LnTo
    ! **********************************************************************************************
    ! Given the user-specified pedo transfer function selected for LnTo, this subroutine adds to a
    ! logical list of required basin predictors, by specifying which predictors it requires.

    subroutine get_BasPred_Td(pedo_tf_Td, req_soils)

        implicit none

        ! declare dummy arguments
        integer, intent(inout)                     :: pedo_tf_Td !pedo_tf_equation to be used
        logical, dimension(*), intent(inout)  :: req_soils    !soils data required logical

        ! define required BasinPredictors for LnTo
        select case(pedo_tf_Td)

            case (0)    !fixed parameter - no basin predictors required

            case (1)    !global parameter - no basin predictors required

            case (2)      !multiplier on ln(T0)
                !print *, "LnTo transfer function selected is Cosby et al. (1984)"
                !req_soils(1) = .true.   !sand
                !req_soils(3) = .true.   !clay
                !req_soils(4) = .true.   !organic content
                !req_soils(15) = .true.  !bulk density
                !req_soils(16) = .true.  !is organic? Logical

            case default
                print *, "WARNING: a valid pedo-tf equation for Td must be specified in the control file"
                print *, "Td will be set to the default of fixed parameter"
                print *, "The following options can be selected in the control file:"
                print *, "pedo_tf_Td = 0, sets it as a fixed parameter"
                print *, "pedo_tf_Td = 1, sets it as a global parameter"
                !print *,

                pedo_tf_Td = 0

        end select

    end subroutine get_BasPred_Td





    ! **********************************************************************************************
    ! subroutine: init_gp_Td : Initialise global parameters for Td
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

    subroutine init_gp_Td(pedo_tf_Td,glob_pms_Td, n_pm_maps, start_seed, n_gp_all,pm_range)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_Td
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_Td      !global param list
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
        select case(pedo_tf_Td)

            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "Td is a fixed parameter, of value ",pm_range(2)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(2)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(2)

            case (1)    !Treat as global parameter

                print *, "Td is a global parameter, in the range ",pm_range(1)," ",pm_range(3)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(1)    !min and max values set from parameter file
                max_gp(1) = pm_range(3)

            case (2)   !multiplier on ln(T0)


                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = 1   !min and max values set from parameter file
                max_gp(1) = 10
                print *, "Td is a multiplier on ln(T0)."


            case default
                print *, "ERROR: pedo_tf_Td incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(6) = n_glob_pms

        allocate(glob_pms_Td(n_glob_pms, n_pm_maps))
        allocate(rand_nums(sum(n_gp_all), n_pm_maps+start_seed-1))

        !generate lists of global parameters between min and max ranges, of length n
        !need to have list such that starting seed = 50 would give same result as 50th number from seed=1
        Call RANDOM_SEED( GET = seed )
        seed(1:12) = (/ 3,3,3,3,3,3,3,3,3,3,3,3/)
        Call RANDOM_SEED( PUT = seed )
        call RANDOM_NUMBER ( rand_nums )



        !now normalise random numbers to within bounds

        DO i = 1,n_glob_pms
            DO j = start_seed, start_seed+n_pm_maps-1


                rn_i = (sum(n_gp_all) - n_glob_pms) + i                 !rn_i is index for the random number
                num = rand_nums(rn_i,j)                                 !the appropriate random number, in range 0-1
                num = (num * (max_gp(i)-min_gp(i)))+min_gp(i)           !normalise the random number to given min/max
                glob_pms_Td(i, j-start_seed+1) = num

            END DO
            IF (n_pm_maps >= 2) THEN
                print *, 'First 2 values for parameter ',i,' : ',glob_pms_Td(:,1:2)
            END IF
        END DO


        !print *,

    end subroutine init_gp_Td





    ! **********************************************************************************************
    ! subroutine: pedotf_LnTo: Pedo-transfer function for LnTo
    ! **********************************************************************************************
    ! This routine takes n_i, signalling that we are now calculating parameter map n_i out of n,
    ! It takes the list of global parameters and required basin predictors, and from that
    ! applies the pedo-transfer functions to produce a parameter map.

    subroutine pedotf_Td(n_i, pedo_tf, glob_pms,bp, pm_map,lnto_map)

        implicit none

        !declare dummy variables
        integer, intent(in)                                          :: n_i         !param map number to do
        integer, intent(in)                                          :: pedo_tf     !equation selection
        double precision, dimension(:,:), allocatable, intent(in)    :: glob_pms    !list of global params
        real, dimension(:,:), allocatable, intent(in)    :: bp          !basin predictor map (for sizing)
        double precision, dimension(:,:), allocatable, intent(out)   :: pm_map      !output parameter map
        double precision, dimension(:,:), allocatable, intent(in)   :: lnto_map      !map for lnto parameter

        !declare local variables
        integer             :: n_glob_pms       !number of global parameters
        integer             :: n_pm_maps        !number of parameter maps - defined by ncols in glob_pms_LnTo
        integer             :: i
        integer             :: j
        integer             :: nrows_bp         !number of rows in basin predictor map
        integer             :: ncols_bp


        !find number of rows and columns in basin predictor files - pm map should be same size
        nrows_bp = size(bp(:,1))
        ncols_bp = size(bp(1,:))
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

            case (2)   !multiplier on ln(T0): Td = a*exp(lnto)

                n_glob_pms = 1
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        !pm_map(i,j) = glob_pms(1, n_i)*(lnto_map(i,j)+20)
                        pm_map(i,j) = glob_pms(1,n_i)*EXP(lnto_map(i,j))

                    END DO
                END DO



            case default
                print *, "ERROR: invalid pedo-transfer setting."
                stop

        end select

         ! set a limit - to stop crazy values of 6x10^25 causing crashes!
        DO i = 1,nrows_bp
            DO j = 1,ncols_bp
                !add in upper and lower caps
                IF (pm_map(i,j)<0.01) THEN
                    pm_map(i,j)=0.01
                ELSEIF (pm_map(i,j)>200) THEN
                    pm_map(i,j)=200
                END IF
            END DO
        END DO


    end subroutine pedotf_Td




end module mpr_Td


