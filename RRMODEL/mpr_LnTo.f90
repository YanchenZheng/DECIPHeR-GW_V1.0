!! Module containing subroutines for MPR.f90
!! This module contains all subroutines involving MPR for the LnTo parameter:
!
!   get_BasPred_LnTo    Complete logical list of basin predictors required for LnTo
!   init_gp_LnTo        Initialise global parameters for LnTo
!   pedotf_LnTo         Pedo-transfer function for LnTo

!! Rosie Lane - 31st August 2017

module mpr_LnTo
contains


    ! **********************************************************************************************
    ! subroutine: get_BasPred_LnTo : Complete logical list of basin predictors required for LnTo
    ! **********************************************************************************************
    ! Given the user-specified pedo transfer function selected for LnTo, this subroutine adds to a
    ! logical list of required basin predictors, by specifying which predictors it requires.

    subroutine get_BasPred_LnTo(pedo_tf_LnTo, req_soils)

        implicit none

        ! declare dummy arguments
        integer, intent(inout)                  :: pedo_tf_LnTo !pedo_tf_equation to be used
        logical, dimension(*), intent(inout)  :: req_soils    !soils data required logical

        ! define required BasinPredictors for LnTo
        select case(pedo_tf_LnTo)

            case (0)    !fixed parameter - no basin predictors required

            case (1)    !global parameter - no basin predictors required

            case(2)     !Use pedo-transfer as Cosby et al. (1984)
                !print *, "LnTo transfer function selected is Cosby et al. (1984)"
                req_soils(1) = .true.   !sand
                req_soils(3) = .true.   !clay
                req_soils(4) = .true.   !organic content
                req_soils(15) = .true.  !bulk density
                req_soils(16) = .true.  !is organic? Logical

            case(3)     !Use pedo-transfer as Cosby et al. (1984), but separate areas with productive geology
                !print *, "LnTo transfer function selected is Cosby et al. (1984)"
                req_soils(1) = .true.   !sand
                req_soils(3) = .true.   !clay
                req_soils(4) = .true.   !organic content
                req_soils(15) = .true.  !bulk density
                req_soils(16) = .true.  !is organic? Logical
                req_soils(18) = .true.  !is highly productive geology?


            case(4)     !Third transfer function has not yet been defined
                print *, "LnTo pedo-transfer 2 has not yet been defined!"
               print *, "Select pedo_tf = 0 to treat it as a fixed parameter"
                print *, "Select pedo_tf = 1 to treat it as a global parameter"
                print *,  "Select pedo_tf = 2 to use the pedo-transfer function of Cosby et al. (1984)"
                stop

            case default
                print *, "WARNING: a valid pedo-tf equation for LnTo must be specified in the control file"
                print *, "LnTo will be set to the default of fixed parameter"
                print *, "pedo_tf_LnTo = 0, sets LnTo as a fixed parameter"
                print *, "pedo_tf_LnTo = 1, sets LnTo as a global parameter"
                print *, "pedo_tf_LnTo = 2, sets the LnTo pedo-transfer eq. to Cosby et al."

                pedo_tf_LnTo = 0

        end select

    end subroutine get_BasPred_LnTo





    ! **********************************************************************************************
    ! subroutine: init_gp_LnTo : Initialise global parameters for LnTo
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

    subroutine init_gp_LnTo(pedo_tf_LnTo,glob_pms_LnTo, n_pm_maps, start_seed, n_gp_all,pm_range)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_LnTo
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_LnTo      !global param list
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
        select case(pedo_tf_LnTo)

            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "LnTo is a fixed parameter, of value ",pm_range(2)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(2)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(2)

            case (1)    !Treat as global parameter

                print *, "LnTo is a global parameter, in the range ",pm_range(1)," ",pm_range(3)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(1)    !min and max values set from parameter file
                max_gp(1) = pm_range(3)

            case(2)     !Use pedo-transfer as Cosby et al. (1984) and HYPRES

                print *, "LnTo will be applied using the pedo-transfer function of Cosby et al. (1984) for non-organic soils."
                print *, " and HYPRES for organic"
                n_glob_pms = 7         !a_const, a_sand, a_clay,
                allocate(min_gp(7))
                allocate(max_gp(7))

                !cosby et al. limits
                min_gp(1) = -3.5 !-2.5!-3.5 increased to -2.5 following dotty plots
                min_gp(2)= 0.006!-0.07
                min_gp(3) = -0.02!-0.07

                max_gp(1) = 0.3!0.5
                max_gp(2)= 0.03!0.01
                max_gp(3) = -0.0032!0.01

                !HYPRES limits Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                !bounds changed following experiments in excel and matlab
                min_gp(4) = 1 !3
                max_gp(4) = 6!11

                min_gp(5) = -1.3!-0.5
                max_gp(5) = -0.5!-1.3

                min_gp(6) = -0.1
                max_gp(6) = -0.003

                min_gp(7) = -0.4!-0.8
                max_gp(7) = 0!-0.24

            case(3)     !Same as case 2 for areas covering low/moderately productive geology,
                        ! higher values for areas covering high productivity geology.

                print *, "LnTo will be applied using the pedo-transfer function of Cosby et al. (1984) for non-organic soils."
                print *, " and HYPRES for organic, with a multiplier for highly productive geology"
                n_glob_pms = 8         !a_const, a_sand, a_clay,
                allocate(min_gp(8))
                allocate(max_gp(8))

                !cosby et al. limits
                min_gp(1) = -3.5
                min_gp(2)= 0.006
                min_gp(3) = -0.02

                max_gp(1) = 0.3
                max_gp(2)= 0.03
                max_gp(3) = -0.0032

                !HYPRES limits Ks = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
                !bounds changed following experiments in excel and matlab
                min_gp(4) = 1 !3
                max_gp(4) = 6!11
                min_gp(5) = -1.3!-0.5
                max_gp(5) = -0.5!-1.3
                min_gp(6) = -0.1
                max_gp(6) = -0.003
                min_gp(7) = -0.4!-0.8
                max_gp(7) = 0!-0.24

                !Addition for productive geology limits
                min_gp(8) = 0 !i.e. no effect.
                max_gp(8) = 10 !Value chosen looking at optimal parameters MC vs MPR for southeast.

            case default
                print *, "ERROR: pedo_tf_LnTo incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(2) = n_glob_pms
        allocate(glob_pms_LnTo(n_glob_pms, n_pm_maps))
        allocate(rand_nums(sum(n_gp_all), n_pm_maps+start_seed-1))

        !3. Generate lists of global parameters between min and max ranges, of length n
        !need to have list such that starting seed = 50 would give same result as 50th number from seed=1
        Call RANDOM_SEED( GET = seed )
        seed(1:12) = (/ 3,3,3,3,3,3,3,3,3,3,3,3/)
        Call RANDOM_SEED( PUT = seed )
        call RANDOM_NUMBER ( rand_nums )

        !now normalise random numbers to within bounds
        DO i = 1, n_glob_pms
            DO j = start_seed, start_seed+n_pm_maps-1

                rn_i = (sum(n_gp_all) - n_glob_pms) + i                 !rn_i is index for the random number
                num = rand_nums(rn_i,j)                                 !the appropriate random number, in range 0-1
                num = (num * (max_gp(i)-min_gp(i)))+min_gp(i)           !normalise the random number to given min/max
                glob_pms_LnTo(i, j-start_seed+1) = num

            END DO
            IF (n_pm_maps >= 2) THEN
                print *, 'First 2 values for parameter ',i,' : ',glob_pms_LnTo(i,1:2)
            END IF
        END DO


        !print *,



    end subroutine init_gp_LnTo





    ! **********************************************************************************************
    ! subroutine: pedotf_LnTo: Pedo-transfer function for LnTo
    ! **********************************************************************************************
    ! This routine takes n_i, signalling that we are now calculating parameter map n_i out of n,
    ! It takes the list of global parameters and required basin predictors, and from that
    ! applies the pedo-transfer functions to produce a parameter map.

    subroutine pedotf_LnTo(n_i, pedo_tf, glob_pms, sand_0_10, clay_0_10, pm_map,bp_maps)

        implicit none

        !declare dummy variables
        integer, intent(in)                                          :: n_i         !param map number to do
        integer, intent(in)                                          :: pedo_tf     !equation selection
        double precision, dimension(:,:), allocatable, intent(in)    :: glob_pms    !list of global params
        real, dimension(:,:), allocatable, intent(in)    :: sand_0_10   !input basin predictor maps
        real, dimension(:,:), allocatable, intent(in)    :: clay_0_10
        double precision, dimension(:,:), allocatable, intent(out)   :: pm_map      !output parameter map
        real, allocatable, dimension(:,:,:), intent(in)              :: bp_maps     !basin predictor maps.(x,y,bp)

        !declare local variables
        integer             :: n_glob_pms       !number of global parameters
        integer             :: nrows_bp         !number of rows in basin predictor map
        integer             :: ncols_bp
        integer             :: n_pm_maps        !number of parameter maps - defined by ncols in glob_pms_LnTo
        integer             :: i
        integer             :: j

        !find number of rows and columns in basin predictor files - pm map should be same size
        nrows_bp = size(sand_0_10(:,1))
        ncols_bp = size(sand_0_10(1,:))
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


            case (2)
            ! if non-organic soils   !Use pedo-transfer as Cosby et al. (1984), LnTo = a1 + a2 * %sand + a3 * %clay
            ! if organic soils       !use HYPRES ptf excluding silt and clay, lnto = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)

             n_glob_pms = 7
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        !check if organic
                        !IF (bp_maps(i,j,16) == 1) THEN
                        IF (bp_maps(i,j,4)>=35) THEN  !now using organic threshold of more than 35%
                            !HYPRES function for organic soils
                            pm_map(i,j) = glob_pms(4,n_i)+(glob_pms(5,n_i)*bp_maps(i,j,15)*bp_maps(i,j,15))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(6,n_i)*(1/(bp_maps(i,j,4))))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(7,n_i)*bp_maps(i,j,15)*bp_maps(i,j,4))
                        ELSE
                            !cosby et al. function for mineral soils
                            pm_map(i,j) = glob_pms(1, n_i) + (glob_pms(2,n_i)*sand_0_10(i,j))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(3,n_i)*clay_0_10(i,j))
                        END IF

                    END DO
                END DO

            case (3)
            ! if non-organic soils   !Use pedo-transfer as Cosby et al. (1984), LnTo = a1 + a2 * %sand + a3 * %clay
            ! if organic soils       !use HYPRES ptf excluding silt and clay, lnto = a + b(BulkDensity)^2 + c(OrganicMatter)^-1 + d(BulkDensity * OrganicMatter)
            ! if highly productive hydrogeology !add an extra global parameter.
             n_glob_pms = 8
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        !check if organic
                        IF (bp_maps(i,j,4)>=35) THEN  !now using organic threshold of more than 35%
                            !HYPRES function for organic soils
                            pm_map(i,j) = glob_pms(4,n_i)+(glob_pms(5,n_i)*bp_maps(i,j,15)*bp_maps(i,j,15))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(6,n_i)*(1/(bp_maps(i,j,4))))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(7,n_i)*bp_maps(i,j,15)*bp_maps(i,j,4))
                        ELSE
                            !cosby et al. function for mineral soils
                            pm_map(i,j) = glob_pms(1, n_i) + (glob_pms(2,n_i)*sand_0_10(i,j))
                            pm_map(i,j) = pm_map(i,j) + (glob_pms(3,n_i)*clay_0_10(i,j))
                        END IF

                        !add extra global parameter if high productivity hydrogeology
                        IF (bp_maps(i,j,18) == 1) THEN
                            pm_map(i,j) = pm_map(i,j) + glob_pms(8,n_i)
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
                IF (pm_map(i,j)<-15) THEN
                    pm_map(i,j)=-15
                ELSEIF (pm_map(i,j)>10) THEN
                    pm_map(i,j)=10
                END IF
            END DO
        END DO

    end subroutine pedotf_LnTo




end module mpr_LnTo


