!! Module containing subroutines for MPR.f90
!! Module for GW model
!! This module contains all subroutines involving MPR for the B&Ks parameter:
!
!   get_BasPred_SRmax    Complete logical list of basin predictors required for LnTo
!   init_gp_SRmax        Initialise global parameters for LnTo
!   pedotf_SRmax         Pedo-transfer function for LnTo

!! Yanchen Zheng - 4th Oct 2023

module mpr_B_KsGW
contains


    ! **********************************************************************************************
    ! subroutine: get_BasPred_LnTo : Complete logical list of basin predictors required for LnTo
    ! **********************************************************************************************
    ! Given the user-specified pedo transfer function selected for LnTo, this subroutine adds to a
    ! logical list of required basin predictors, by specifying which predictors it requires.

    subroutine get_soilclass(sand, silt,clay,soil_class)

        implicit none
        real, allocatable, dimension(:,:),intent(in)    :: sand
        real, allocatable, dimension(:,:),intent(in)    :: silt
        real, allocatable, dimension(:,:),intent(in)    :: clay
        real, allocatable, dimension(:,:),intent(out)    :: soil_class
        real, allocatable, dimension(:,:)    :: sand_0_10
        real, allocatable, dimension(:,:)    :: silt_0_10
        real, allocatable, dimension(:,:)    :: clay_0_10
        integer i,j

        allocate(soil_class(size(sand,1),size(sand,2)))
        allocate(sand_0_10(size(sand,1),size(sand,2)))
        allocate(silt_0_10(size(sand,1),size(sand,2)))
        allocate(clay_0_10(size(sand,1),size(sand,2)))

        do i = 1, size(sand,1)
            do j = 1, size(sand,2)
                if (sand(i, j)<0) then
                sand_0_10(i, j)=-9999
                silt_0_10(i, j)=-9999
                clay_0_10(i, j)=-9999
                else
                
                    sand_0_10(i, j) = sand(i, j) / 100.0
                    silt_0_10(i, j) = silt(i, j) / 100.0
                    clay_0_10(i, j) = clay(i, j) / 100.0
                endif
            end do
        end do

       

        do i=1,size(sand_0_10,1)
        do j=1,size(sand_0_10,2)
            if (sand_0_10(i,j)<0) then
                soil_class(i,j)=sand_0_10(i,j) !-9999
            else if (silt_0_10(i,j)+1.5*clay_0_10(i,j)<0.15) then
                soil_class(i,j)=1 !'Sand'
            else if ((silt_0_10(i,j)+1.5*clay_0_10(i,j))>=0.15.AND.(silt_0_10(i,j)+2*clay_0_10(i,j))<0.3) then
                soil_class(i,j)=2 !'Loamy Sand'
            else if (clay_0_10(i,j)>=0.07 .AND. (clay_0_10(i,j)<=0.2) .AND. (sand_0_10(i,j)>0.52)&
             .AND. (silt_0_10(i,j)+2*clay_0_10(i,j))>=0.3) then
                soil_class(i,j)=3 !'Sandy Loam'
            else if (clay_0_10(i,j)<0.07 .AND. (silt_0_10(i,j) < 0.5) .AND. (silt_0_10(i,j)+2*clay_0_10(i,j))>=0.3) then
                soil_class(i,j)=3 !'Sandy Loam'
            else if (clay_0_10(i,j)>=0.07 .AND. (clay_0_10(i,j)<=0.27) .AND. (silt_0_10(i,j)>=0.28).AND. & 
                  (silt_0_10(i,j)<0.5) .AND. sand_0_10(i,j)<=0.52) then
                soil_class(i,j)=5 !'Loam'
            else if ((silt_0_10(i,j)>=0.5 .AND. clay_0_10(i,j)>=0.12 .AND. clay_0_10(i,j)<0.27) .or. &
                (silt_0_10(i,j)>=0.5 .AND. silt_0_10(i,j)<0.8 .AND. clay_0_10(i,j)<0.12)) then
                soil_class(i,j)=4 !'SILT LOAM'
            else if (silt_0_10(i,j)>=0.8 .AND. clay_0_10(i,j)<0.12) then
                !'SILT' !might be useless
                print*,'silt cells found! setting to silt loam'
                soil_class(i,j)=4 !'SILT' !might be useless
            else if (clay_0_10(i,j)>=0.2 .AND. clay_0_10(i,j)<0.35 .AND. silt_0_10(i,j)<0.28 .AND. sand_0_10(i,j)>0.45) then
                soil_class(i,j)=6 !'SANDY CLAY LOAM'
            else if (clay_0_10(i,j)>=0.27 .AND. clay_0_10(i,j) <0.4 .AND.sand_0_10(i,j)>0.2 .AND. sand_0_10(i,j)<=0.45) then
                soil_class(i,j)=8 !'CLAY LOAM'
            else if (clay_0_10(i,j)>=0.27 .AND. clay_0_10(i,j)<0.4 .AND. sand_0_10(i,j)<=0.2) then
                soil_class(i,j)=7 !'SILTY CLAY LOAM'
            else if (clay_0_10(i,j)>=0.35 .AND. sand_0_10(i,j)>=0.45) then
                soil_class(i,j)=9 !'SANDY CLAY'
            else if (clay_0_10(i,j)>=0.4 .AND. silt_0_10(i,j)>=0.4) then
                soil_class(i,j)=10 !'Silty CLAY'
            else if (clay_0_10(i,j)>= 0.4 .AND. sand_0_10(i,j)<=0.45.AND. silt_0_10(i,j)<0.4) then
                soil_class(i,j)=11 !'CLAY'
            end if
          end do    
         end do

    end subroutine get_soilclass

     ! **********************************************************************************************
    ! subroutine: init_gp_SRmax : Initialise global parameters for SRmax
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

    subroutine init_gp_B(pedo_tf_B,glob_pms_B, n_pm_maps, start_seed, n_gp_all,pm_range,&
        B_Ks_table)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_B
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_B      !global param list
        integer, intent(in)                                         :: n_pm_maps
        integer, intent(in)                                         :: start_seed   !to start rand num generator
        integer, dimension(11), intent(inout)                        :: n_gp_all     !vector of number of global params for all params
        double precision, dimension(3), intent(in)                  :: pm_range    !min, fixed, max values for this parameter
        double precision, allocatable,dimension(:,:),intent(in)     :: B_Ks_table

        !declare local variables
        integer                                        :: n_glob_pms   !number of global parameters
        double precision, allocatable, dimension(:)    :: min_gp    !min values for global parameters
        double precision, allocatable, dimension(:)    :: max_gp
        real, dimension(:,:),allocatable               :: rand_nums
        double precision                               :: num
        integer, dimension(12)                         :: seed
        integer                                        :: i,j, rn_i


        !1. Define how many global parameters should exist and min/max ranges, based on selected transfer function
        select case(pedo_tf_B)
            
            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "B is a fixed parameter, of value ",pm_range(2)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(2)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(2)

            !2024 March sensitivity analysis
             case (1)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "B is a global parameter, in the range ",pm_range(1)," ",pm_range(3)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(1)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(3)    

            !B lookup table
            case (2)    !Pedotransfer function!!!

                print *, "B will use a pedotransfer function based on soil texture table"
                print *, "B is a if soil class=1, b if soil class=2, *c "

                !allocate number of global parameters

                n_glob_pms = size(B_Ks_table(:,1))
                allocate(min_gp(n_glob_pms))
                allocate(max_gp(n_glob_pms))

                !set global parameter ranges - from soil rooting zone table
                DO i = 1,n_glob_pms
                    min_gp(i) = B_Ks_table(i,2)
                    max_gp(i) = B_Ks_table(i,3)
                 END DO

            case default
                print *, "ERROR: pedo_tf_B incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(8) = n_glob_pms

        allocate(glob_pms_B(n_glob_pms, n_pm_maps))
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

                rn_i = (sum(n_gp_all) - n_glob_pms )+ i                 !rn_i is index for the random number
                num = rand_nums(rn_i,j)                                 !the appropriate random number, in range 0-1
                num = (num * (max_gp(i)-min_gp(i)))+min_gp(i)           !normalise the random number to given min/max
                glob_pms_B(i, j-start_seed+1) = num

            END DO
            IF (n_pm_maps >= 2) THEN
                print *, 'First 2 values for parameter ',i,' : ',glob_pms_B(i,1:2)
            END IF
        END DO


    end subroutine init_gp_B


    subroutine init_gp_Ks(pedo_tf_Ks,glob_pms_Ks, n_pm_maps, start_seed, n_gp_all,pm_range,&
        B_Ks_table)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_Ks
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_Ks      !global param list
        integer, intent(in)                                         :: n_pm_maps
        integer, intent(in)                                         :: start_seed   !to start rand num generator
        integer, dimension(11), intent(inout)                        :: n_gp_all     !vector of number of global params for all params
        double precision, dimension(3), intent(in)                  :: pm_range    !min, fixed, max values for this parameter
        double precision, allocatable,dimension(:,:),intent(in)     :: B_Ks_table

        !declare local variables
        integer                                        :: n_glob_pms   !number of global parameters
        double precision, allocatable, dimension(:)    :: min_gp    !min values for global parameters
        double precision, allocatable, dimension(:)    :: max_gp
        real, dimension(:,:),allocatable               :: rand_nums
        double precision                               :: num
        integer, dimension(12)                         :: seed
        integer                                        :: i,j, rn_i


        !1. Define how many global parameters should exist and min/max ranges, based on selected transfer function
        select case(pedo_tf_Ks)
            
            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "Ks is a fixed parameter, of value ",pm_range(2)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(2)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(2)

            !2024 March sensitivity analysis
             case (1)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "Ks is a global parameter, in the range ",pm_range(1)," ",pm_range(3)
                n_glob_pms = 1
                allocate(min_gp(1))         !min and max values have length 1, as we have 1 global parameter
                allocate(max_gp(1))
                min_gp(1) = pm_range(1)    !min and max values set to same fixed parameter value.
                max_gp(1) = pm_range(3)     

            !Ks lookup table
            case (2)    !Pedotransfer function!!!

                print *, "Ks will use a pedotransfer function based on soil texture table"
                print *, "Ks is a if soil class=1, b if soil class=2, *c "

                !allocate number of global parameters

                n_glob_pms = size(B_Ks_table(:,1))
                allocate(min_gp(n_glob_pms))
                allocate(max_gp(n_glob_pms))

                !set global parameter ranges - from soil rooting zone table
                DO i = 1,n_glob_pms
                    min_gp(i) = B_Ks_table(i,4)
                    max_gp(i) = B_Ks_table(i,5)
                 END DO

            case default
                print *, "ERROR: pedo_tf_Ks incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(9) = n_glob_pms

        allocate(glob_pms_Ks(n_glob_pms, n_pm_maps))
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

                rn_i = (sum(n_gp_all) - n_glob_pms )+ i                 !rn_i is index for the random number
                num = rand_nums(rn_i,j)                                 !the appropriate random number, in range 0-1
                num = (num * (max_gp(i)-min_gp(i)))+min_gp(i)           !normalise the random number to given min/max
                glob_pms_Ks(i, j-start_seed+1) = num

            END DO
            IF (n_pm_maps >= 2) THEN
                print *, 'First 2 values for parameter ',i,' : ',glob_pms_Ks(i,1:2)
            END IF
        END DO

    end subroutine init_gp_Ks


    ! **********************************************************************************************
    ! subroutine: pedotf_B: Pedo-transfer function for SRmax
    ! **********************************************************************************************
    ! This routine takes n_i, signalling that we are now calculating parameter map n_i out of n,
    ! It takes the list of global parameters and required basin predictors, and from that
    ! applies the pedo-transfer functions to produce a parameter map.

    subroutine pedotf_B(n_i, pedo_tf, glob_pms,bp, pm_map)

        implicit none

        !declare dummy variables
        integer, intent(in)                                          :: n_i         !param map number to do
        integer, intent(in)                                          :: pedo_tf     !equation selection
        double precision, dimension(:,:), allocatable, intent(in)    :: glob_pms    !list of global params
        real, dimension(:,:), allocatable, intent(in)                :: bp          !basin predictor - soil class map
        double precision, dimension(:,:), allocatable, intent(out)   :: pm_map      !output parameter map
        

        !declare local variables
        integer             :: n_glob_pms       !number of global parameters
        integer             :: n_pm_maps        !number of parameter maps - defined by ncols in glob_pms_LnTo
        integer             :: i
        integer             :: j
        integer             :: nrows_bp         !number of rows in basin predictor map
        integer             :: ncols_bp
        integer             :: lu_class
        logical             :: warning_sent


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

             case (1)    !global parameter

                n_glob_pms = 1
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        pm_map(i,j) = glob_pms(1, n_i)

                    END DO
                END DO   


            case (2)   !pedotransfer! soil class table

                n_glob_pms = 11
                warning_sent = .false.
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        lu_class = int(bp(i,j))

                        !deal with any -9999 values without breaking code
                        IF (lu_class == -9999) THEN
                            IF (warning_sent .eqv. .false.) THEN
                                print *,'WARNING: some soil classes are -9999, setting B to 0.01 here'
                                warning_sent = .true.
                            END IF
                            pm_map(i,j) = 0.01

                         ELSEIF (lu_class == 0) THEN
                            IF (warning_sent .eqv. .false.) THEN
                                print *,'WARNING: some soil classes are 0, setting B to 0.01 here'
                                warning_sent = .true.
                            END IF
                            pm_map(i,j) = 0.01
        
                        ELSE

                       !if no problems with missing data then apply pedotransfer function!
                            pm_map(i,j) = glob_pms(lu_class,n_i) 
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
                IF (pm_map(i,j)<0.01) THEN
                    pm_map(i,j)=0.01
                ELSEIF (pm_map(i,j)>23.75) THEN
                    pm_map(i,j)=23.75
                END IF

                
            END DO
        END DO

    end subroutine pedotf_B


     subroutine pedotf_Ks(n_i, pedo_tf, glob_pms,bp, pm_map,pm_map_Drz,pm_map_szm)

        implicit none

        !declare dummy variables
        integer, intent(in)                                          :: n_i         !param map number to do
        integer, intent(in)                                          :: pedo_tf     !equation selection
        double precision, dimension(:,:), allocatable, intent(in)    :: glob_pms    !list of global params
        real, dimension(:,:), allocatable, intent(in)                :: bp          !basin predictor - soil class map
        double precision, dimension(:,:), allocatable, intent(in)                :: pm_map_Drz
        double precision, dimension(:,:), allocatable, intent(in)                :: pm_map_szm
        double precision, dimension(:,:), allocatable, intent(out)   :: pm_map      !output parameter map
        

        !declare local variables
        integer             :: n_glob_pms       !number of global parameters
        integer             :: n_pm_maps        !number of parameter maps - defined by ncols in glob_pms_LnTo
        integer             :: i
        integer             :: j
        integer             :: nrows_bp         !number of rows in basin predictor map
        integer             :: ncols_bp
        integer             :: lu_class
        logical             :: warning_sent


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

            case (1)    !global parameter

                n_glob_pms = 1
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        pm_map(i,j) = glob_pms(1, n_i)

                    END DO
                END DO

            case (2)   !pedotransfer! soil class table

                n_glob_pms = 11
                warning_sent = .false.
                DO i = 1, nrows_bp
                    DO j = 1, ncols_bp

                        lu_class = int(bp(i,j))

                        !deal with any -9999 values without breaking code
                        IF (lu_class == -9999) THEN
                            IF (warning_sent .eqv. .false.) THEN
                                print *,'WARNING: some soil classes are -9999, setting Ks to 0.001 here'
                                warning_sent = .true.
                            END IF
                            pm_map(i,j) = 0.001

                         ELSEIF (lu_class == 0) THEN
                            IF (warning_sent .eqv. .false.) THEN
                                print *,'WARNING: some soil classes are 0, setting Ks to 0.001 here'
                                warning_sent = .true.
                            END IF
                            pm_map(i,j) = 0.001
        
                        ELSE

                       !if no problems with missing data then apply pedotransfer function!
                            !Ks look up table
                            pm_map(i,j) = glob_pms(lu_class,n_i) 
                            !Ks root zone bottom
                            !pm_map(i,j) = glob_pms(lu_class,n_i)*exp(-pm_map_Drz(i,j)/pm_map_szm(i,j))
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
                IF (pm_map(i,j)<0.00001) THEN
                    pm_map(i,j)=0.00001
                ELSEIF (pm_map(i,j)>15.21) THEN
                    pm_map(i,j)=15.21
                END IF

                
            END DO
        END DO

    end subroutine pedotf_Ks





end module mpr_B_KsGW