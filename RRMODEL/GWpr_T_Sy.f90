!! Module containing subroutines for MPR.f90
!! Module for GW model
!! This module contains all subroutines involving MPR for the T&Sy parameter:
!
!   get_BasPred_SRmax    Complete logical list of basin predictors required for LnTo
!   init_gp_SRmax        Initialise global parameters for LnTo
!   pedotf_SRmax         Pedo-transfer function for LnTo

!! Yanchen Zheng - 21th Feb 2024

module GWpr_T_Sy
contains

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

    subroutine init_gp_T(pedo_tf_T,glob_pms_T, n_pm_maps, start_seed, n_gp_all,&
        T_Sy_table,GW_T_clip,GW_T)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_T
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_T      !global param list
        integer, intent(in)                                         :: n_pm_maps
        integer, intent(in)                                         :: start_seed   !to start rand num generator
        integer, dimension(11), intent(inout)                        :: n_gp_all     !vector of number of global params for all params
        !double precision, dimension(3), intent(in)                  :: pm_range    !min, fixed, max values for this parameter
        double precision, allocatable,dimension(:,:),intent(in)     :: T_Sy_table
        double precision, allocatable, dimension(:,:) :: GW_T_clip, GW_T
        !declare local variables
        integer                                        :: n_glob_pms   !number of global parameters
        double precision, allocatable, dimension(:)    :: min_gp    !min values for global parameters
        double precision, allocatable, dimension(:)    :: max_gp
        real, dimension(:,:),allocatable               :: rand_nums
        double precision                               :: num
        integer, dimension(12)                         :: seed
        integer                                        :: i,j, rn_i


        !1. Define how many global parameters should exist and min/max ranges, based on selected transfer function
        select case(pedo_tf_T)
            
            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "T is a fixed parameter (estimated clip map) "
                GW_T=GW_T_clip
                n_glob_pms = 1

                return

            !B lookup table
            case (2)    !Pedotransfer function!!!

                print *, "T will use a pedotransfer function based on geology table"
                print *, "T is a if geology class=1, b if geology class=2, *c "

                !allocate number of global parameters

                n_glob_pms = size(T_Sy_table(:,1))
                allocate(min_gp(n_glob_pms))
                allocate(max_gp(n_glob_pms))

                !set global parameter ranges - from soil rooting zone table
                DO i = 1,n_glob_pms
                    min_gp(i) = T_Sy_table(i,2)
                    max_gp(i) = T_Sy_table(i,3)
                 END DO

            case default
                print *, "ERROR: pedo_tf_T incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(10) = n_glob_pms

        allocate(glob_pms_T(n_glob_pms, n_pm_maps))
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
                glob_pms_T(i, j-start_seed+1) = num

            END DO
            ! IF (n_pm_maps >= 2) THEN
            !     print *, 'First 2 values for parameter ',i,' : ',glob_pms_T(i,1:2)
            ! END IF
        END DO



    end subroutine init_gp_T




subroutine init_gp_Sy(pedo_tf_Sy,glob_pms_Sy, n_pm_maps, start_seed, n_gp_all,&
        T_Sy_table,GW_Sy_clip,GW_Sy)

        use dta_utility

        implicit none

        !declare dummy variables
        integer, intent(in)                                         :: pedo_tf_Sy
        double precision, allocatable, dimension(:,:),intent(inout) :: glob_pms_Sy      !global param list
        integer, intent(in)                                         :: n_pm_maps
        integer, intent(in)                                         :: start_seed   !to start rand num generator
        integer, dimension(11), intent(inout)                        :: n_gp_all     !vector of number of global params for all params
        !double precision, dimension(3), intent(in)                  :: pm_range    !min, fixed, max values for this parameter
        double precision, allocatable,dimension(:,:),intent(in)     :: T_Sy_table
        double precision, allocatable, dimension(:,:) :: GW_Sy_clip, GW_Sy
        !declare local variables
        integer                                        :: n_glob_pms   !number of global parameters
        double precision, allocatable, dimension(:)    :: min_gp    !min values for global parameters
        double precision, allocatable, dimension(:)    :: max_gp
        real, dimension(:,:),allocatable               :: rand_nums
        double precision                               :: num
        integer, dimension(12)                         :: seed
        integer                                        :: i,j, rn_i


        !1. Define how many global parameters should exist and min/max ranges, based on selected transfer function
        select case(pedo_tf_Sy)
            
            case (0)    !Treat as fixed parameter (for ease this is a global parameter with a 0 range)

                print *, "Sy is a fixed parameter (estimated clip map) "
                GW_Sy=GW_Sy_clip
                n_glob_pms = 1

                return

            !Sy lookup table
            case (2)    !Pedotransfer function!!!

                print *, "Sy will use a pedotransfer function based on geology table"
                print *, "Sy is a if geology class=1, b if geology class=2, *c "

                !allocate number of global parameters

                n_glob_pms = size(T_Sy_table(:,1))
                allocate(min_gp(n_glob_pms))
                allocate(max_gp(n_glob_pms))

                !set global parameter ranges - from soil rooting zone table
                DO i = 1,n_glob_pms
                    min_gp(i) = T_Sy_table(i,4)
                    max_gp(i) = T_Sy_table(i,5)
                 END DO

            case default
                print *, "ERROR: pedo_tf_Sy incorrectly defined. "

        end select

        !2. From above, set sizes of the global parameter array and number of random numbers to produce.
        n_gp_all(11) = n_glob_pms

        allocate(glob_pms_Sy(n_glob_pms, n_pm_maps))
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
                glob_pms_Sy(i, j-start_seed+1) = num

            END DO
            ! IF (n_pm_maps >= 2) THEN
            !     print *, 'First 2 values for parameter ',i,' : ',glob_pms_Sy(i,1:2)
            ! END IF
        END DO



    end subroutine init_gp_Sy




    subroutine Tglob_pmsTomap(glob_pms_T,i_mc,Geo_index,GW_T,T_Sy_table)
        !YZ 2024
        !Generate T/Sy geo spatial map GW_T/Sy
        !according to glob_pms_t/Sy
        double precision, allocatable, dimension(:,:),intent(in) :: glob_pms_T
        integer :: i_mc,i,j,ind,t
        double precision, allocatable, dimension(:,:) :: GEO_index
        double precision, allocatable, dimension(:,:)   :: T_Sy_table !YZ 2024
        double precision, allocatable, dimension(:,:):: GW_T


        !loop through the GEO_index maps and gives correponding par values
        do i=1,size(GEO_index,1)
            do j=1,size(GEO_index,2)

                if (GEO_index(i,j)>=0) then
                   
                   !find the rows of the index
                   do t=1,size(glob_pms_T,1)
                      if (T_Sy_table(t,1)==GEO_index(i,j)) then
                         ind = t 
                         !return
                      end if
                   end do
                   
                   GW_T(i,j)=glob_pms_T(ind,i_mc)
                else
                   !if Geo_index is -9999 missing value
                   GW_T(i,j)=1e-06
                end if
            end do 
        end do 

    end subroutine Tglob_pmsTomap


subroutine Syglob_pmsTomap(glob_pms_Sy,i_mc,Geo_index,GW_Sy,T_Sy_table)
        !YZ 2024
        !Generate T/Sy geo spatial map GW_T/Sy
        !according to glob_pms_t/Sy
        double precision, allocatable, dimension(:,:),intent(in) :: glob_pms_Sy
        integer :: i_mc,i,j,ind,t
        double precision, allocatable, dimension(:,:) :: GEO_index
        double precision, allocatable, dimension(:,:)   :: T_Sy_table !YZ 2024
        double precision, allocatable, dimension(:,:):: GW_Sy


        !loop through the GEO_index maps and gives correponding par values
        do i=1,size(GEO_index,1)
            do j=1,size(GEO_index,2)

                if (GEO_index(i,j)>=0) then
                   
                   !find the rows of the index
                   do t=1,size(glob_pms_Sy,1)
                      if (T_Sy_table(t,1)==GEO_index(i,j)) then
                         ind = t 
                         !return
                      end if
                   end do
                   
                   GW_Sy(i,j)=glob_pms_Sy(ind,i_mc)
                else
                   !if Geo_index is -9999 missing value
                   GW_Sy(i,j)=1e-08
                end if
            end do 
        end do 

    end subroutine Syglob_pmsTomap






end module GWpr_T_Sy