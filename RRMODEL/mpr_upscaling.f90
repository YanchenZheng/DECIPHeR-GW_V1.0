!! Module containing subroutines for MPR.f90
!! This module contains all subroutines involving in parameter upscaling
!
!   extract_pm_HRU :    extracts list of params from all cells of selected HRU
!   upscale_harmonic :  upscales a list of parameters using the harmonic mean
!   upscale_arithmetic: upscales a list of parameters using the arithmetic mean

!! Rosie Lane - 31st August 2017

module mpr_upscaling
contains



    ! **********************************************************************************************
    ! subroutine: extract_pm_HRU : extracts list of params from all cells of selected HRU
    ! **********************************************************************************************
    !  extract_pm_HRU aims to:
    !  1. Loop through all values in a HRU and parameter value file
    !  2. For a specified HRU number, extract a vector of all parameter values within that HRU
    !  3. Output vector of parameter values

    subroutine extract_pm_HRU(HRU, HRU_map, pm_map, pms_HRU)

        implicit none

        !declare dummy variables
        integer, intent(in)                          :: HRU                 !id of current HRU
        double precision, dimension(:,:), intent(in) :: HRU_map             !map of HRUs
        double precision, dimension(:,:), intent(in) :: pm_map              !map of MPR parameters
        double precision, allocatable,dimension(:), intent(out) :: pms_HRU  !output list of params in HRU

        !declare local variables
        integer                 :: nrows_HRU    !dimensions of HRU & ASSUMED dimensions of pm maps
        integer                 :: ncols_HRU
        integer                 :: i
        integer                 :: j
        integer                 :: n_pms_HRU     !number of cells in current HRU


        !define local variables from HRU map
        nrows_HRU = size(HRU_map(:,1))
        ncols_HRU = size(HRU_map(1,:))

        !initialise
        n_pms_HRU = 0

        ! Gemma - n_pms_HRU2 = sum(HRU_map[, HRU_map.eq.HRU])
        ! Find the number of HRU squares matching the given HRU value.
!        DO i = 1, nrows_HRU
!            DO j = 1, ncols_HRU
!                IF (HRU_map(i,j)==HRU) THEN     !When you index a variable of the correct HRU...
!                    n_pms_HRU = n_pms_HRU + 1
!                END IF
!            END DO
!        END DO

        ! Find the number of HRU squares matching the given HRU value
        n_pms_HRU = int(SUM(HRU_map,MASK=HRU_map==HRU))/HRU

        !Allocate size of pms_HRU to prevent segmentation fault
        allocate(pms_HRU(n_pms_HRU))

!        !Loop through HRU_map again, extracting all pm_map values in the correct HRU
!        DO i = nrows_HRU,1,-1
!            DO j = ncols_HRU,1,-1
!                IF (HRU_map(i,j)==HRU) THEN     !When you index a variable of the correct HRU...
!                    pms_HRU(n_pms_HRU) = pm_map(i,j)  !...add the parameter value into pms_HRU.
!                    n_pms_HRU = n_pms_HRU -1
!                END IF
!            END DO
!        END DO

        !Extract all pm values within the correct HRU
        pms_HRU = PACK(pm_map, HRU_map==HRU)


    end subroutine extract_pm_HRU





    ! **********************************************************************************************
    ! subroutine: upscale_harmonic : upscales a list of parameters using the harmonic mean
    ! **********************************************************************************************
    !   upscale_harmonic aims to:
    !   take input of a vector of parameter values
    !   output the harmonic mean of those values

    subroutine upscale_harmonic(pms_HRU, upscaled_pm)

        implicit none

        double precision, dimension(:), intent(in)  :: pms_HRU      !List of params from extract_pm_HRU
        double precision, intent(out)               :: upscaled_pm  !output upscaled param
        double precision, dimension(:), allocatable :: recip_pms
        double precision                            :: sum_pms

        recip_pms = 1 / pms_HRU     !take reciprocal of all params
        sum_pms = SUM(recip_pms)    !take sum of recipocal params
        upscaled_pm = SIZE(pms_HRU) / sum_pms  !H = n params / sum recip params

    end subroutine upscale_harmonic


    ! **********************************************************************************************
    ! subroutine: upscale_arithmetic: upscales a list of parameters using the arithmetic mean
    ! **********************************************************************************************

    subroutine upscale_arithmetic(pms_HRU, upscaled_pm)

        implicit none

        double precision, dimension(:), intent(in) :: pms_HRU       !List of params from extract_pm_HRU
        double precision, intent(out)              :: upscaled_pm   !output upscaled param
        double precision                           :: sum_pms
        integer :: npoints,i

        npoints = SIZE(pms_HRU)
        !do i = 1,SIZE(pms_HRU)
        !    if(pms_HRU(i)==-9999)then
        !        pms_HRU(i)=0
        !        npoints = npoints-1
        !    end if
        !end

        sum_pms = SUM(pms_HRU)                  !take sum of params

        !if (npoints>0) then
            upscaled_pm = sum_pms / npoints  !A = sum params / n params
        !else
        !    print *, 'Error in upscale arithmetic - all points were NaN!!!!!!!    pm set to 1'
        !    upscaled_pm = 1
        !end

    end subroutine upscale_arithmetic

    ! **********************************************************************************************
    ! subroutine: upscale_geometric: upscales a list of parameters using the geometric mean
    ! **********************************************************************************************

    subroutine upscale_geometric(pms_HRU, upscaled_pm)

        implicit none

        double precision, dimension(:), intent(in) :: pms_HRU       !List of params from extract_pm_HRU
        double precision, intent(out)              :: upscaled_pm   !output upscaled param
        double precision                           :: prod_pms

        prod_pms = PRODUCT(pms_HRU)                     !take product of params
        upscaled_pm = prod_pms**(1/SIZE(pms_HRU))        !G = (product(params))^(1/n)

    end subroutine upscale_geometric

    ! **********************************************************************************************
    ! subroutine: upscale_majority: upscales selecting the most frequently occuring value
    ! **********************************************************************************************

    subroutine upscale_majority(pms_HRU, upscaled_pm)

        implicit none

        double precision, dimension(:), intent(in) :: pms_HRU       !List of params from extract_pm_HRU
        double precision, intent(out)              :: upscaled_pm   !output upscaled param

        !local variables
        double precision                           :: prod_pms
        integer,dimension(:),allocatable           :: freq_occ  !frequency of occurance of each parameter value
        integer                                    :: i,j

        allocate(freq_occ(size(pms_HRU)))

        !find the frequency of occurance of all parameter values
        DO i = 1,size(pms_HRU)
            freq_occ(i) = 0
            DO j = 1,size(pms_HRU)
                IF (pms_HRU(i)==pms_HRU(j)) THEN
                    freq_occ(i) = freq_occ(i)+1
                END IF
            END DO
        END DO

        !return parameter value that occurs most frequently
        DO i = 1,size(pms_HRU)
            IF (freq_occ(i) == MAXVAL(freq_occ)) THEN
                upscaled_pm = pms_HRU(i)
            END IF
        END DO

    end subroutine upscale_majority






end module mpr_upscaling


