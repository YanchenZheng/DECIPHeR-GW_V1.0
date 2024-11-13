!! Module containing subroutines for MPR.f90
!! This module contains all subroutines for reading, clipping and writing Basin predictors

!! Rosie Lane - 31st August 2017
!
!   Subroutines:
!   Read_clip_bp        Read basin predictor file, and clip to the same area as HRU file
!   Read_req_soils       Reads in all required soils, writes clipped grid if required


module mpr_extract_BasPred
    contains


    ! **********************************************************************************************
    ! subroutine: Read_clip_bp : Read basin predictor file, and clip to the same area as HRU file
    ! **********************************************************************************************
    ! This subroutine reads in a basin predictor grid, specified by filename, and clips it to the
    ! required area covered by the HRU file given in stats_HRU, returning grid_clipped.

    subroutine Read_clip_bp(filename,grid_clipped, stats_HRU, HRU_map)

    use dta_utility

    implicit none

        !declare dummy variables
        character(len=*), intent(in)         :: filename     !filepath to stored .asc file of basin predictor
        double precision, dimension(:,:), allocatable, intent(out)  :: grid_clipped !output ascii grid
        double precision, dimension(4), intent(in)   :: stats_HRU    !HRU xll, yll, cellsize and nodata value
        double precision, dimension(:,:), intent(in)    :: HRU_map


        !declare local variables

        !for reading an ascii grid
        double precision, dimension(:,:), allocatable   :: grid
        double precision                                :: xll
        double precision                                :: yll
        double precision                                :: cellsize
        double precision                                :: nodata
        integer                                         :: ncol
        integer                                         :: nrow

        !HRU map stats
        integer                                         :: nrows_HRU
        integer                                         :: ncols_HRU

        !calculation of row/col equivalents for national grid locations
        integer     :: i_xll
        integer     :: i_yll
        integer     :: i_yul
        integer     :: i_xlr
        integer     :: x, y

        ! get HRU map stats
        nrows_HRU = size(HRU_map(:,1))
        ncols_HRU = size(HRU_map(1,:))

        ! read basin predictor file
        call read_ascii_grid(trim(filename),grid,ncol,nrow,xll,yll,cellsize,nodata)

        ! make a check to see if cellsizes differ - problems if so!!
        IF (cellsize /= stats_HRU(3)) THEN
            print *, "ERROR: Cellsize differs between HRU and basin predictor maps"
            stop
        END IF

        ! get row/col numbers of corners where HRU map and basin predictor map overlap
        i_xll = ((stats_HRU(1) - xll) / cellsize) + 1
        i_yll = nrow - ((stats_HRU(2) - yll) / cellsize)
        i_yul = i_yll - (nrows_HRU - 1)
        i_xlr = i_xll + ncols_HRU - 1

        ! allocate output ascii grid to be same size as HRU grid
        allocate(grid_clipped(nrows_HRU, ncols_HRU))

        ! loop through area where HRU map overlays basin predictor map, and extract basin pred map values
        DO x = 1, ncols_HRU-1
            DO y = 1, nrows_HRU-1
                grid_clipped(y,x) = grid(i_yul+y-1,i_xll+x-1)
            END DO
        END DO


        ! make sure nodata value is set to the same as HRU nodata value
        IF (nodata /= stats_HRU(4)) THEN
            print *, "Converting basin predictor nodata value to that of HRU map: ", int(stats_HRU(4))
            DO x = 1, ncols_HRU-1
                DO y = 1, nrows_HRU-1
                    IF (grid_clipped(y,x) == nodata) THEN
                        grid_clipped(y,x) = stats_HRU(4)
                    END IF
                END DO
            END DO
        END IF

        ! set all areas which are nodata in HRU file to nodata in basin predictor file
        ! - this can crop the map to the correct extent if later importing into arcmap etc.
        DO x = 1,ncols_HRU
            DO y = 1,nrows_HRU
                IF (HRU_map(y,x) == stats_HRU(4)) THEN !if nodata in HRU map
                    grid_clipped(y,x) = stats_HRU(4)
                END IF
            END DO
        END DO

        !print *, "Successfully read basin predictor: ",filename


    end subroutine Read_clip_bp



    ! **********************************************************************************************
    ! subroutine: Read_req_soils
    ! **********************************************************************************************
    ! This routine reads in all the required basin predictors (as specified in req_bp)
    ! All basin predictors are returned in the variable bp_maps
    ! If the user has specified writing of basin predictors, then all basin predictors required will
    !   also be written as an ascii file in the OUTPUT folder.
    subroutine read_req_soils(input_folder, &
        out_folder, &
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

        !use dta_MPR
        use dta_utility

        implicit none

        !declare dummy variables
        character(len=1024), intent(in)                 :: input_folder    !MPR input folder
        character(len=1024), intent(in)                 :: out_folder      !Folder to store output - currently MPRfolder
        character(len=1024), dimension(27),intent(in)   :: fnames_baspred
        logical, dimension(*), intent(in)               :: req_bp    !logical declaring whether each basin predictor is required
        double precision, dimension(4), intent(in)      :: stats_HRU    !HRU map xll, yll, cellsize and nodata value
        double precision, dimension(:,:), intent(in)    :: HRU_map
        integer, intent(in)                             :: save_bp_maps !value of 1 if clipped basin predictor maps should be saved
        real, allocatable, dimension(:,:,:), intent(out)    :: bp_maps  !array of all basin predictor ascii maps
        double precision, allocatable, dimension(:,:,:), intent(out)      :: soilmusiddata !array of all soils depth data tables.
        real, allocatable, dimension(:,:), intent(out) :: bp_root_depths    !landcover class rooting depth table
        double precision, allocatable, dimension(:,:), intent(out) :: B_Ks_table    !landcover class rooting depth table
        double precision, allocatable, dimension(:,:), intent(out) :: T_Sy_table    !Geology class

        !declare local variables
        double precision, allocatable, dimension(:,:)   :: sand_0_10
        double precision, allocatable, dimension(:,:)   :: silt_0_10
        double precision, allocatable, dimension(:,:)   :: clay_0_10
        double precision, allocatable, dimension(:,:)   :: orgm_0_10
        real, allocatable, dimension(:,:)   :: ksat_decline
        double precision, allocatable, dimension(:,:)   :: MUSIDs
        double precision, allocatable, dimension(:,:)   :: n_entries
        double precision, allocatable, dimension(:,:)   :: tempdata
        double precision, allocatable, dimension(:,:)   :: map
        character(len=1024)                             :: temp_fn
        double precision                                :: xllcorner_HRU
        double precision                                :: yllcorner_HRU
        double precision                                :: cellsize_HRU
        double precision                                :: nodata_HRU
        integer                                         :: ncols_HRU
        integer                                         :: nrows_HRU
        integer                                         :: i
        integer                                         :: nrows
        integer                                         :: ncols
        integer                                         :: nmusids

        !
        xllcorner_HRU = stats_HRU(1)
        yllcorner_HRU = stats_HRU(2)
        cellsize_HRU = stats_HRU(3)
        nodata_HRU = stats_HRU(4)
        ncols_HRU = size(HRU_map(1,:))
        nrows_HRU = size(HRU_map(:,1))

        !bp maps dimensions refer to (x, y, req_bp)
        allocate(bp_maps(nrows_HRU, ncols_HRU,18))

        ! Read in required basin predictor maps - and clip to size of HRU map
        !   bp 1    surface % sand
        !   bp 2    surface % silt
        !   bp 3    surface % clay
        !   bp 4    surface % organic
        !   bp 5    map of soil depth profiles
        !   bp 7-11 tables of soil depth profiles
        !   bp 12   Land cover map
        !   bp 13   Root depth table
        !   bp 14   Depth to bedrock map
        !   bp 15   surface bulk density
        !   bp 16   locations of organic (binary with 1=organic soil) (NOT USED - now threshold of 35% for organic soils)
        !   bp 17   porosity map
        !   bp 18   locations of high productive hydrogeology (binary, 1=highly productive bedrock)

        !Read in surface percentage sand - basin predictor 1
        IF (req_bp(1)) THEN
            print *, "Reading file ",trim(fnames_baspred(1))//".asc"
            call Read_clip_bp(trim(input_folder)//trim(fnames_baspred(1))//".asc",sand_0_10, stats_HRU, HRU_map)
            IF (save_bp_maps == 1) THEN
                temp_fn = trim(input_folder) // trim(fnames_baspred(1))//"_clipped.asc"
                print *, 'Writing clipped sand basin predictor file in input folder'
                !print *,
                call write_ascii_grid(temp_fn,sand_0_10,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
            bp_maps(:,:,1) = sand_0_10 !store all soils data in bp_maps for easy passing between functions
        END IF

        !Read in surface percentage silt - basin predictor 2
        IF (req_bp(2)) THEN
            print *, "Reading file ",trim(fnames_baspred(2))//".asc"
            call Read_clip_bp(trim(input_folder)//trim(fnames_baspred(2))//".asc",silt_0_10, stats_HRU, HRU_map)
            IF (save_bp_maps == 1) THEN
                temp_fn = trim(input_folder) // trim(fnames_baspred(2))//"_clipped.asc"
                print *, 'Writing clipped silt basin predictor file in input folder'
                !print *,
                call write_ascii_grid(temp_fn,silt_0_10,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
            bp_maps(:,:,2) = silt_0_10 !store all soils data in bp_maps for easy passing between functions
        END IF

        !Read in surface percentage clay - basin predictor 3
        IF (req_bp(3)) THEN
            print *, "Reading file ",trim(fnames_baspred(3))//".asc"
            call Read_clip_bp(trim(input_folder)//trim(fnames_baspred(3))//".asc",clay_0_10, stats_HRU, HRU_map)
            IF (save_bp_maps == 1) THEN
                temp_fn = trim(input_folder) // trim(fnames_baspred(3))//"_clipped.asc"
                print *, 'Writing clipped clay basin predictor file in input folder'
                !print *,
                call write_ascii_grid(temp_fn,clay_0_10,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
            bp_maps(:,:,3) = clay_0_10 !store all soils data in bp_maps for easy passing between functions
        END IF

        !Read in surface percentage organic matter - basin predictor 4
        IF (req_bp(4)) THEN
            print *, "Reading file ",trim(fnames_baspred(4))//".asc"
            call Read_clip_bp(trim(input_folder)//trim(fnames_baspred(4))//".asc",orgM_0_10, stats_HRU, HRU_map)
            IF (save_bp_maps == 1) THEN
                temp_fn =  trim(input_folder) // trim(fnames_baspred(4))//"_clipped.asc"
                print *, 'Writing clipped organic matter basin predictor file in input folder'
                !print *,
                call write_ascii_grid(temp_fn,orgm_0_10,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
            bp_maps(:,:,4) = orgm_0_10 !store all soils data in bp_maps for easy passing between functions
        END IF

!        !Read in form of ksat declines
!        IF (req_soils(5,1)) THEN
!            print *, "Readingfile ksat_bfit.asc"
!            call Read_clip_bp(trim(input_folder)//"INPUT/Ksat_bfit.asc",ksat_decline, stats_HRU, HRU_map)
!            IF (save_bp_maps == 1) THEN
!                temp_fn = trim(input_folder) // "OUTPUT/ksat_bfit.asc"
!                print *, 'Writing clipped ksat decline basin predictor file in output folder'
!                print *,
!                call write_ascii_grid(temp_fn,ksat_decline,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU, &
!                cellsize_HRU,nodata_HRU,6)
!            END IF
!            bp_maps(:,:,5) = ksat_decline !store all soils data in bp_maps for easy passing between functions
!        END IF

        !Read in ascii information describing SZM soils depth profiles - basin predictor 5
        IF (req_bp(5)) THEN
            print *, "Reading file ",trim(fnames_baspred(5))//".asc"
            temp_fn =(trim(input_folder) // trim(fnames_baspred(5))//".asc")
            call Read_clip_bp(trim(temp_fn),MUSIDs, stats_HRU, HRU_map)
            bp_maps(:,:,6) = MUSIDs !store all soils data in bp_maps for easy passing between functions
            bp_maps(:,:,5) = MUSIDs
            !print *, "Reading file ",trim(input_folder),trim(fnames_baspred(6))
            !temp_fn =  trim(input_folder) // trim(fnames_baspred(6))
            !call Read_clip_bp(trim(temp_fn),n_entries, stats_HRU, HRU_map)
            !bp_maps(:,:,7) = n_entries !store all soils data in bp_maps for easy passing between functions

            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped MUSIDs in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(5))//"_clipped.asc"
                call write_ascii_grid(temp_fn,MUSIDs,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
                !temp_fn =  trim(out_folder) // trim(fnames_baspred(6))//"_clipped.asc"
                !call write_ascii_grid(temp_fn,n_entries,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in tables of soils depth profiles - basin predictor 7-11
        IF (req_bp(7)) THEN
            i = 7
            !print *, "Readingfile ", trim(input_folder),trim(fnames_baspred(i)),".txt"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(i))//".txt"
            open (99, file = temp_fn, status = 'old')
            read(99,*) nrows, ncols, nmusids
            close(99)
            !print *, 'nrows = ',nrows,' ncols = ',ncols,' nmusids = ',nmusids
            allocate(soilmusiddata(nrows,ncols,5))

            DO i = 7,11
                print *, "Reading file ", trim(fnames_baspred(i))
                temp_fn =  trim(input_folder) // trim(fnames_baspred(i)) // ".txt"
                call read_numeric_list(temp_fn, ncols, 2,tempdata )
                soilmusiddata(:,:,i-6) = tempdata
                deallocate(tempdata)    !because this is allocated within read_numeric_list
            END DO
             if (save_bp_maps == 1) THEN
                !print *,
            end if
        END IF

        !Read in SRmax land cover map - basin predictor 12
        IF (req_bp(12)) THEN

            print *, "Reading file ",trim(fnames_baspred(12))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(12))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,12) = map !store all soils data in bp_maps for easy passing between functions

            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped Land Cover Map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(12))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
                !temp_fn =  trim(out_folder) // trim(fnames_baspred(6))//"_clipped.asc"
                !call write_ascii_grid(temp_fn,n_entries,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in SRmax rooting depth table  - basin predictor 13
        IF (req_bp(13)) THEN
            print *, "Reading file ", trim(fnames_baspred(13)),".txt"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(13))//".txt"
            open (99, file = temp_fn, status = 'old')
            read(99,*) nrows, ncols
            close(99)
            call read_numeric_list(temp_fn, 4, 3,tempdata )
            allocate(bp_root_depths(nrows,4))
            bp_root_depths= tempdata
            deallocate(tempdata)
            if (save_bp_maps == 1) THEN
                !print *,
            end if

        END IF


        !Read in Smax depth to bedrock map  - basin predictor 14
        IF (req_bp(14)) THEN
            print *, "Reading file ",trim(fnames_baspred(14))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(14))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,14) = map !store all soils data in bp_maps for easy passing between functions

            IF (save_bp_maps == 1) THEN

                print *, 'Writing clipped soil depth map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(14))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
                !temp_fn =  trim(out_folder) // trim(fnames_baspred(6))//"_clipped.asc"
                !call write_ascii_grid(temp_fn,n_entries,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in lnTo bulk density  - basin predictor 15
        IF (req_bp(15)) THEN
            print *, "Reading file ",trim(fnames_baspred(15))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(15))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,15) = map !store all soils data in bp_maps for easy passing between functions
            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped bulk density map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(15))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF


        !Read in ISorganic? map  - basin predictor 16
        IF (req_bp(16)) THEN
            print *, "Reading file ",trim(fnames_baspred(16))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(16))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,16) = map !store all soils data in bp_maps for easy passing between functions
            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped identify organic soils map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(16))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in porosity map  - basin predictor 17
        IF (req_bp(17)) THEN
            print *, "Reading file ",trim(fnames_baspred(17))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(17))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,17) = map !store all soils data in bp_maps for easy passing between functions
            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped porosity map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(17))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in hydrogeology map - basin predictor 18
        IF (req_bp(18)) THEN
            print *, "Reading file ",trim(fnames_baspred(18))//".asc"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(18))//".asc"
            call Read_clip_bp(trim(temp_fn),map, stats_HRU, HRU_map)
            bp_maps(:,:,18) = map !store all soils data in bp_maps for easy passing between functions
            IF (save_bp_maps == 1) THEN
                print *, 'Writing clipped hydrogeology map in input folder'
                !print *,
                temp_fn =  trim(input_folder) // trim(fnames_baspred(18))//"_clipped.asc"
                call write_ascii_grid(temp_fn,map,ncols_HRU,nrows_HRU,xllcorner_HRU,yllcorner_HRU,cellsize_HRU,nodata_HRU,6)
            END IF
        END IF

        !Read in B&Ks lookup table YZ Oct 2023
        IF (req_bp(19)) THEN
            print *, "Reading file ", trim(fnames_baspred(19)),".txt"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(19))//".txt"
            open (99, file = temp_fn, status = 'old')
            read(99,*) nrows, ncols
            close(99)
            call read_numeric_list(temp_fn, ncols, 3,tempdata )
            allocate(B_Ks_table(nrows,ncols))
            B_Ks_table= tempdata
            deallocate(tempdata)

        END IF

        !Read in T&Sy lookup table YZ Feb 2024
        IF (req_bp(20)) THEN
            print *, "Reading file ", trim(fnames_baspred(20)),".txt"
            temp_fn =  trim(input_folder) // trim(fnames_baspred(20))//".txt"
            open (100, file = temp_fn, status = 'old')
            read(100,*) nrows, ncols
            close(100)
            call read_numeric_list(temp_fn, ncols, 3,tempdata )
            allocate(T_Sy_table(nrows,ncols))
            T_Sy_table= tempdata

        END IF





    end subroutine read_req_soils


end module mpr_extract_BasPred


