!! Module containing subroutines for MPR.f90
!! This module contains all subroutines handling errors
!
!

!! Rosie Lane - 31st August 2017

module mpr_errors
contains


    ! **********************************************************************************************
    ! subroutine: check_files_exist
    ! **********************************************************************************************
    ! Checks that all the required files for MPR exist.
    ! Prints out appropriate error messages and stops the program if not.

    subroutine check_files_exist(MPRfolder,fpath_HRUmap,folder_input,folder_output,fname_control,fpath_PMSettings)

        implicit none

        !declare dummy variables
        character(len=1024), intent(in)     :: MPRfolder        !location of MPR folder containing INPUT and OUTPUT folders, and MPR_control file
        character(len=1024), intent(in)     :: fpath_HRUmap     !filepath pointing to HRU .asc file created in DTA
        character(len=1024), intent(in)     :: folder_input   !folder containing Basin Predictor files - INPUT folder within MPR main folder
        character(len=1024), intent(in)     :: folder_output   !output folder to store parameter input files - OUTPUT folder within MPR main folder
        character(len=1024), intent(in)     :: fname_control    !full path to control file
        character(len=1024), intent(in)     ::  fpath_PmSettings   !file containing parameter details

        !declare local variables
        logical, dimension(6)               :: LEXIST           !Logical - does file exist?


        ! initialise variables
        LEXIST(:) = .false.

        !check that each file required for MPR exists
        INQUIRE (FILE = MPRfolder, EXIST = LEXIST(1))
        INQUIRE (FILE = fpath_HRUmap, EXIST = LEXIST(2))
        INQUIRE (FILE = folder_input, EXIST = LEXIST(3))
        INQUIRE (FILE = folder_output, EXIST = LEXIST(4))
        INQUIRE (FILE = fname_control, EXIST = LEXIST(5))
        INQUIRE (FILE = fpath_PMSettings, EXIST = LEXIST(6))

        !set the appropriate error message and stop the program if any files dont exist
        !print *,

        IF (LEXIST(1) .EQV. .FALSE.) THEN
            print *, "ERROR: specified MPR folder does not exist!"
            print *, "Incorrect folder name given : ",trim(MPRfolder)
            !print *,
        END IF

        IF (LEXIST(2) .EQV. .FALSE.) THEN
            print *, "ERROR: specified fpath_HRUmap does not exist!"
            print *, "Incorrect HRU map folder given : ",trim(fpath_HRUmap)
            !print *,
        END IF

        IF (LEXIST(3) .EQV. .FALSE.) THEN
            print *, "ERROR: cannot find INPUT folder"
            print *, "The input folder MUST be located within the MPRfolder"
            print *, "Please create folder ", trim(folder_input)
            print *, "This must contain input Basin Predictors"
           ! print *,
        END IF

        IF (LEXIST(4) .EQV. .FALSE.) THEN
            print *, "ERROR: cannot find OUTPUT folder"
            print *, "The output folder MUST be located within the MPRfolder"
            print *, "Please create folder ", trim(folder_output)
            !print *,
        END IF

        IF (LEXIST(5) .EQV. .FALSE.) THEN
            print *, "ERROR: Control file does not exist!"
            print *, "The control file 'MPR_control.dat' MUST be saved in the MPRfolder/SETTINGS"
            print *, "Please create control file ", trim(fname_control)
            print *, "Format should be identical to previous control files."
            !print *,
        END IF

        IF (LEXIST(6) .EQV. .FALSE.) THEN
            print *, "ERROR: Parameter file does not exist!"
            print *, "The parameter file 'dynatop_parameters.dat' MUST be saved in the MPRfolder/SETTINGS"
            print *, "Please create parameter file ", trim(fpath_PMSettings)
            print *, "Format should be identical to previous parameter files."
            !print *,
        END IF

        IF (LEXIST(1) .eqv. .FALSE.) THEN
            ! stop the program because we are missing crucial files!!
            stop
        ELSEIF (LEXIST(2) .EQV. .FALSE.) THEN
            stop
        ELSEIF (LEXIST(3) .EQV. .FALSE.) THEN
            stop
        ELSEIF (LEXIST(4) .EQV. .FALSE.) THEN
            stop
        ELSEIF (LEXIST(5) .EQV. .FALSE.) THEN
            stop
        ELSEIF (LEXIST(6) .EQV. .FALSE.) THEN
            stop
        END IF







    end subroutine check_files_exist





end module mpr_errors


