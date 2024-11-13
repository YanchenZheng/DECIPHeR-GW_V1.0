!! Module containing subroutines for MPR.f90
!! This module contains subroutines to interact with the control file
!
!   read_control_file_int   Read integer variables from the control file
!   read_param_file         reads parameter min/max values from file
!   file_open_err           Open file with error checking
!   read_mpr_commands       Reads the MPR folder and fpath_HRUmap from the command line

!!  Rosie Lane - 31st August 2017

module mpr_control
contains



    ! **********************************************************************************************
    ! subroutine: read_mpr_commands : Reads the MPR folder and fpath_HRUmap from the command line
    ! **********************************************************************************************
    ! This function reads the MPR folder and filepath to the HRU map from the command line
    ! It then makes the strings folder_input, folder_output and fname_control based on
    ! the location of the MPR main folder.

    subroutine read_mpr_commands(MPRfolder,fpath_HRUmap,folder_input,folder_output,fname_control)

        use dyna_utility

        implicit none

        !declare dummy variables
        character(len=1024), intent(out)        :: MPRfolder
        character(len=1024), intent(out)        :: fpath_HRUmap
        character(len=1024), intent(out)        :: folder_input
        character(len=1024), intent(out)        :: folder_output
        character(len=1024), intent(in)        :: fname_control

        !declare local variables
        character(len=1024)                     :: arg
        logical                                 :: input_is_valid
        integer                                 :: i

        ! initialise
        i = 0
        fpath_HRUmap = ""


        ! search for command line input:
        do
            CALL get_command_argument(i, arg)
            if(len_trim(arg) == 0) exit
            if (are_equal(arg, '-MPRfolder')) then
                CALL get_command_argument(i+1, MPRfolder)
                input_is_valid = .true.
            elseif (are_equal(arg, '-fpath_HRUmap')) then
                CALL get_command_argument(i+1, fpath_HRUmap)
            endif
            i = i + 1
        enddo

        if (len_trim(MPRfolder) == 0) then
            print *, '-MPRfolder not specified'
            input_is_valid = .false.
        elseif (len_trim(fpath_HRUmap) == 0) then
            print *, '-fpath_HRUmap not specified'
            input_is_valid = .false.
        endif

        if(input_is_valid .eqv. .false.) then
            print *, 'ERROR: MPR requires the following command options:'
            print *, "-MPRfolder <Full path to MPR folder>   e.g. '/home/rl1023/2017_07_20_MPR_setup/ '"
            print *, "-fpath_HRUmap <Full path to HRU map file>   "
            print *, "e.g. '/home/rl1023/2017_07_20_MPR_setup/INPUT/54057_classarray.asc'"
            stop
        endif

        folder_input = trim(MPRfolder) // "INPUT"
        folder_output = trim(MPRfolder) // "OUTPUT"

        print *, 'MPR folder given: ',trim(MPRfolder)
        !print *, 'MPR control file given: ', trim(fname_control)
        print *, 'HRU map location given: ',trim(fpath_HRUmap)
        !print *, 'Basin predictor files should be stored in: ',trim(folder_input)
        !print *, 'Parameter files will be stored in: ', trim(folder_output)

    end subroutine read_mpr_commands






    ! **********************************************************************************************
    ! subroutine: read_param_file : reads parameter min/max values from file
    ! **********************************************************************************************
    ! This function reads the parameter file.

    subroutine read_param_file(fpath,SZM_range,LnTo_range,SRmax_range,SRinit_range,CHV_range,Td_range,Smax_range,&
            B_range,Ks_range)

        implicit none

        !declare dummy variables
        character(len=1024), intent(in)                 :: fpath   !Full filename
        double precision, dimension(3), intent(out)     :: SZM_range
        double precision, dimension(3), intent(out)     :: LnTo_range   !min, best and max param values
        double precision, dimension(3), intent(out)     :: SRmax_range
        double precision, dimension(3), intent(out)     :: SRinit_range
        double precision, dimension(3), intent(out)     :: CHV_range
        double precision, dimension(3), intent(out)     :: Td_range
        double precision, dimension(3), intent(out)     :: Smax_range
        !YZ 2023
        double precision, dimension(3), intent(out)     :: B_range
        double precision, dimension(3), intent(out)     :: Ks_range

        !declare local variables
        character(len=1024)             :: temp_fn      !filename of parameter file
        character(len=64)               :: name         !name of parameter in the parameter file
        double precision, dimension(3)  :: dth1_range
        double precision, dimension(3)  :: fracdir_range


        !open the file
        temp_fn = fpath
        CALL file_open_err(temp_fn,2)

        !skip header lines
        Read(2,*)
        Read(2,*)

        !read parameters from each line
        Read(2,*) SZM_range(1), SZM_range(2), SZM_range(3), name
        IF (trim(name) /= "SZM") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : SZM expected on line 2"
            stop
        END IF

        Read(2,*) LnTo_range(1), LnTo_range(2), LnTo_range(3), name
        IF (trim(name) /= "LnTo") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : LnTo expected on line 3"
            stop
        END IF

        Read(2,*) SRmax_range(1), SRmax_range(2), SRmax_range(3), name
        IF (trim(name) /= "SRmax") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : SRmax expected on line 4"
            stop
        END IF

        Read(2,*) SRinit_range(1), SRinit_range(2), SRinit_range(3), name
        IF (trim(name) /= "SRinit") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : SRinit expected on line 5"
            stop
        END IF

        Read(2,*) CHV_range(1), CHV_range(2), CHV_range(3), name
        IF (trim(name) /= "CHV") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : CHV expected on line 6"
            stop
        END IF

        Read(2,*) Td_range(1), Td_range(2), Td_range(3), name
        IF (trim(name) /= "Td") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : Td expected on line 7"
            stop
        END IF

        Read(2,*) dth1_range(1), dth1_range(2), dth1_range(3), name
        IF (trim(name) /= "dth1") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : dth1 expected on line 8"
            stop
        END IF

        Read(2,*) Smax_range(1), Smax_range(2), Smax_range(3), name
        IF (trim(name) /= "Smax") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : Smax expected on line 9"
            stop
        END IF

        Read(2,*) fracdir_range(1), fracdir_range(2), fracdir_range(3), name
        IF (trim(name) /= "fracdir") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : fracdir expected on line 10"
            stop
        END IF

        Read(2,*) B_range(1), B_range(2), B_range(3), name
        IF (trim(name) /= "B") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : B expected on line 11"
            stop
        END IF

        Read(2,*) Ks_range(1), Ks_range(2), Ks_range(3), name
        IF (trim(name) /= "Ks") THEN
            print *, "ERROR IN PARAMETER RANGE FILE : Ks expected on line 12"
            stop
        END IF


        print *, "Parameter file successfully read-in. "
        !print *, "SZM min: ",SZM_range(1)," SZM max: ",SZM_range(3)
        !print *,


    end subroutine





    ! **********************************************************************************************
    ! subroutine: read_parammap_perturbations file
    ! **********************************************************************************************
    ! This function reads values from the parammap_perturbations file.
    ! These are the perturbations to the parameter maps if re-running a map for a set parameter.

    subroutine read_parammap_pert_file(fpath,parammap_pert)

            implicit none

        !declare dummy variables
        character(len=1024), intent(in)                 :: fpath   !Full filename
        double precision, dimension(7,4), intent(out)   :: parammap_pert !returned ranges for additions and multiplications to parameter map.

        !declare local variables
        character(len=64)               :: name         !name of parameter in the parameter file
        integer                         :: i

        !open the file
        CALL file_open_err(fpath,2)

        !skip header lines
        Read(2,*)
        Read(2,*)

        !read parameters from each line
        DO i = 1,7
            READ(2,*) parammap_pert(i,1), parammap_pert(i,2), parammap_pert(i,3), parammap_pert(i,4), name
            write(*,"(A28,A6,A3,F5.2,F5.2,F5.2,F5.2)") 'Set perturbation ranges for ',&
            trim(name),' = ',parammap_pert(i,:)
        END DO


    end subroutine



    ! **********************************************************************************************
    ! subroutine: read_control_file_int : Read integer variables from the control file
    ! **********************************************************************************************
    ! Given a variable name, this function checks the control file for that variable, and returns its
    ! associated integer vale.

    subroutine read_control_file_int(fname_control, varname, value,default_value)

        use dta_utility

        implicit none

        !declare dummy variables
        character(len=1024),intent(in)              :: fname_control       !filepath to control file
        character(len=*),intent(in)                 :: varname      !name of variable
        integer,intent(out)                         :: value        !variable values
        integer, intent(in)                         :: default_value  !default value to set if variable missing from control file

        !declare local variables
        character(len=1024)                         :: varname_read !variable name read from control file
        integer                                     :: value_read
        logical                                     :: var_unfound
        logical                                     :: var_missing

        !initialise local variables
        var_unfound = .true.
        var_missing = .true.

        CALL file_open_err(fname_control,2)

        READ(2,*) !read 1 header lines and ignore content

        DO WHILE (var_unfound )  !read all lines until finding the variable name of interest
            READ(2,*) varname_read, value_read
            if (are_equal(trim(varname), trim(varname_read))) then
                value = value_read
                print *, 'User defined ',trim(varname),' = ',value
                !print *, value ,' has been selected for ',trim(varname)
                var_unfound = .false.
                var_missing = .false.
            elseif (are_equal(trim(varname_read),'end')) then

                var_unfound = .false.
            endif


        END DO

        close(2)

        !If variable not declared in control file, then set to default value.
        IF (var_missing) THEN
            value = default_value
            !print *, 'IMPORTANT: variable ', trim(varname),' not specified in control file'
            print *, 'Default defn ',trim(varname),' = ',value
        ENDIF


    end subroutine read_control_file_int



    ! **********************************************************************************************
    ! subroutine: read_control_file_int : Read integer variables from a file
    ! **********************************************************************************************
    ! Given a variable name, this function checks the control file for that variable, and returns its
    ! associated integer vale.

    subroutine read_control_file_char(fname_control, varname, value,default_value)

        use dta_utility

        implicit none

        !declare dummy variables
        character(len=1024),intent(in)              :: fname_control       !filepath to control file
        character(len=*),intent(in)                 :: varname      !name of variable
        character(len=1024),intent(out)             :: value        !variable values
        character(len=1024), intent(in)             :: default_value  !default value to set if variable missing from control file

        !declare local variables
        character(len=1024)                         :: varname_read !variable name read from control file
        character(len=1024)                         :: value_read
        logical                                     :: var_unfound
        logical                                     :: var_missing

        !initialise local variables
        var_unfound = .true.
        var_missing = .true.

        CALL file_open_err(fname_control,2)

        READ(2,*) !read 1 header lines and ignore content

        DO WHILE (var_unfound )  !read all lines until finding the variable name of interest
            READ(2,*) varname_read, value_read
            if (are_equal(trim(varname), trim(varname_read))) then
                value = value_read
                print *, 'User defined ',trim(varname),' = ',trim(value)
                !print *, value ,' has been selected for ',trim(varname)
                var_unfound = .false.
                var_missing = .false.
            elseif (are_equal(trim(varname_read),'end')) then

                var_unfound = .false.
            endif
            !print *, 'Read varname: ',trim(varname_read)
            !print *, 'Read value: ',trim(value_read)
            !print *,

        END DO

        close(2)

        !If variable not declared in control file, then set to default value.
        IF (var_missing) THEN
            value = default_value
            !print *, 'IMPORTANT: variable ', trim(varname),' not specified in control file'
            print *, 'Default defn ',trim(varname),' = ',trim(value)
        ENDIF


    end subroutine read_control_file_char






    ! **********************************************************************************************
    ! subroutine: file_open_err : Open file with error checking
    ! **********************************************************************************************
    ! Checks that the file to be opened exists, and is not already in use.

    subroutine file_open_err(infile,unt)

        implicit none

        !declare dummy variables
        character(*),intent(in)             :: infile      ! filename
        integer,intent(in)                  :: unt         ! file unit
        character(len=1024)                 :: message     ! error message

        !declare local variables
        integer                             :: err
        logical                             :: xist        ! .TRUE. if the file exists
        logical                             :: xopn        ! .TRUE. if the file is already open

        ! initialize errors
        err=0; message="f-file_open/"

        ! check if the file exists
        inquire(file=trim(infile),exist=xist)
        if(.not.xist)then
            message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
            err = 1
            print *, message
            return
        endif

        ! check if the file is already open
        inquire(file=trim(infile),opened=xopn)
        if(xopn)then
            message=trim(message)//"FileIsAlreadyOpen[file='"//trim(infile)//"']"
            print *, message
            err = 1
            return
        endif

        ! open file
        open(unt,file=trim(infile),status="old",action="read",iostat=err)
        if(err/=0)then
            message=trim(message)//"OpenError['"//trim(infile)//"']"
            print *, message
            return
        endif

    end subroutine file_open_err



end module mpr_control


