module dyna_file_open
contains
    !
    !===============================================================
    !  ROUTINES FOR OPENING THE OUTPUT FILES NEEDED PER SIMULATION RUN
    !===============================================================
    !
    subroutine file_open (i_mc, &
        out_dir_path, pf_data,GW_print_control)

        use dyna_common_types

        implicit none

        ! Argument Declares
        integer :: i_mc
        character(900) :: out_dir_path

        ! Local Declares
        character(900) :: filename
        character(20) :: mc_id_val

        type(pf_data_type) :: pf_data
        integer:: GW_print_control

        ! End declares

        ! Now open all the files needed
        ! First get the padded i_mc values into characters for
        ! all the filenames to open
        write (mc_id_val, "('mc_id_',I8.8)") i_mc

        ! HYDROL FILES
        ! 1) Then get the filename for the .res file
        filename = trim(out_dir_path) // &
            trim(mc_id_val) // &
            '.res'
        ! Then open this .res file
        open(unit=40,file=filename,status='unknown')

        ! 2) Then get the filename for the .flow file
        filename = trim(out_dir_path) // &
            trim(mc_id_val) // &
            '.flow'
        ! Then open this .flow file
        open(unit=41,file=filename,status='unknown')

        !YZ 2024 print out GW table
        !3) Write .gw file
        if (GW_print_control==1) then
            filename = trim(out_dir_path) // &
                trim(mc_id_val) // &
                '.gw'
            ! Then open this .gw file
            open(unit=43,file=filename,status='unknown')
        end if


        if(pf_data%pf_enabled) then
            ! 3) Then get the filename for the .stor file
            filename = trim(out_dir_path) // &
                trim(mc_id_val) // &
                '_stor' // &
                '.stor'
            ! Then open this .flow file
            open(unit=42,file=filename,status='unknown')
        endif

        return
    end subroutine file_open

end module dyna_file_open
