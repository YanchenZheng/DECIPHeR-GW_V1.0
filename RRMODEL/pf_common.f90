module pf_common
! provides the common processing for processing functions
! this file shouldn't need to be modified to add additional processing functions
    use dyna_utility
    use pf_extend
    use dyna_route_riv_tree_node
    use dyna_common_types
    implicit none

contains

    subroutine pf_init_disabled(pf_data)
        implicit none
        type(pf_data_type), intent(inout) :: pf_data

        allocate(pf_data%location(0))
        allocate(pf_data%param_set(0))

    end subroutine

    subroutine pf_read_param_sets(pf_data)
        implicit none
        type(pf_data_type), intent(inout) :: pf_data
        integer :: ioerr

        ! read parameter set file
        open (99, file = 'pf_param_sets_fortran_clatworthy.txt', iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print*,'pf_read_param_sets could not open file: ' !, trim(param_set_filename)
            stop
        endif

        call pf_read_param_sets_fid(99, pf_data)

        close (99)
    end subroutine pf_read_param_sets

    subroutine pf_read_param_sets_fid(fid, pf_data)
    implicit none
        integer :: fid
        type(pf_data_type), intent(inout) :: pf_data
        character(len=1024) :: tmp, structure_type_name

        integer :: i
        integer :: n_param_sets, n_groups
        integer :: group_i, param_set_i
        integer :: structure_type_id, n_param_set_rows, n_param_set_cols
        integer :: id

        character(len=20), dimension(:), allocatable :: param_names
        character(len=64), dimension(:), allocatable :: param_values

        ! N_PARAM_SETS n
        read (fid, *) tmp, n_param_sets
        ! N_GROUPS n
        read (fid, *) tmp, n_groups

        allocate(pf_data%param_set(n_param_sets))

        param_set_i = 1

        do group_i=1,n_groups

            !read blank line before group
            !read (fid, *) tmp

            !print *, 'blank'
            !print *, tmp

            ! TYPE reservoir
            read (fid, *) tmp, structure_type_name
            structure_type_id = get_structure_type_id(trim(structure_type_name))

            print*,trim(structure_type_name), structure_type_id

            read (fid, *) tmp, n_param_set_rows
            read (fid, *) tmp, n_param_set_cols

            call pf_init_list(pf_data, structure_type_id, n_param_set_rows)

            allocate(param_names(n_param_set_cols))
            allocate(param_values(n_param_set_cols))

            read (fid, *) tmp, param_names

            ! Read in whole parameter values table
            do i = 1, n_param_set_rows

                read(fid, *) id, param_values(:)
                pf_data%param_set(param_set_i)%id = id
                pf_data%param_set(param_set_i)%structure_type_id = structure_type_id
                pf_data%param_set(param_set_i)%structure_param_index = i

                call pf_init_function(pf_data, structure_type_id, n_param_set_cols, i, param_names, param_values)

                param_set_i = param_set_i + 1
            end do

            deallocate(param_names)
            deallocate(param_values)
        enddo

        ! end read parameter set file
    end subroutine

    subroutine pf_read_location(pf_data)
        implicit none
        type(pf_data_type), intent(inout) :: pf_data
        integer :: ioerr
        ! read location set file
        open (99, file = 'pf_locations_fortran_clatworthy.txt', iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print*,'pf_read_location could not open file: ' !, trim(location_filename)
            stop
        endif

        call pf_read_location_fid(99, pf_data)

        close (99)

    end subroutine

    function pf_get_param_index(param_set_id, pf_data) result(param_set_index)
    implicit none
        integer :: param_set_id
        type(pf_data_type), intent(in) :: pf_data
        integer :: param_set_index
        integer :: i

        param_set_index = -1

        do i=1,size(pf_data%param_set)
            !print *, 'test', pf_data%param_set(j)%id, param_set_id
            if(pf_data%param_set(i)%id == param_set_id) then
                param_set_index = i
                exit
            endif
        end do

        if(param_set_index == -1) then
            print *, 'param id not found in param set file ', param_set_id
            stop
        endif

    end function

    subroutine pf_read_location_fid(fid, pf_data)
    implicit none
        integer :: fid
        type(pf_data_type), intent(inout) :: pf_data
        integer :: n_locs, i, j, n_in, n_out
        integer :: param_set_id, param_set_index, gauge
        integer :: structure_type_id

        integer :: state_ncols
        ! headers for state
        character(len=20), dimension(:), allocatable :: param_names
        ! values for state
        character(len=64), dimension(:), allocatable :: param_values

        character(len=1024) :: tmp

        ! N_LOCS n
        read (fid, *) tmp, n_locs
        if(are_equal(tmp,'NUM_LOCS').eqv..false.) then
            print *, 'Location file expected NUM_LOCS got ', trim(tmp)
            stop
        endif

        allocate (pf_data%location(n_locs))

        do i = 1, n_locs

            read (fid, *) tmp, param_set_id
            if(are_equal(tmp,'PF_PARAM_ID').eqv..false.) then
                print *, 'Location file expected N_OUT got ', trim(tmp)
                stop
            endif

            param_set_index = pf_get_param_index(param_set_id, pf_data)

            read (fid, *) tmp, state_ncols
            if(are_equal(tmp,'STATE_NCOLS').eqv..false.) then
                print *, 'Location file expected STATE_NCOLS got ', trim(tmp)
                stop
            endif

            if(state_ncols > 0) then
                allocate(param_names(state_ncols))
                allocate(param_values(state_ncols))
                read (fid, *) param_names
                read (fid, *) param_values
            else
                allocate(param_names(0))
                allocate(param_values(0))
            endif

            pf_data%location(i)%param_set_index = param_set_index

            structure_type_id = pf_data%param_set(param_set_index)%structure_type_id

            call pf_structure_init_state(pf_data, structure_type_id, i, state_ncols, param_names, param_values)

            deallocate(param_names)
            deallocate(param_values)

            read (fid, *) tmp, n_out
            !print *, tmp
            !print *, n_out
            if(are_equal(tmp,'N_OUT').eqv..false.) then
                print *, 'Location file expected N_OUT got ', trim(tmp)
                stop
            endif

            ! allocate list of outputs for location
            allocate(pf_data%location(i)%output_gauge(n_out))
            pf_data%location(i)%n_outputs = n_out

            do j=1,n_out
                read (fid, *) gauge
                pf_data%location(i)%output_gauge(j) = gauge
            end do

            read (fid, *) tmp, n_in
            if(are_equal(tmp,'N_IN').eqv..false.) then
                print *, 'Location file expected N_IN got ', trim(tmp)
                stop
            endif

            ! allocate list of inputs for location
            allocate(pf_data%location(i)%input_gauge(n_in))
            pf_data%location(i)%n_inputs = n_in

            do j=1,n_in
                read (fid, *) gauge
                pf_data%location(i)%input_gauge(j) = gauge
            end do
        end do

    end subroutine pf_read_location_fid

    subroutine pf_node_init(riv, pf_data)
        implicit none
        type(route_river_info_type) :: riv
        type(pf_data_type) :: pf_data

        integer :: i, j
        integer :: gauge_id, node_index


        ! loop over pf data locations
        ! record indexes

        do i = 1, size(pf_data%location)

            allocate(pf_data%location(i)%output_node_index(pf_data%location(i)%n_outputs))
            allocate(pf_data%location(i)%input_node_index(pf_data%location(i)%n_inputs))

            do j = 1, pf_data%location(i)%n_outputs
                gauge_id = pf_data%location(i)%output_gauge(j)
                node_index = get_node_index(riv%node_list, size(riv%node_list), gauge_id)
                pf_data%location(i)%output_node_index(j) = node_index
            end do

            do j = 1, pf_data%location(i)%n_inputs
                gauge_id = pf_data%location(i)%input_gauge(j)
                node_index = get_node_index(riv%node_list, size(riv%node_list), gauge_id)
                pf_data%location(i)%input_node_index(j) = node_index
            end do

        end do

    end subroutine

    subroutine pf_process(tdh, riv, pf_data)
        implicit none
        type(route_river_info_type) :: riv
        type(route_time_delay_hist_type) :: tdh
        type(pf_data_type) :: pf_data

        integer :: i, j, node_index
        double precision :: q_val, out_flow

        ! loop over pf data locations
        do i = 1, size(pf_data%location)
        ! Added this if to skip locations that aren't in the catchment being simulated 
            if (pf_data%location(i)%input_node_index(1) < 0) then 
            !print*, 'Skipping', pf_data%location(i)
            cycle !I want it to go to the next location in the loop 
            endif
            q_val = 0
            do j = 1, pf_data%location(i)%n_inputs

                node_index = pf_data%location(i)%input_node_index(j)

                ! record the flow from the histogram
                q_val = q_val + get_flow(riv, tdh, node_index, 1)

                ! zero the flow taken from histogram
                call set_flow(riv, tdh, node_index, 1, 0.0d0)
            end do

            out_flow = pf_process_location(pf_data, i, q_val)

            ! note - only one output supported (always 1)
            node_index = pf_data%location(i)%input_node_index(1)

            call set_flow(riv, tdh, node_index, 1, out_flow)
            end do

    end subroutine pf_process

end module pf_common
