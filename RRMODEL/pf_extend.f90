module pf_extend
use dyna_common_types
use pf_reservoir
    implicit none

    contains

    function get_structure_type_id(structure_type_name) result(structure_type_id)
        implicit none
        character(*) :: structure_type_name
        integer :: structure_type_id

        select case (structure_type_name)
            case ('reservoir')
                structure_type_id = 1
            !---------------------
            !case ('example')
            !    get_group_type_id = n
            !---------------------
            case default
                print *, 'invalid structure name', trim(structure_type_name)
                stop
        end select
    end function get_structure_type_id

    subroutine pf_init_list(pf_data, structure_type_id, n_param_set_rows)
        implicit none
        type(pf_data_type) :: pf_data
        integer :: structure_type_id, n_param_set_rows

        select case (structure_type_id)
            !-------------
            case (1)
            allocate(pf_data%reservoir_param(n_param_set_rows))
            !-------------
            ! case n ! example
            ! allocate(pf_data.pf_example_list(n_param_set_rows))
        end select
    end subroutine pf_init_list

    subroutine pf_structure_init_state(pf_data, structure_type_id, location_index, n_cols, param_names, param_values)
        implicit none
        integer :: structure_type_id
        integer :: location_index
        integer :: n_cols
        character(len=20), dimension(n_cols) :: param_names
        character(len=64), dimension(n_cols) :: param_values
        type(pf_data_type) :: pf_data

        select case (structure_type_id)
            !-------------
            case (1)

            if(allocated(pf_data%reservoir_state).eqv..false.) then
                allocate(pf_data%reservoir_state(0))
            endif
            pf_data%location(location_index)%structure_state_index = size(pf_data%reservoir_state) + 1

            pf_data%reservoir_state = [pf_data%reservoir_state, &
                pf_reservoir_state_init(n_cols, param_names, param_values)]
            !-------------

        end select

    end subroutine pf_structure_init_state

    subroutine pf_init_function(pf_data, structure_type_id, n_cols, i, param_names, param_values)
        implicit none
        type(pf_data_type) :: pf_data
        integer :: structure_type_id, i, n_cols
        character(len=20), dimension(n_cols) :: param_names
        character(len=64), dimension(n_cols) :: param_values

        select case (structure_type_id)
            !-------------
            case (1)
                pf_data%reservoir_param(i) = pf_reservoir_param_init(n_cols, param_names, param_values)
            !-------------
        end select
    end subroutine pf_init_function

    function pf_process_location(pf_data, location_index, q_value) result(flow_out)
        implicit none
        type(pf_data_type) :: pf_data
        integer :: location_index
        double precision :: q_value
        double precision :: flow_out

        integer :: param_set_index
        integer :: structure_type_id

        integer :: structure_param_index
        integer :: structure_state_index

        structure_state_index = pf_data%location(location_index)%structure_state_index
        param_set_index = pf_data%location(location_index)%param_set_index

        structure_param_index = pf_data%param_set(param_set_index)%structure_param_index

        structure_type_id = pf_data%param_set(param_set_index)%structure_type_id

        flow_out = 0.0

        select case (structure_type_id)
            !-------------
            case (1)
                flow_out = pf_reservoir_process(pf_data, param_set_index, q_value)
            !-------------
        end select

    end function pf_process_location

end module pf_extend
