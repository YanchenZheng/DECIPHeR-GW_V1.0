module dyna_common_types
    implicit none

    ! write out the lisflood output as in the old version
    integer, parameter :: ROUTE_MODE_NON_ROUTED = 1
    ! route to the reach outlets only (not enabled yet)
    integer, parameter :: ROUTE_MODE_TDH_REACH_ONLY = 2
    ! route to the reach outlets and down the river
    integer, parameter :: ROUTE_MODE_TDH_RIV_NET = 3

    interface checked_allocate
        module procedure checked_allocate_i
        module procedure checked_allocate_i_2
        module procedure checked_allocate_i_3

        module procedure checked_allocate_r
        module procedure checked_allocate_r_2
        module procedure checked_allocate_r_3
    end interface

    ! to be used as array of length nac
    type dyna_hru_type
        double precision :: ac

        integer :: ippt
        integer :: ipet
        integer :: ipar
        integer :: ims
        integer :: ipriv

        ! Model structure choices
        integer :: ims_rz
        integer :: ims_uz

        ! HRU flux weightings
        integer :: ntrans
        integer, allocatable, dimension(:) :: itrans
        double precision, allocatable, dimension(:) :: wtrans
        integer :: ntransriv

        double precision :: st
        double precision :: mslp
        !double precision :: q0
        double precision :: melev
        double precision :: numcells

        doubleprecision :: pex
        doubleprecision :: psrz
        !doubleprecision :: psuz
        doubleprecision :: quz
        doubleprecision :: p
        doubleprecision :: ep
        doubleprecision :: ea
        !doubleprecision :: uz
        doubleprecision :: ex
        ! doubleprecision :: exus
        ! doubleprecision :: exs

        !doubleprecision :: sd
        !doubleprecision :: suz
        doubleprecision :: srz
        !doubleprecision :: qin
        !doubleprecision :: qint1
        doubleprecision :: qbf
        !doubleprecision :: qbft1
        doubleprecision :: pc

       ! Summary stats for each HRU
        doubleprecision :: sum_pptin
        doubleprecision :: sum_pein
        doubleprecision :: sum_aeout
        doubleprecision :: sumex
        doubleprecision :: sumqbf
        doubleprecision :: sumquz
        ! doubleprecision :: sumexus
        !doubleprecision :: sdini
        !doubleprecision :: suzini
        doubleprecision :: srzini
         doubleprecision :: qbfini

        ! YZ GW model
        doubleprecision :: hgw,htopo

    end type dyna_hru_type

    type dyna_riv_type

        doubleprecision :: qb
        doubleprecision :: qof
        doubleprecision :: qout
        doubleprecision :: sumqb
        doubleprecision :: sumqof
        doubleprecision :: sumqout

    end type dyna_riv_type

    ! to be used to read in the initialisation file
    ! Gauge ID, number of downstreams, gauge area, downstream ID and area
    type route_init_reaches
        integer :: gauge_id
        integer :: num_ds
        double precision :: gauge_area
        integer, allocatable, dimension(:) :: ds_id
        double precision, allocatable, dimension(:) :: ds_area
    end type

    ! Added types from dyna_route_riv_tree_node.f90 and dyna_route_processing.f90

    type riv_tree_node
        !% 1=sea outlet
        !% 2=gauge
        !% 3=river reach increment
        integer :: node_type
        integer :: gauge_id
        !%list of cells directly connected upstream
        integer, allocatable, dimension(:) :: upstream_indexes
        integer upstream_count

        !%list of all upstream indexes from the entire tree
        integer, allocatable, dimension(:) :: upstream_tree_indexes
        integer :: upstream_tree_count

        integer :: row
        integer :: col

        logical :: enabled

        !% single index of the directly connected downstream
        integer :: downstream_index

        !% extra physical properties
        double precision :: point_h !% height of point on dem
        double precision :: point_dist !% distance of point on river to outlet

        double precision :: downstream_dist
        double precision :: downstream_delta_h
        double precision :: downstream_slope

        double precision :: reach_delay
        double precision :: total_downstream_delay
        double precision :: frac_sub_area
        double precision :: catch_area

    end type riv_tree_node

    type route_river_info_type
        logical:: is_enabled

        ! from flow_conn and point_con
        type(riv_tree_node), allocatable :: node_list(:)
        ! from river_data
        !% riv_id, area, dist, section_dist, slope, elevation
        double precision, allocatable, dimension(:,:) :: river_data


        ! only write output for node_type=NODE_TYPE_GAUGE
        integer :: flow_output_count
        integer, allocatable, dimension(:) :: flow_output_indexes
        character(64), allocatable, dimension(:) :: flow_output_headers

        character(1024) :: output_file_prefix

    end type

    ! routing time delay - processed from:
    ! * route_river_info_type
    ! * timestep
    ! * velocity params
    type route_time_delay_hist_type
        double precision :: timestep
        integer ::v_mode
        double precision, allocatable ::v_param(:)

        ! each row represents the histogram for the corresponding node_list item
        double precision, allocatable, dimension(:,:) :: route_hist_table
        ! each row represents the histogram index of flow output (:,1) for the corresponding node_list item
        ! (:,2) is the max index of the row histogram
        integer, allocatable, dimension(:,:) :: node_hist_indexes

        ! this is where all the magic happens
        ! flow is added to this table each timestep and is shifted one step left
        double precision, allocatable, dimension(:,:) :: flow_matrix

    end type

        ! --------------- pf_types ----------------------
        ! --------------- USER TYPES --------------------

    ! type 1 reservoir
    type pf_reservoir_param_type
        double precision :: s_max
        double precision :: s_0
        double precision :: demand
        double precision :: evap
        double precision :: env_min
        double precision :: catch_area
        double precision :: storage
        double precision :: frac_sub_area
    end type pf_reservoir_param_type

    type pf_reservoir_state_type
        double precision :: volume
    end type pf_reservoir_state_type

    type pf_policy_param_type
        double precision, allocatable, dimension(:) :: x0
        double precision, allocatable, dimension(:) :: x1
        double precision, allocatable, dimension(:) :: x2
        double precision, allocatable, dimension(:) :: x3

    end type pf_policy_param_type


    ! --------------- END USER TYPES --------------------

    type pf_location_type
        ! index into param_set array
        integer :: param_set_index

        ! index into the structure specific state array
        integer :: structure_state_index

        ! single output only - could be extended to support multiple
        integer :: n_outputs
        integer, allocatable, dimension(:) :: output_gauge
        integer, allocatable, dimension(:) :: output_node_index

        ! multiple inputs supported
        integer :: n_inputs
        integer, allocatable, dimension(:) :: input_gauge
        integer, allocatable, dimension(:) :: input_node_index
    end type pf_location_type

    type pf_param_set_type
        ! id from param set file
        integer :: id
        ! determines which type array to find this structure
        ! 1 = reservoir_param
        ! n = pf_example_list
        integer :: structure_type_id
        ! index in the structure specific param array
        integer :: structure_param_index
    end type pf_param_set_type



    type pf_data_type
        ! false if no processing functions used
        logical :: pf_enabled

        ! one per param set
        type(pf_param_set_type), allocatable, dimension(:) :: param_set

        ! one per processing function location
        type(pf_location_type), allocatable, dimension(:) :: location


        ! ------------ MEMORY FOR USER TYPES ---------------

        ! an entry for each reservoir param set
        type(pf_reservoir_param_type), allocatable, dimension(:) :: reservoir_param
        ! an entry for each instance of reservoir (reservoir at a location)
        type(pf_reservoir_state_type), allocatable, dimension(:) :: reservoir_state
        ! ---------------

        ! ---------------
        ! ------------ END MEMORY FOR USER TYPES ---------------
    end type pf_data_type

contains

    subroutine checked_allocate_r(ARR, d1)
        double precision, DIMENSION(:), ALLOCATABLE :: ARR
        integer :: d1
        integer :: AllocateStatus

        allocate(ARR (d1), stat = AllocateStatus)
        if (allocatestatus /= 0) stop "*** Not enough memory ***"

    end subroutine checked_allocate_r

    subroutine checked_allocate_i(ARR, d1)
        integer, DIMENSION(:), ALLOCATABLE :: ARR
        integer :: d1
        integer :: AllocateStatus

        allocate(ARR (d1), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i

    subroutine checked_allocate_r_2(ARR, d1, d2)
        double precision, DIMENSION(:,:), ALLOCATABLE :: ARR
        integer :: d1, d2
        integer :: AllocateStatus

        allocate(ARR (d1, d2), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_r_2

    subroutine checked_allocate_i_2(ARR, d1, d2)
        implicit none
        integer, DIMENSION(:,:), ALLOCATABLE :: ARR
        integer :: d1, d2
        integer :: AllocateStatus

        allocate(ARR (d1, d2), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i_2

    subroutine checked_allocate_r_3(ARR, d1, d2, d3)
        implicit none
        double precision, DIMENSION(:,:,:), ALLOCATABLE :: ARR
        integer :: d1, d2, d3
        integer :: AllocateStatus

        allocate(ARR (d1, d2, d3), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_r_3

    subroutine checked_allocate_i_3(ARR, d1, d2, d3)
        implicit none
        integer, DIMENSION(:,:,:), ALLOCATABLE :: ARR
        integer :: d1, d2, d3
        integer :: AllocateStatus

        allocate(ARR (d1, d2, d3), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i_3

end module dyna_common_types
