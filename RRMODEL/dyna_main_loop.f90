module dyna_main_loop
    implicit none
contains
    !
    !===========================================================
    !  MainLoop Subroutine to set parameters, initialise model and run Model
    !===========================================================
    !
    subroutine mainloop (all_pm_names, &
        nac, &
        nstep, &
        num_rivers, &
        acc, &
        dt, &
        dtt, &
        dyna_hru, &
        mcpar, &
        ntt, &
        node_to_flow_mapping, &
        num_par_types, &
        pe_step, &
        qobs_riv_step_start, &
        r_gau_step, &
        rivers, &
        route_riv, &
        route_tdh, &
        sum_ac_riv, &
        wt,&
        GW_data,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff,&
        stats_hru,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w,&
        GW_print_control,GW_maxiteration,tolerance)

        use dyna_common_types
        use dyna_param_setup
        use dyna_initialise_run
        use dyna_init_satzoneGW
        use dyna_topmod
        use dyna_write_output
        use dyna_route_processing
        use pf_common

        implicit none

        ! Argument Declares

        integer :: nac
        integer :: nstep
        integer :: num_rivers
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)

        ! Parameters
        character(len=20), dimension(7) :: all_pm_names !names of all possible parameters
        doubleprecision, allocatable, dimension(:,:) :: mcpar
        integer :: num_par_types
        integer :: ntt
        doubleprecision :: wt
        doubleprecision :: acc
        doubleprecision, allocatable, dimension(:) :: srmax
        doubleprecision, allocatable, dimension(:) :: chv
        doubleprecision, allocatable, dimension(:) :: chvdt
        doubleprecision, allocatable, dimension(:) :: lnto
        doubleprecision, allocatable, dimension(:) :: smax
        doubleprecision, allocatable, dimension(:) :: srinit
        doubleprecision, allocatable, dimension(:) :: szm
        doubleprecision, allocatable, dimension(:) :: t0dt
        doubleprecision, allocatable, dimension(:) :: td
        !YZ 2023
        doubleprecision, allocatable, dimension(:) :: Drz
        doubleprecision, allocatable, dimension(:) :: B 
        doubleprecision, allocatable, dimension(:) :: Ks

        ! Input and flow data
        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step
        double precision :: qobs_riv_step_start(num_rivers)

        ! River Data
        doubleprecision :: rivers(nac, num_rivers)
        double precision, dimension(:,:), allocatable :: q
        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh
        integer :: node_to_flow_mapping(:)
        doubleprecision :: sum_ac_riv (num_rivers)
        type(pf_data_type) :: pf_data
        double precision, dimension(:,:), allocatable :: s_full

        !YZ 2023 GW declares
        double precision, allocatable, dimension(:,:) :: GW_data,GW_Re,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff
        double precision, allocatable, dimension(:,:) :: GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w
        double precision, allocatable, dimension(:,:) :: GWtable_out
        double precision, dimension(4)  :: stats_HRU 
        integer::GW_print_control,GW_maxiteration
        double precision :: tolerance

        ! End declares

        !
        !=========================================================
        !  PUT THE PARAM VALUES INTO THE VARIABLE NAMES FOR THE MODEL TYPE
        !=========================================================
        !
        call param_setup (chv, &
            chvdt, &
            dt, &
            lnto, &
            num_par_types, &
            mcpar, &
            smax, &
            srinit, &
            srmax, &
            szm, &
            t0dt, &
            td,&
            B,&
            Ks,&
            Drz)
        
        !
        !=========================================================
        !  INITIALISE VARIABLES, ROUTING AND STORES
        !=========================================================
        !
        call initialise_run (dyna_hru, &
            dyna_riv, &
            nac, &
            num_rivers)

            ! todo
            ! initialise pf_param_sets_file and pf_location_file
            ! set pf_enabled = .true.
            pf_data%pf_enabled = .false.


            if(pf_data%pf_enabled) then
                print *, 'init processing functions'
                call pf_read_param_sets(pf_data)
                call pf_read_location(pf_data)
                call pf_node_init(route_riv, pf_data)
            else
                call pf_init_disabled(pf_data)
            endif


            call init_satzoneGW(nac, &
                nstep, &
                num_rivers, &
                chvdt, &
                dt, &
                dyna_hru, &
                mcpar, &
                pf_data, &
                q, &
                qobs_riv_step_start, &
                rivers, &
                route_riv, &
                route_tdh, &
                s_full, &
                smax, &
                srinit, &
                srmax, &
                sum_ac_riv, &
                szm, &
                t0dt)

        !
        !=========================================================
        !  CALL THE MAIN TOPMODEL STRUCTURE
        !=========================================================
        !
        call Topmod (nac, &
            nstep, &
            num_rivers, &
            acc, &
            dt, &
            dtt, &
            dyna_hru, &
            dyna_riv, &
            pf_data, &
            node_to_flow_mapping, &
            ntt, &
            pe_step, &
            q, &
            r_gau_step, &
            rivers, &
            route_riv, &
            route_tdh, &
            s_full, &
            smax, &
            srmax, &
            szm, &
            t0dt, &
            td, &
            wt,&
            GW_data,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff,&
            stats_hru,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w,&
            Drz,Ks,B,&
            GWtable_out,GW_maxiteration,tolerance)

        !
        !=========================================================
        !  WRITE OUTPUTS
        !=========================================================

        call write_output(dyna_hru, &
            mcpar, &
            nac, &
            all_pm_names, &
            num_par_types, &
            num_rivers, &
            pf_data, &
            q, &
            s_full,&
            Ks, B,&
            GW_print_control,GWtable_out)

        deallocate(mcpar)
        deallocate(Ks)
        deallocate(B)

        return
                                                                        
    end subroutine mainloop

end module dyna_main_loop
