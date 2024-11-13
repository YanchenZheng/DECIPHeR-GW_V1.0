module dyna_topmod
contains
    !
    !=================================================================
    !  DynaTOP 1-2 (more modular structure for development)
    !=================================================================
    !  Keith Beven, Lausanne, 1997: Jim Freer Lancaster 12/1/00
    !  Updated by Gemma Coxon and Toby Dunne - update to Fortran 2003, made more modular
    !  Works with new routing routines and routes flows through a regional river netowrk
    !  add GW model YZ

    subroutine Topmod (nac, &
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
        GWtable_out,max_iteration,tolerance)
        
        use dyna_common_types
        use dyna_rain
        use dyna_dynamic_distGW
        use dyna_root_zoneGW
        use dyna_results_ac
        use dyna_river
        use dyna_results_tsGW
        use dyna_route_processing
        use GW_modelNEW
        !use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        implicit none

        ! Argument Declares
        integer :: nac
        integer :: nstep,GWbal_it
        integer :: num_rivers
        doubleprecision :: acc
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        integer :: ntt
        integer :: node_to_flow_mapping(:)

        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step

        double precision :: q(num_rivers, nstep)
        double precision :: rivers(nac, num_rivers)
        type(pf_data_type) :: pf_data
        double precision:: s_full(size(pf_data%location),nstep)

        ! Parameters
        doubleprecision, dimension(:) :: smax
        doubleprecision, dimension(:) :: srmax
        doubleprecision, dimension(:) :: szm
        doubleprecision, dimension(:) :: t0dt
        doubleprecision, dimension(:) :: td
        doubleprecision, dimension(:) :: Drz
        doubleprecision,  dimension(:) :: Ks
        doubleprecision,  dimension(:) :: B
        doubleprecision :: wt
        !YZ 2024
        double precision, allocatable, dimension(:,:) ::  GWtable_out
        


        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh

        ! Local Declares
        integer :: ia
        integer :: it,idGW
        integer :: max_iteration
        double precision :: tolerance,Solverresidual,GWbal_sumR,GWbal_sumQ, GWbal_sumQoutside
        double precision :: storeChange,GWbal_errors,GWbal_index
        doubleprecision :: sat_ac

        !YZ 2023 GW declares
        double precision, allocatable, dimension(:,:) :: GW_data,GW_Re,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff
        double precision, allocatable, dimension(:,:) :: GWbal_Hini,GWbal_Hend
        double precision, allocatable, dimension(:,:) :: GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w
        double precision, allocatable, dimension(:,:) :: GW_Re_temp
        double precision, allocatable, dimension(:) :: GW_hi_ID,GW_topo_ID,GW_runoff_ID
        integer :: id,GW_cellsize,i,j,Cell_Area
        double precision :: hru_re,hgw,htopo,hru_q,GWtablev,GWwbal,GWwbalHm
        double precision :: totalRunoff,sumGW_outsideCat,GW_outsideCat,sumGWwbal,sumGWtablev
        double precision, allocatable, dimension(:,:) :: it_p,it_pet,it_ea,it_pex,it_quz,it_srz
        double precision, allocatable, dimension(:,:) :: it_hgw,it_ex,it_qbf,it_stats
        double precision, dimension(4)  :: stats_HRU 
        character(900) :: filename
        character(8) :: tsname
        double precision, allocatable, dimension(:) :: TS_Recharge,TS_Runoff,TS_errors,TS_storechange,TS_index
        
        
        allocate(GW_Re(size(GW_data,1),size(GW_data,2)))
        allocate(GW_Re_temp(size(GWid,1),nac))
        allocate(GW_hi_ID(size(GWid,1)))
        allocate(GW_topo_ID(size(GWid,1)))
        allocate(GW_runoff_ID(size(GWid,1)))
        allocate(GWbal_Hini(size(GW_data,1),size(GW_data,2)))
        allocate(GWbal_Hend(size(GW_data,1),size(GW_data,2)))

        allocate(it_p(nstep,nac))
        allocate(it_pet(nstep,nac))
        allocate(it_ea(nstep,nac))
        allocate(it_pex(nstep,nac))
        allocate(it_quz(nstep,nac))
        allocate(it_srz(nstep,nac))
        allocate(it_hgw(nstep,nac))
        allocate(it_ex(nstep,nac))
        allocate(it_qbf(nstep,nac))
        allocate(it_stats(6,nac))

        allocate(GWtable_out(int(sum(GWid(:,7))), nstep))
        
        allocate(TS_Recharge(nstep))
        allocate(TS_Runoff(nstep))
        allocate(TS_errors(nstep))
        allocate(TS_index(nstep))
        allocate(TS_storechange(nstep))
        
        GW_Re_temp=0.D0
        GW_Re=0.D0
        GW_hi_ID=0.D0 !GW inital water tabel in GWid grid format
        GW_topo_ID=0.D0
        GW_runoff_ID=0.D0 !GW runoff in GWid grid format
        GW_cellsize=1000 !GW model spatial resolution 1km
        sumGW_outsideCat=0.D0
        totalRunoff=0.D0
        sumGWwbal=0.D0
        sumGWtablev=0.D0
        GWtable_out=0.D0
        Cell_Area=GW_cellsize*GW_cellsize
        GWbal_sumR=0.D0
        GWbal_sumQ=0.D0   
        GWbal_sumQoutside=0.D0
        GWbal_Hini=0.D0
        GWbal_Hend=0.D0
        storeChange=0.D0
        GWbal_errors=0.D0
        GWbal_index=0.D0

        TS_Recharge=0.D0
        TS_Runoff=0.D0
        TS_errors=0.D0
        TS_index=0.D0
        TS_storechange=0.D0
      

        ! End declares

        dtt = dt / dble (ntt)
        
        !YZ June 2023
        !Before start loop time step
        
        do id = 1, size(GWid,1)
            
            !GW water table in GWid format
            GW_hi_ID(id)=GW_hi(int(GWid(id,5)),int(GWid(id,6)))
            GW_topo_ID(id)=GW_data(int(GWid(id,5)),int(GWid(id,6)))

        end do 


        !==========================================================
        !  START LOOP ON TIME STEPS (MAIN MODEL LOOP)
        !=========================================================
        !For long-term GWbal: 20% timeseries data for warm up
        !calculate long-term GWbal for the rest 80% data
        GWbal_it= floor(nstep*0.2)
        
        do it = 1, nstep
            
            !  INITIALISE QB FOR EVERY RIVER REACH FOR THIS TIMESTEP
            !
            dyna_riv%qof = 0.D0
            dyna_riv%qb = 0.D0

            GW_Re_temp=0.D0
            GW_Re=0.D0
            GW_runoff=0.D0
            GW_runoff_ID=0.D0

            !
            !  CALL THE RAIN ROUTINE TO SPLIT UP THE RAIN INPUTS
            !  AND DISTRIBUTE THE EVAPORATION DATA
            !
            call rain (nac, &
                dyna_hru, &
                it, &
                pe_step, &
                r_gau_step)
            !
            !================================================================
            !
            !  FIRST PART OF THE DYNAMIC TOPMODEL FORMULATION
            !
            !  Distribute current QBF values to other elements
            !  Old time step QBF used to form new time step QIN values
            !
            !  CALL DYNAMIC_DIST() TO DO THESE CALCULATIONS
            !===========================================================
            !
             ! !UPdate dyna_riv(iriv)%qb?? with new hru%qbf
            ! !Update this below
            ! call dynamic_distGW (nac, &
            !     num_rivers, &
            !     dtt, &
            !     dyna_hru, &
            !     dyna_riv, &
            !     rivers)
            

            !=================================================================
            !  START LOOP ON Fractional Area Cells
            !=================================================================

            sat_ac = 0.D0
            hru_re=0.D0 !YZ track total volum of recharge from hru
            

            do ia = 1, nac
                

                dyna_hru(ia)%quz = 0.D0
                dyna_hru(ia)%ex = 0.D0
                

                !YZ 2023
                !initialize transfer GW table to HRU
                !calculate dyna_hru(ia)%hgw (m)
                !vol mapping
                dyna_hru(ia)%hgw = sum(GW_hi_ID*H_average_w(:,ia))
                dyna_hru(ia)%htopo = sum(GW_topo_ID*H_average_w(:,ia))
               
                !
                !=============================================================
                !  CALL THE APPROPRIATE ROOT ZONE FORMULATION (1 is standard)
                !  THIS CALCULATES THE MAIN P OR EVAP INPUTS/OUTPUTS TO SRZ
                !=============================================================
                if (dyna_hru(ia)%ims_rz.eq.1) then
                    
                    call root_zoneGW (nac, &
                        dyna_hru, &
                        ia, &
                        it, &
                        srmax,&
                        Drz,&
                        Ks,&
                        B)
                 
                end if
                !
                !=============================================================
                !  CALL THE APPROPRIATE SUZ ZONE FOMULATION (1 is standard)
                !  THIS CALCULATES QUZ AND CHANGES IN THE SUZ STORE
                !=============================================================
                
             
                !remove unsat zone
                !YZ June 2023
                !Within each time step, run GW model
                !for each HRU, the proportion of quz on every GW grid
                !GW_Re in V volume m3
                GW_Re_temp(:,ia)=dyna_hru(ia)%quz*Recharge_input_w(:,ia)
                !Check NaN values
                if (dyna_hru(ia)%quz/=dyna_hru(ia)%quz) then
                    print*, 'Check NaN values:',it, ia, dyna_hru(ia)%quz
                    stop
                end if 

                !count all HRU recharge total 
                !YZ mapping
                !hru_re=hru_re+dyna_hru(ia)%quz
                !volume mapping
                hru_re=hru_re+dyna_hru(ia)%quz*dyna_hru(ia)%numcells*((stats_hru(3)*stats_hru(3))/Cell_Area)
                dyna_hru(ia)%sumquz=dyna_hru(ia)%sumquz+dyna_hru(ia)%quz
                !
                !=============================================================
                !  CALL THE FRACTIONAL AREA RESULTS SUBROUTINE
                !=============================================================
                !only Qof (updat hru%ex to river%qof)
                call results_ac (nac, &
                    num_rivers, &
                    dyna_hru, &
                    dyna_riv, &
                    ia)

                !==========================================================
                !  END OF Fractional area LOOP
                !==========================================================
            
            !HRU loop finish
            end do

            !YZ 2023 GW
            !For each time step, Run GW model
            !Total all HRU contributions to every GW grid
            !Transfer GW recharge volumn m3 to m
            do id = 1, size(GWid,1)
                GW_Re(int(GWid(id,5)),int(GWid(id,6))) = real(sum(GW_Re_temp(id,:)),8)
            end do 
            
            !print*, 'Recharge from HRU: ', hru_re,' at TS ',it
            !Check the recharge total volumn same from HRU to GW model
            if ( abs(hru_re-sum(GW_Re))/hru_re>0.000001) then
                print*, 'Recharge from HRU: ', hru_re,' at TS ',it
                print*, "Error! Recharge from HRU is not equal to GW model "
                print*,sum(GW_Re)
                !pause
                stop
            end if 

            !if too big recharge
            !Assume the recharge can't be 1000 mm/d (i.e.,1 m/d) on average over the domain
            !Need to make sure (sumRecharge/A) < 1 m/d
            !A is the number of the cells (A, Grid domain)
             !if ( abs(hru_re-sum(GW_Re))>1) then
             if ( sum(GW_Re)/(size(GW_Re,1)*size(GW_Re,2))>1) then
                print*, 'Recharge from HRU: ', hru_re,' at TS ',it
                print*, "Numerical explosion! Recharge from HRU is too big "
                print*,sum(GW_Re)
                !pause !windows
                exit
            end if 
            

            !Long-term GW water balance check
            !Record the initial GW ini head at GWbal_it, i.e.,after warm up period
            if (it == GWbal_it) then !GWbal_it
               GWbal_Hini=GW_hi
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !For each time step
            !Run GW model
            !May 2024 YZ
            !max_iteration and tolerance are two parameter for GW solver
            !These two are user defined parameter
            !max_iteration default 1500, but can be set to any number, large number will increase the computation time
            !tolerance can be set to 1e-6 (strict,default) or 1e-4 (loose)
            !dt in GWmodel =1, run the model daily
            call GW_modelNEW_implicit(GW_data,size(GW_data,1),size(GW_data,2),&
                            1,GW_cellsize,GW_cellsize,GW_Sy,GW_Re,GW_hi,GW_T,GW_hNew,GW_runoff,&
                            max_iteration,tolerance,it,Solverresidual)
            
            !13th May 2024 YZ
            !Check GW solver calculation results
            !dot_product(r,r) can be NAN values
            !if (IEEE_IS_NAN(Solverresidual) ) then !for gfortran version 5
            if (is_nan(Solverresidual) ) then
                 print*,'At time step:',it
                 print*,'GW Solver cannot converge after',max_iteration,'iteration',' Residual:',Solverresidual
                 print*,'Skip this simulation'
                 exit
            elseif ( Solverresidual>=tolerance) then
                 print*,'At time step:',it
                 print*,'GW Solver cannot converge after',max_iteration,'iteration',' Residual:',Solverresidual
                 print*,'Skip this simulation'
                 exit
            end if
    
            
           !Long-term GW water balance check
           !Record the GW end head at final time step
            if (it == nstep) then !nstep
               GWbal_Hend=GW_hNew
            end if

            !GW table variation
            TS_storechange(it)=sum(GW_hNew*GW_Sy*Cell_Area)-sum(GW_hi*GW_Sy*Cell_Area)

            !Update GW head
            GW_hi=GW_hNew
            !Update GW
            idGW=1
            do id = 1, size(GWid,1)
                
                !GW runoff transfer into GWid format
                GW_Runoff_ID(id)=GW_Runoff(int(GWid(id,5)),int(GWid(id,6))) 
                !Nov 2023 update hru GW table
                GW_hi_ID(id)=GW_hnew(int(GWid(id,5)),int(GWid(id,6)))

                !YZ 2024 print out GW tables
                if (GWid(id,7)==1) then
                  GWtable_out(idGW,it)=GW_hi(int(GWid(id,8)),int(GWid(id,9)))
                  idGW=idGW+1
                end if

            end do 

            TS_Recharge(it)=sum(GW_Re*Cell_Area)
            TS_Runoff(it)=sum(GW_runoff*Cell_Area)
            TS_errors(it)=TS_Recharge(it)-TS_Runoff(it)-TS_storechange(it)
            TS_index(it)=TS_errors(it)/TS_Recharge(it)

            !Calculated accumulated recharge/runoff for the rest 80% time period
            !GWbal_it
            if (it .ge. GWbal_it) then
              GWbal_sumR = GWbal_sumR+sum(GW_Re*Cell_Area)
              GWbal_sumQ = GWbal_sumQ+sum(GW_runoff*Cell_Area)
              GWbal_sumQoutside = GWbal_sumQoutside+sum(GW_runoff*Cell_Area)-sum(GW_runoff_ID*Cell_Area)
               
           end if


            !!!!!!!!!!!!!!!!!!!
            !loop HRU, transfer GW_runoff to HRU%qbf
            hru_q=0.D0
            do ia= 1, nac
               !GW runoff transfer to HRU sat flow Qsat to river
               dyna_hru(ia)%qbf=sum(GW_Runoff_ID*Runoff_output_w(:,ia))
               !count all HRU runoff total volume
               hru_q=hru_q+dyna_hru(ia)%qbf*dyna_hru(ia)%numcells*((stats_hru(3)*stats_hru(3))/Cell_Area)
               dyna_hru(ia)%sumqbf=dyna_hru(ia)%sumqbf+dyna_hru(ia)%qbf

          
            end do
            
            
            !Check runoff with cats boudary total volumn same from GW to HRU model
            !print*, 'Runoff from GW model: ', sum(GW_Runoff_ID),' at TS ',it
            !/hru_q>0.000001
            if ( abs(hru_q-sum(GW_Runoff_ID))/hru_q>0.000001) then
                print*, 'Runoff from GW model: ', sum(GW_Runoff_ID),' at TS ',it
                print*, "Error! Runoff from GW model is not equal to HRU"
                print*,hru_q
                stop
            end if 
            
            !if too big runoff
            !Assume the runoff can't be 1000 mm/d (i.e.,1 m/d) on average over the domain
            !Need to make sure (sumrunoff/A) < 1 m/d
            !A is the number of the cells (A, Grid domain)
            !if ( abs(hru_q-sum(GW_Runoff_ID))>1) then
            if ( sum(GW_Runoff)/(size(GW_Runoff,1)*size(GW_Runoff,2))>1) then
                print*, 'Runoff from GW model: ', sum(GW_Runoff_ID),' at TS ',it
                print*, "Numerical explosion! Runoff from GW model is too big"
                print*,hru_q
                print*,sum(GW_Runoff)
                !pause !windows
                exit
            end if 

           
            !May 2024 long-term GWbal check
            if (it == nstep) then !nstep
               storeChange = sum(GWbal_Hend*GW_Sy*Cell_Area)-sum(GWbal_Hini*GW_Sy*Cell_Area)
               GWbal_errors=GWbal_sumR-GWbal_sumQ-storeChange
               GWbal_index=GWbal_errors/GWbal_sumR
               
               if (abs(GWbal_index) > 0.01) then
                   print *,'Warning:GWbal is larger than 1%',GWbal_index*100,'%'
               else 
                   print *,'GWbal is within 1%',GWbal_index*100,'%'
               end if

            end if

             
            ! !UPdate dyna_riv(iriv)%qb with new hru%qbf
            ! !should use dyna_dynamic_distGW.f90
            call dynamic_distGW (nac, &
                num_rivers, &
                dtt, &
                dyna_hru, &
                dyna_riv, &
                rivers)   


                !=============================================================
                !  CALL THE RIVER ROUTING SUBROUTINE
                !=============================================================
                call river (dyna_riv, &
                    nstep, &
                    num_rivers, &
                    it, &
                    node_to_flow_mapping, &
                    q, &
                    s_full, &
                    route_riv, &
                    route_tdh, &
                    pf_data)

                !=============================================================
                !  CALL THE TIMESTEP RESULTS SUBROUTINE - CHECKS WATER BALANCE
                !=============================================================

                    call results_tsGW (it,dtt, &
                        nac, &
                        num_rivers, &
                        dyna_hru, &
                        dyna_riv)
                !endif

                !==========================================================
                !    END OF TIME STEP LOOP
                !==========================================================
                    if (mod(it,5000) == 0) then
                        print*, it
                    end if 


                end do

                return

            end subroutine

            logical function is_nan(value)
                 double precision, intent(in) :: value
                 is_nan = (value /= value)
             end function is_nan

        end module dyna_topmod
