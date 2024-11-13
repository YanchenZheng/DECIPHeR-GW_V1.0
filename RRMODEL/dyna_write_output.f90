module dyna_write_output
contains
!
!===============================================================
!  ROUTINE TO WRITE DYMOND OUTPUT
!===============================================================
!
    subroutine write_output(dyna_hru, &
        mcpar, &
        nac, &
        all_pm_names, &
        num_par_types, &
        num_rivers, &
        pf_data, &
        q, &
        s_full,&
        Ks,B,&
        GW_print_control,GWtable_out)

        use dyna_common_types

        ! Function to write output for dynamic topmodel

        implicit none


        double precision :: q(:,:)
        double precision :: s_full(:,:)
        double precision :: mcpar(:,:)
        integer :: i, ia
        integer :: num_par_types
        integer :: num_rivers
        integer :: nac,GW_print_control
        type(dyna_hru_type) :: dyna_hru(nac)
        type(pf_data_type) pf_data
        character(len=1024) :: line_format
        character(len=20), dimension(:) :: all_pm_names !names of all possible parameters
        character(len=20) :: pm_class = 'Param_class'
        character(len=20) :: Drz_class = 'Drz'
        double precision, allocatable,dimension(:) :: Ks
        double precision, allocatable,dimension(:) :: B
        character(len=20) :: Ks_class = 'Ks_new'
        character(len=20) :: B_class = 'B_new'
        double precision, allocatable, dimension(:,:) :: GWtable_out

        double precision, dimension(num_rivers) :: tot_reach_ppt, tot_reach_pet
        double precision, dimension(num_rivers) :: tot_reach_aet, tot_reach_frac
        !YZ GW 2024
        double precision, dimension(num_rivers) :: tot_reach_qof,tot_reach_qbf,tot_reach_quz


        !format
        !CHARACTER(LEN=200), PARAMETER :: FMT0 = "(T5,T5,T5,T5,T5,T5,T5)"
        CHARACTER(LEN=200), PARAMETER :: FMT0 = "(A16,A11,  A10,  A11 , A12,  A12,  A11,  A12,  A11,  A14,  A10)"
        !add Drz
        CHARACTER(LEN=200), PARAMETER :: FMT1 = "(I12,F10.4,F12.4,F10.4,F10.4,F13.4,F11.4,F10.4,F12.4,F14.7,F10.4)"

        !writing to the new .res file - Parameters first

        all_pm_names(1) = "szm_def"
        all_pm_names(2) = "lnto_def"
        all_pm_names(3) = "srmax_def"
        all_pm_names(4) = "srinit_def"
        all_pm_names(5) = "chv_def"
        all_pm_names(6) = "td_def"
        all_pm_names(7) = "smax_def"

        write(40, '(A)') '! PARAMETERS'
        write(40,FMT0) [pm_class,all_pm_names(1:7),Drz_class,Ks_class,B_class]

        !SRinit sets to the half of the SRmax
        do i = 1, num_par_types
            write(40,FMT1) i, [mcpar(i,1:3),mcpar(i,3)*0.5,mcpar(i,5:7),mcpar(i,10),Ks(i),B(i)] ! output Drz
        end do

        write(40, '(A)') ''

        ! writing to the .res file - summary stats

        tot_reach_pet = 0
        tot_reach_aet = 0
        tot_reach_frac = 0
        tot_reach_qof = 0
        tot_reach_quz = 0
        tot_reach_qbf = 0
        ! Calculate summary stats

        do ia = 1, nac

            tot_reach_ppt(dyna_hru(ia)%ipriv) = tot_reach_ppt(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_pptin*dyna_hru(ia)%ac)
            tot_reach_pet(dyna_hru(ia)%ipriv) = tot_reach_pet(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_pein*dyna_hru(ia)%ac)
            tot_reach_aet(dyna_hru(ia)%ipriv) = tot_reach_aet(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_aeout*dyna_hru(ia)%ac)

            tot_reach_frac(dyna_hru(ia)%ipriv) = tot_reach_frac(dyna_hru(ia)%ipriv) + dyna_hru(ia)%ac
 
            tot_reach_qof(dyna_hru(ia)%ipriv) = tot_reach_qof(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sumex*dyna_hru(ia)%ac)
            tot_reach_quz(dyna_hru(ia)%ipriv) = tot_reach_quz(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sumquz*dyna_hru(ia)%ac)
            tot_reach_qbf(dyna_hru(ia)%ipriv) = tot_reach_qbf(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sumqbf*dyna_hru(ia)%ac)

        end do

        do ia = 1, num_rivers

            tot_reach_ppt(ia) = (tot_reach_ppt(ia) / tot_reach_frac(ia))*1000
            tot_reach_pet(ia) = (tot_reach_pet(ia) / tot_reach_frac(ia))*1000
            tot_reach_aet(ia) = (tot_reach_aet(ia) / tot_reach_frac(ia))*1000
            
            tot_reach_qof(ia) = (tot_reach_qof(ia) / tot_reach_frac(ia))*1000
            tot_reach_quz(ia) = (tot_reach_quz(ia) / tot_reach_frac(ia))*1000
            tot_reach_qbf(ia) = (tot_reach_qbf(ia) / tot_reach_frac(ia))*1000

        end do

        write (line_format, '(A,I0, A)') '(',num_rivers,'(1X, f0.2), A)'

        write(40, '(A)') '! SUMMARY STATS PER RIVER REACH'

        write(40, '(A)') ''
        write(40, '(A)') '! PRECIP TOTALS (mm)'
        write(40, line_format) tot_reach_ppt

        write(40, '(A)') ''
        write(40, '(A)') '! PET TOTALS (mm)'
        write(40, line_format) tot_reach_pet

        write(40, '(A)') ''
        write(40, '(A)') '! ET TOTALS (mm)'
        write(40, line_format) tot_reach_aet

        !YZ GW 2024
        write(40, '(A)') ''
        write(40, '(A)') '! Recharge TOTALS (mm)'
        write(40, line_format) tot_reach_quz

        write(40, '(A)') ''
        write(40, '(A)') '! GW discharge TOTALS (mm)'
        write(40, line_format) tot_reach_qbf

        write (line_format, '(A,I0, A)') '(',num_rivers,'(1X, f0.5), A)'
        write(40, '(A)') ''
        write(40, '(A)') '! Overland Flow TOTALS (mm)'
        write(40, line_format) tot_reach_qof
        

        write(40, '(A)') ''
        write(40, '(A)') '! Reach fraction (%)'
        write(40, line_format) tot_reach_frac



        write (line_format, '(A,I0, A)') '(',size(q, 1),'(1X, f0.8), A)'

        do i = 1, size(q, 2)
            write(41, line_format) q(:, i)
        end do

        !YZ 2024 print out GW tables
        if (GW_print_control==1) then
            write (line_format, '(A,I0, A)') '(',size(GWtable_out, 1),'(1X, f0.3), A)'
            do i = 1, size(GWtable_out, 2)
                write(43, line_format) GWtable_out(:, i)
            end do
        end if 


        ! Write out storage file 
        if(pf_data%pf_enabled) then
            write (line_format, '(A,I0, A)') '(',size(s_full, 1),'(1X, f0.8), A)'

            do i = 1, size(s_full, 2)
                write(42, line_format) s_full(:, i)
            end do
        endif

    end subroutine

end module dyna_write_output
