module dyna_results_tsGW
contains
    !
    !===============================================================
    !  ROUTINES FOR CALCULATING THE TS DEPENDENT RESULTS
    !  FOR DYNAMIC VERSION
    !===============================================================
    !GW version
    subroutine results_tsGW (it,dtt, &
        nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv)

        use dyna_common_types
        implicit none
        ! Argument Declares
        integer :: nac,it
        integer :: num_rivers
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        double precision :: dtt

        ! Local Declares
        integer :: ia
        double precision :: sumqbf, sumsrz, sumqbfini
        !double precision :: sumqbf, sumqbfini, sumsd, sumsrz, sumsuz
        double precision :: sumallini,  sumsrzini !, sumsdini,sumsuzini
        double precision :: sumcatin, sumeout
        double precision :: sumex,sumquz!, sumexus
        double precision :: wbal,GWwbal,Fwbal,GWtablev,sumGW_outsideCat,GW_outsideCat
        double precision :: sumqb, sumqof, sumfout
        

        ! End declares

        !
        !  SUM UP THE HRU STORES AND CHECK THE WBAL PER TIME STEP
        !
        sumqbf = 0
        sumsrz = 0
        sumquz = 0
        sumallini = 0
        sumcatin =0
        sumex = 0
        sumqb = 0
        sumqof = 0
        sumfout = 0
        sumsrzini = 0
        sumeout = 0
        
                                                                        
        do ia = 1, nac

            sumqbf = sumqbf + dyna_hru(ia)%sumqbf * dyna_hru(ia)%ac
            sumex = sumex + dyna_hru(ia)%sumex * dyna_hru(ia)%ac
            
            
            sumsrz = sumsrz + dyna_hru(ia)%srz * dyna_hru(ia)%ac
            sumsrzini = sumsrzini + (dyna_hru(ia)%srzini*dyna_hru(ia)%ac)
            

            sumcatin = sumcatin + dyna_hru(ia)%sum_pptin * dyna_hru(ia)%ac
            sumeout = sumeout + dyna_hru(ia)%sum_aeout * dyna_hru(ia)%ac
            ! sumcatin = sumcatin + dyna_hru(ia)%p * dyna_hru(ia)%ac
            ! sumeout = sumeout + dyna_hru(ia)%ea * dyna_hru(ia)%ac
            
            
            sumquz = sumquz + dyna_hru(ia)%sumquz* dyna_hru(ia)%ac

        end do

        sumallini = sumsrzini !+ sumsuzini + sumsdini

        do ia = 1, num_rivers

            sumqb = sumqb + dyna_riv(ia)%sumqb
            sumqof = sumqof + dyna_riv(ia)%sumqof
            sumfout = sumfout + dyna_riv(ia)%sumqout
            ! sumqb = sumqb + dyna_riv(ia)%qb
            ! sumqof = sumqof + dyna_riv(ia)%qof
            ! sumfout = sumfout + dyna_riv(ia)%qout


        end do
        
        
        
        !wbal = sumcatin - sumeout - sumfout - (sumsd+sumsrz+sumsuz)+ sumallini + (sumqbfini- sumqbf)
        !wbal = sumallini + sumcatin - sumeout - sumfout - (sumsrz) + ((sumqbfini-sumqbf)*dtt)-sumGWtablev
        
        !print*,sumquz,sumqb,GWwbal
        !wbal = sumallini + sumcatin - sumeout +GWtablev*150/15.7-GW_outsideCat*122/15.7 - sumsrz- sumfout
        
        wbal = sumallini + sumcatin - sumeout  - sumsrz -sumqof-sumquz
        !wbal = sumallini + sumcatin - sumeout  - sumsrz -sumqof-(GWtablev*150/15.7+GW_outsideCat*122/15.7+sumqb)
        
        !print*,wbal
         !print*,'-------'
        
        if (abs(sumex + sumqbf - sumfout) > 0.0000001) then
            print * , ' Sum on outputs wrong in results_ts', sumex,      &
                 sumqof, sumqbf, sumfout
            print*,it
            print *, abs(sumex + sumqbf - sumfout)
            !pause
            stop
        endif

        if (abs(wbal).gt.0.0000001) then
            print *, 'wbal is wrong - check results'
            print *, it, wbal 
            print *, sumallini, sumcatin, sumeout
            print *, sumsrz, sumqof, sumquz
            !pause
            STOP
        end if

    end subroutine

end module dyna_results_tsGW
