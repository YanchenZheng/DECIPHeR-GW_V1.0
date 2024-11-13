module dyna_root_zoneGW
contains
!                                                                       
!===============================================================        
!  ROUTINES FOR CALCULATING THE ROOT ZONE FUNCTIONS                     
!  ROOT_ZONE1 = BASIC ROOT ZONE FORMULATION                             
!===============================================================        
!                                                                       
      subroutine root_zoneGW (nac, &
      dyna_hru, &
      ia, &
      it, &
      srmax,&
      Drz,&
      Ks,&
      B)

          use dyna_common_types
          implicit none

! Argument Declares
          integer :: nac
          type(dyna_hru_type) :: dyna_hru(nac)
          integer :: ia
          integer :: it
          integer :: s_it
          double precision, dimension(:) :: srmax
          !double precision, dimension(:) :: Td
          double precision, dimension(:) :: Drz

! Local Declares
          double precision :: storebal
          double precision :: storeend
          double precision :: storein
          double precision :: storeini
          double precision :: storeout

          !YZ add variables
          double precision,  dimension(:) :: Ks
          double precision,  dimension(:) :: B
          double precision :: s_quz,s_sumquz,sdt
          


! End declares

          dyna_hru(ia)%pex = 0.D0 
          storein = dyna_hru(ia)%p
          storeini = dyna_hru(ia)%srz 
          storeout = 0.D0 


          !  Explicit ROOT ZONE CALCULATIONS given P in this time step

          dyna_hru(ia)%srz = dyna_hru(ia)%srz + dyna_hru(ia)%p
          dyna_hru(ia)%sum_pptin = dyna_hru(ia)%sum_pptin + dyna_hru(ia)%p
          

          !YZ 2023
          !Three secnaries
          !check current water table
          !a. GW head below the bottom of root zone
          if (dyna_hru(ia)%hgw .lt. (dyna_hru(ia)%htopo-Drz(dyna_hru(ia)%ipar))) then

                !root zone is full
                if (dyna_hru(ia)%srz .gt.SRmax(dyna_hru(ia)%ipar) ) then

                    dyna_hru(ia)%pex = (dyna_hru(ia)%srz - SRmax(dyna_hru(ia)%ipar) )
                    dyna_hru(ia)%srz = SRmax(dyna_hru(ia)%ipar)
                    dyna_hru(ia)%ex = dyna_hru(ia)%pex
                    dyna_hru(ia)%sumex = dyna_hru(ia)%sumex + dyna_hru(ia)%ex

                    storeout = dyna_hru(ia)%ex 

                endif 

                !  CALCULATE EVAPOTRANSPIRATION FROM ROOT ZONE DEFICIT                  
                  dyna_hru(ia)%ea = 0

                  if (dyna_hru(ia)%ep .gt.0) then 

                      dyna_hru(ia)%ea = dyna_hru(ia)%ep * (dyna_hru(ia)%srz / SRmax(dyna_hru(ia)%ipar) )

                      if (dyna_hru(ia)%ea.gt.dyna_hru(ia)%ep ) dyna_hru(ia)%ea = dyna_hru(ia)%ep
                      if (dyna_hru(ia)%ea.ge.dyna_hru(ia)%srz ) dyna_hru(ia)%ea = dyna_hru(ia)%srz

                      dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%ea
                      dyna_hru(ia)%sum_pein = dyna_hru(ia)%sum_pein + dyna_hru(ia)%ep
                      dyna_hru(ia)%sum_aeout = dyna_hru(ia)%sum_aeout + dyna_hru(ia)%ea

                      storeout = storeout + dyna_hru(ia)%ea

                  endif 

                  !GW table below root zone
                  !new recharge equation
                 
                  dyna_hru(ia)%quz = Ks(dyna_hru(ia)%ipar)*(dyna_hru(ia)%srz / SRmax(dyna_hru(ia)%ipar))&
                  **((2+3*B(dyna_hru(ia)%ipar))/B(dyna_hru(ia)%ipar))
                  
                  
                ! !To avoid too big recharge
                if (dyna_hru(ia)%quz .ge. SRmax (dyna_hru(ia)%ipar)*0.05) then
                    !1 is one day, need to change if it runs in hour scale
                    !10 represents smaller time step
                        sdt=0.01 !(1/100)
                        s_quz=0
                        s_sumquz=0
                        do s_it=1,100
                            s_quz=Ks(dyna_hru(ia)%ipar)*(dyna_hru(ia)%srz / SRmax(dyna_hru(ia)%ipar))&
                                **((2+3*B(dyna_hru(ia)%ipar))/B(dyna_hru(ia)%ipar))
                            !quz can't be bigger than srz, otherwise cause negative srz
                             if (s_quz .ge.dyna_hru(ia)%srz) then
                                 s_quz=dyna_hru(ia)%srz
                             end if
                            dyna_hru(ia)%srz=dyna_hru(ia)%srz-s_quz*sdt
                            s_sumquz=s_sumquz+s_quz*sdt
                        end do
                        dyna_hru(ia)%quz =s_sumquz
                    
                else
                    
                    dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%quz

                end if
                                    
                ! !Old 
                ! if (dyna_hru(ia)%quz .ge. dyna_hru(ia)%srz) then
                !     dyna_hru(ia)%quz = dyna_hru(ia)%srz
                !     dyna_hru(ia)%srz = 0.D0
                ! else
                !     !print*,it, ia, 'quz > srz'
                !     dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%quz
                
                ! end if 

            storeout = storeout + dyna_hru(ia)%quz
                

          endif 
          
          !b. GW head reach & within the bottom of root zone
          !generate ex if root zone is full
           if ((dyna_hru(ia)%hgw .ge. (dyna_hru(ia)%htopo-Drz(dyna_hru(ia)%ipar))&
                  .AND.(dyna_hru(ia)%hgw<dyna_hru(ia)%htopo))) then
             
             !recharge=0
             dyna_hru(ia)%quz = 0.D0

               if (dyna_hru(ia)%srz .gt.SRmax (dyna_hru(ia)%ipar) ) then

                    dyna_hru(ia)%pex = (dyna_hru(ia)%srz - SRmax (dyna_hru(ia)%ipar) )
                    dyna_hru(ia)%srz = SRmax (dyna_hru(ia)%ipar)
                    dyna_hru(ia)%ex = dyna_hru(ia)%pex
                    dyna_hru(ia)%sumex = dyna_hru(ia)%sumex + dyna_hru(ia)%ex
                    
                    storeout = dyna_hru(ia)%ex 
              endif
               
             
                !  CALCULATE EVAPOTRANSPIRATION FROM ROOT ZONE DEFICIT                  
                dyna_hru(ia)%ea = 0

                if (dyna_hru(ia)%ep .gt.0) then 

                    dyna_hru(ia)%ea = dyna_hru(ia)%ep * (dyna_hru(ia)%srz / SRMAX (dyna_hru(ia)%ipar) )

                    if (dyna_hru(ia)%ea.gt.dyna_hru(ia)%ep ) dyna_hru(ia)%ea = dyna_hru(ia)%ep
                    if (dyna_hru(ia)%ea.ge.dyna_hru(ia)%srz ) dyna_hru(ia)%ea = dyna_hru(ia)%srz

                    dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%ea
                    dyna_hru(ia)%sum_pein = dyna_hru(ia)%sum_pein + dyna_hru(ia)%ep
                    dyna_hru(ia)%sum_aeout = dyna_hru(ia)%sum_aeout + dyna_hru(ia)%ea

                    storeout = storeout + dyna_hru(ia)%ea

                endif 

          endif

          !c. GW head reach topography
          !pex=p, root zone is not nesscery full
          if (dyna_hru(ia)%hgw .ge. dyna_hru(ia)%htopo) then

               !recharge=0
               dyna_hru(ia)%quz = 0.D0

               !rainfall goes to river channel directly
               dyna_hru(ia)%pex =  dyna_hru(ia)%p
               dyna_hru(ia)%ex = dyna_hru(ia)%pex
               dyna_hru(ia)%sumex = dyna_hru(ia)%sumex + dyna_hru(ia)%ex
               !2024 Jan rainfall directly goes to river
               dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%p

               storeout = dyna_hru(ia)%ex 

                 !  CALCULATE EVAPOTRANSPIRATION FROM ROOT ZONE DEFICIT                  
                dyna_hru(ia)%ea = 0

                if (dyna_hru(ia)%ep .gt.0) then 

                    dyna_hru(ia)%ea = dyna_hru(ia)%ep * (dyna_hru(ia)%srz / SRMAX (dyna_hru(ia)%ipar) )

                    if (dyna_hru(ia)%ea.gt.dyna_hru(ia)%ep ) dyna_hru(ia)%ea = dyna_hru(ia)%ep
                    if (dyna_hru(ia)%ea.ge.dyna_hru(ia)%srz ) dyna_hru(ia)%ea = dyna_hru(ia)%srz

                    dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%ea
                    dyna_hru(ia)%sum_pein = dyna_hru(ia)%sum_pein + dyna_hru(ia)%ep
                    dyna_hru(ia)%sum_aeout = dyna_hru(ia)%sum_aeout + dyna_hru(ia)%ea

                    storeout = storeout + dyna_hru(ia)%ea

                endif 

          endif


          ! Water balance calculations

          storeend = dyna_hru(ia)%srz 
          storebal = storeini + storein - storeout - storeend 
          if (storebal.gt.0.0000000001.or.storebal.lt. - 0.0000000001) then 
             print * , ' Storebal SRZ wrong ', it, ia, storebal 
             !pause
             STOP
          endif 

      end subroutine

end module
