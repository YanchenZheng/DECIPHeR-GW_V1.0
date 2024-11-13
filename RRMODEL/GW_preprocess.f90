module GW_preprocess

contains
      
      ! Yanchen 12th May 2023
    subroutine subset_GWdata(dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                    catch_mask, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata, &
                    buffer,GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)

            

            implicit none

            double precision, intent(in) :: dta_data(:,:), catch_mask(:,:)
            double precision, intent(in) :: xllcorner, yllcorner, nodata, cellsize
            double precision, intent(in) :: catch_xllcorner, catch_yllcorner, catch_nodata, catch_cellsize
            integer, intent(in) :: ncols, nrows, catch_ncols, catch_nrows, buffer

            ! Locals
            double precision :: GW_xllcorner, GW_yllcorner, GW_nodata
            double precision :: GW_cellsize, ix(2), iy(2)
            integer :: GW_ncols, GW_nrows, xlength, ylength
            double precision, allocatable :: GW_data(:,:)

            
            iy(1) = nrows-(floor((catch_yllcorner+10*catch_cellsize)/cellsize)-buffer)
            ix(1) = floor((catch_xllcorner+10*catch_cellsize)/cellsize)-buffer+1
            
            !right-left
            !10*50m buffer in hru_array file
            xlength=floor((catch_xllcorner+(catch_ncols-10)*catch_cellsize)/cellsize)-&
                 floor((catch_xllcorner+10*catch_cellsize)/cellsize)+1
            !top-bottom
            ylength=floor((catch_yllcorner+(catch_nrows-10)*catch_cellsize)/cellsize)-&
                 floor((catch_yllcorner+10*catch_cellsize)/cellsize)+1

            iy(2) = iy(1) - buffer - ylength -buffer +1
            ix(2) = ix(1) + buffer + xlength +buffer -1

            GW_ncols      = ix(2) - ix(1) +1
            GW_nrows      = iy(1) - iy(2) +1
            
            GW_xllcorner  = (floor((catch_xllcorner+10*catch_cellsize)/cellsize)-buffer)*cellsize
            GW_yllcorner  = (floor((catch_yllcorner+10*catch_cellsize)/cellsize)-buffer)*cellsize
            GW_cellsize   = cellsize
            GW_nodata     = nodata

            allocate (GW_data(GW_nrows,GW_ncols))

            GW_data = dta_data(int(iy(2)):int(iy(1)), int(ix(1)):int(ix(2)))


    end subroutine subset_GWdata

    
    ! Yanchen 24th May 2023
    subroutine calculate_SGW_weight(GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata, &
                    hru_mask, hru_ncols, hru_nrows, hru_xllcorner, hru_yllcorner, hru_cellsize, hru_nodata,&
                    GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w)
            
            use dyna_utility

            implicit none

            double precision, intent(in) :: GW_data(:,:)
            double precision, intent(in) :: GW_xllcorner, GW_yllcorner, GW_nodata, GW_cellsize
            double precision, intent(in) :: hru_xllcorner, hru_yllcorner,  hru_cellsize
            integer, intent(in) :: GW_ncols, GW_nrows, hru_ncols, hru_nrows
            double precision,allocatable :: hru_nodata
            double precision ::  hru_mask(:,:)
    
            

            ! Locals
            integer :: i, j,maxhru,t,n,count,ix,iy
            double precision, allocatable :: GWid_mask(:,:),Temp_GWgrid(:,:),GWhru(:,:),GWid(:,:)
            double precision, allocatable :: Recharge_input_w(:,:),Runoff_output_w(:,:),H_average_w(:,:)
            real(4), allocatable :: cyi(:),cxi(:),Ncyi(:),Ncxi(:)
            integer , dimension(:),allocatable :: guyi,guxi
            double precision :: GTxy(2),CTxy(2)
            double precision, allocatable :: subgw(:)
            double precision , allocatable :: Num_hru(:)

            hru_nodata=-9999

            !GWdata left top cell's bottom left corner xy
            GTxy(1)=GW_yllcorner+(GW_nrows-1)*GW_cellsize
            GTxy(2)=GW_xllcorner

            !hrumask left top cell's bottom left corner xy
            CTxy(1)=hru_yllcorner+(hru_nrows-1)*hru_cellsize
            CTxy(2)=hru_xllcorner

            allocate (cyi(hru_nrows))
            allocate (cxi(hru_ncols))
            allocate (Ncyi(hru_nrows))
            allocate (Ncxi(hru_ncols))
            cyi=0
            cxi=0
            Ncyi=0
            Ncxi=0
            allocate (guyi(hru_nrows))
            allocate (guxi(hru_ncols))
            guyi=0
            guxi=0

            j=1
            do i=1,hru_nrows
                cyi(i)=CTxy(1)/GW_cellsize-0.05*(i-1)
                if (i==1) then
                   guyi(j)=floor(cyi(i))
                   j=j+1
                else if (i>1 .AND. floor(cyi(i))/=guyi(j-1)) then
                   guyi(j)=floor(cyi(i))
                   j=j+1
                end if 
                Ncyi(i)=i
            enddo
            guyi=guyi(1:j-1)

            j=1
            do i=1,hru_ncols
                cxi(i)=CTxy(2)/GW_cellsize+0.05*(i-1)
                if (i==1) then
                   guxi(j)=floor(cxi(i))
                   j=j+1
                else if (i>1 .AND. floor(cxi(i))/=guxi(j-1)) then
                   guxi(j)=floor(cxi(i))
                   j=j+1
        
                end if 
                Ncxi(i)=i
            enddo
            guxi=guxi(1:j-1)
            
            maxhru=maxval(hru_mask)
            
            !calcualte weight through the loop
            allocate(subgw(maxhru))
            allocate(GWhru(gw_ncols*gw_nrows,maxhru))
            allocate(Runoff_output_w(gw_ncols*gw_nrows,maxhru))
            allocate(H_average_w(gw_ncols*gw_nrows,maxhru))
            allocate(GWid(gw_ncols*gw_nrows,10))
            GWid=0
            

            !convert nodata value to 0
            !GW calculate weight code need the nodata set to 0
            if (hru_nodata<0) then
                call set_nodata_value(hru_mask,hru_nodata,0.D0)
                hru_nodata = 0.D0
            end if 

            t=1
            do i=1,size(guyi)
              do j=1,size(guxi)
                
                ! print*,floor(sum(hru_mask(int(pack(Ncyi,cyi>=guyi(size(guyi)-i+1).AND.cyi<guyi(size(guyi)-i+1)+1)),&
                !         int(pack(Ncxi,cxi>=guxi(j).AND.cxi<guxi(j)+1)))))
                if (floor(sum(hru_mask(int(pack(Ncyi,cyi>=guyi(size(guyi)-i+1).AND.cyi<guyi(size(guyi)-i+1)+1)),&
                        int(pack(Ncxi,cxi>=guxi(j).AND.cxi<guxi(j)+1)))))/=0) then
                    

                    allocate(Temp_GWgrid(floor(GW_cellsize/hru_cellsize),floor(GW_cellsize/hru_cellsize)))
                    !GW grid
                    Temp_GWgrid=hru_mask(int(pack(Ncyi,cyi>=guyi(size(guyi)-i+1).AND.cyi<guyi(size(guyi)-i+1)+1)),&
                        int(pack(Ncxi,cxi>=guxi(j).AND.cxi<guxi(j)+1)))
                    
                    !Num_hru, which hrus are in this GW grid
                    
                    call find_unique_matrix(Temp_GWgrid,Num_hru)
                    !call get_unique_matrix(Temp_GWgrid,Num_hru) 
                    
                    subgw=0
                    do n=1,size(Num_hru)
                    count=0
                    do iy=1,size(Temp_GWgrid,1)
                        do ix=1,size(Temp_GWgrid,2)
                            if (Temp_GWgrid(iy,ix)==Num_hru(n)) then 
                                count=count+1
                            end if
                        end do
                    end do
                    if (int(Num_hru(n))>0) then
                        subgw(int(Num_hru(n)))=count
                    end if 
                    end do
                    GWhru(t,:)=subgw
                    !YZ mapping
                    !Runoff_output_w(t,:)=subgw/sum(subgw)
                    !volume mapping (first step not completed)
                    Runoff_output_w(t,:)=(GW_cellsize*GW_cellsize/(hru_cellsize*hru_cellsize))*subgw/sum(subgw)
                    H_average_w(t,:)=subgw
                    GWid(t,:)=[real(guyi(size(guyi)-i+1),8),real(guxi(j),8),real(t,8),sum(subgw),real(-99,8),real(-99,8)]
                            
                    t=t+1
                    deallocate(Temp_GWgrid)
                    deallocate(Num_hru)

                end if
               end do
            end do
           GWid=GWid(1:t-1,:)
           GWhru=GWhru(1:t-1,:)

           !weight for GW runoff output
           Runoff_output_w=Runoff_output_w(1:t-1,:)
           H_average_w=H_average_w(1:t-1,:)

           !weight for HRU recharge input
           allocate(Recharge_input_w(t-1,maxhru))
           do i=1,maxhru
              !YZ mapping
              !Recharge_input_w(:,i)=GWhru(:,i)/sum(GWhru(:,i))
              !volume mapping
               Recharge_input_w(:,i)=GWhru(:,i)/(GW_cellsize*GW_cellsize/(hru_cellsize*hru_cellsize))
               Runoff_output_w(:,i)=Runoff_output_w(:,i)/sum(GWhru(:,i))
               H_average_w(:,i)=H_average_w(:,i)/sum(GWhru(:,i))
           end do
           
           
           allocate(GWid_mask(GW_nrows,GW_ncols))
           GWid_mask=0
           do i=1,size(GWid,1)
              GWid_mask(int(GTxy(1)/GW_cellsize-GWid(i,1)+1),int(GWid(i,2)-GTxy(2)/GW_cellsize+1))=GWid(i,3)
              GWid(i,5)=real(GTxy(1)/GW_cellsize-GWid(i,1)+1,8)
              GWid(i,6)=real(GWid(i,2)-GTxy(2)/GW_cellsize+1,8)
           end do
           !allocate (GWid_mask(GW_nrows,GW_ncols))
            
            !set hru_nodata back to -999, just in case MPR needed
            if (hru_nodata /= -9999) then
                call set_nodata_value(hru_mask,hru_nodata,real(-9999,8))
                hru_nodata = -9999
            end if 


    end subroutine calculate_SGW_weight
    

    !23th May YZ calculate unique value of 2D matrix
    subroutine get_unique_matrix(matrix,unique_val) 
     
     implicit none
     double precision, intent(in) :: matrix(:,:)
     integer :: num_rows, num_columns
     integer :: i,j

     double precision, allocatable :: unique_array(:),array(:),unique_val(:)
     double precision :: temp

     num_rows=size(matrix,1)
     num_columns=size(matrix,2)
     
      ! Flatten the matrix into a 1D array
        array = reshape(matrix, [num_rows * num_columns])

        ! Sort the 1D array
        do i = 1, size(array)-1
            do j = 1, size(array)-i
                if (array(j) > array(j+1)) then
                    ! Swap elements
                    temp = array(j)
                    array(j) = array(j+1)
                    array(j+1) = temp
                end if
            end do
        end do
        
        allocate(unique_array(size(array)))
        allocate(unique_val(size(array)))
        unique_array=array
        
        ! Iterate through the sorted array and find unique values
        j=1
        if (unique_array(1)>0) then
           unique_val(j)=unique_array(1)
           j=j+1
        end if 
        
        do i = 2, size(unique_array)
            if (unique_array(i) /= unique_array(i-1)) then
                unique_val(j)=unique_array(i)
                j=j+1
            end if
        end do
        unique_val=pack(unique_val,unique_val>0)
    
     END subroutine get_unique_matrix



    !2024 March YZ calculate unique value of 2D matrix
    subroutine find_unique_matrix(matrix,unique_val) 
     
     implicit none
     double precision, intent(in) :: matrix(:,:)
     integer :: num_rows, num_columns
     integer :: i,j,count

     double precision  :: temp(5000)
     double precision, allocatable :: unique_val(:)

     num_rows=size(matrix,1)
     num_columns=size(matrix,2)
     temp=0

     count = 0
    do i = 1, num_rows
        do j = 1, num_columns
            if (matrix(i, j) /= 0.0) then
                count = count + 1
                temp(count) = matrix(i, j)
            endif
        end do
    end do
    
    allocate(unique_val(count))
    unique_val(1:count)=temp(1:count)

    END subroutine find_unique_matrix


end module GW_preprocess




