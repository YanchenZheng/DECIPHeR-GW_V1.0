module GW_read
  
  
  contains

     subroutine GW_read_data(cats_name,hru_map,stats_HRU,folder_input,&
                 save_GW_maps,fnames_baspred,&
                 GW_data,GW_Sy,GW_T,GW_hi,GW_hNew,GW_runoff,GEO_index,&
                 GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w,&
                 GWgauge_meta,GW_print_control,buffer)
       
       use GW_preprocess
       use dyna_utility
      
       implicit none

       double precision, intent(in) :: hru_map(:,:)
       character(len=1024), intent(in)::folder_input
       character(len=6), intent(in)::cats_name
       double precision, dimension(4)  :: stats_HRU
       character(len=1024) :: GW_alldata_fn,catch_fn,GWmaskfile
       character(len=1024), dimension(27)   :: fnames_baspred 
       !clip data
       double precision, allocatable, dimension(:,:) :: GW_data,GW_Sy,GW_T
       !initialize
       double precision, allocatable, dimension(:,:) :: GW_hi,GEO_index
       !double precision, allocatable, dimension(:,:) :: GW_Re
       double precision, allocatable, dimension(:,:) :: GW_hNew,GW_runoff
       double precision, allocatable, dimension(:,:) ::GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w
       double precision, allocatable, dimension(:,:) :: dta_data,dta_Sy,dta_T,dta_head,dta_index
       integer :: ncols, nrows, GW_ncols, GW_nrows,I,J,NR,NC,buffer
       integer :: save_GW_maps
       double precision :: xllcorner, yllcorner, cellsize, nodata
       double precision :: GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata
       double precision,allocatable :: hru_nodata
       character(len=1024) :: tmp,filename
       integer :: num_cols,ngauge,GW_print_control,row,col
       integer :: exist_flag

       character(len=1024), allocatable, dimension(:) :: GWgauge_meta_cols
       double precision, allocatable, dimension(:,:) :: GWgauge_meta,GWtable_idxy
      
      !User defined parameter
      !buffer
      !!No.1
       !clip GW topography data
        if (save_GW_maps==1) then

          GW_alldata_fn=trim(folder_input)//'GB_TopoR1km.asc'
          call read_ascii_grid(GW_alldata_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
          
          call subset_GWdata(dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                          hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), stats_hru(4), &
                          buffer,GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
          
          !no data in topography set to 0, representing sea level
          if (GW_nodata<0) then
                  call set_nodata_value(GW_data,GW_nodata,0.D0)
                  GW_nodata = 0.D0
              end if 

          NR=size(GW_data,1)
          NC=size(GW_data,2)

          catch_fn = trim(folder_input)//trim(cats_name)//'_GWtopo_clipped.asc'
          
          call write_ascii_grid(catch_fn, GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, real(0,8), 8)
        
        elseif (save_GW_maps==0) then

            !Read in clipped topo data for faster runs
            GW_alldata_fn=trim(folder_input)//trim(fnames_baspred(21))//".asc"
            call read_ascii_grid(GW_alldata_fn, GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)

            NR=size(GW_data,1)
            NC=size(GW_data,2)

        end if

        
        !calculate GW input/output weight
        call calculate_SGW_weight(GW_data, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata, &
                    hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), hru_nodata,&
                    GWid_mask,GWid,GWhru,Recharge_input_w,Runoff_output_w,H_average_w)
        

        if (save_GW_maps==1) then
          !Only save the mask data when it runs for the first time
          GWmaskfile = trim(folder_input)//trim(cats_name)//'_GWmask_clipped.asc'
                  call write_ascii_grid(GWmaskfile, GWid_mask, GW_ncols, GW_nrows, &
                          GW_xllcorner, GW_yllcorner, GW_cellsize, real(0,8), 3)
        end if
        
        ! !Initial GW table (Topography-10m)
        !  allocate(GW_hi(NR,NC))
        !  do I=1,NR
        !   do J=1,NC
        !     if (GW_data(I,J)<0) then
        !         GW_hi(I,J)=GW_data(I,J)
        !       else 
        !         GW_hi(I,J)=GW_data(I,J)-5
        !     end if
        !   enddo
        ! ENDDO



        allocate(GW_hNew(NR,NC))
        allocate(GW_runoff(NR,NC))
        !allocate(GW_Re(NR,NC))
        GW_hNew=0
        GW_runoff=0
        !GW_Re=0



        !!No.2
        !clip GW Sy data
        if (save_GW_maps==1) then

            GW_alldata_fn=trim(folder_input)//'GB_Sy_Rm.asc'
            call read_ascii_grid(GW_alldata_fn, dta_Sy, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
            
            call subset_GWdata(dta_Sy, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                            hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), stats_hru(4), &
                            buffer,GW_Sy, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
            

            catch_fn = trim(folder_input)//trim(cats_name)//'_GWSy_clipped.asc'
            call write_ascii_grid(catch_fn, GW_Sy, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, nodata, 6)

        elseif (save_GW_maps==0) then

            !Read in clipped Sy data for faster runs
            GW_alldata_fn=trim(folder_input)//trim(fnames_baspred(23))//".asc"
            call read_ascii_grid(GW_alldata_fn, GW_Sy, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)

        end if

        !!No.3
        !clip GW T data
        if (save_GW_maps==1) then

            GW_alldata_fn=trim(folder_input)//'GB_T_Rm.asc'
            call read_ascii_grid(GW_alldata_fn, dta_T, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
            
            call subset_GWdata(dta_T, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                            hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), stats_hru(4), &
                            buffer,GW_T, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
            

            catch_fn = trim(folder_input)//trim(cats_name)//'_GWT_clipped.asc'
            call write_ascii_grid(catch_fn, GW_T, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, nodata, 3)
       
        elseif (save_GW_maps==0) then

            !Read in clipped T data for faster runs
            GW_alldata_fn=trim(folder_input)//trim(fnames_baspred(24))//".asc"
            call read_ascii_grid(GW_alldata_fn, GW_T, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)

        end if


        !!No.4
        !clip GW initial head data
        if (save_GW_maps==1) then

            GW_alldata_fn=trim(folder_input)//'GB_GW_inihead.asc'
            call read_ascii_grid(GW_alldata_fn, dta_head, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
            
            call subset_GWdata(dta_head, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                            hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), stats_hru(4), &
                            buffer,GW_hi, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
            
            !April 2024 replace -9999 to 0 in GW_hi
            do I=1,NR
              do J=1,NC
                if (GW_hi(I,J)==-9999) then
                    GW_hi(I,J)=0
                end if
              enddo
            ENDDO

            catch_fn = trim(folder_input)//trim(cats_name)//'_GWiniH_clipped.asc'
            call write_ascii_grid(catch_fn, GW_hi, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, 0.D0, 8)
          
          elseif (save_GW_maps==0) then

            !Read in clipped GW initial table data for faster runs
            GW_alldata_fn=trim(folder_input)//trim(fnames_baspred(25))//".asc"
            call read_ascii_grid(GW_alldata_fn, GW_hi, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)

         end if

        
         !!No.5
        !clip GW geo index map
        if (save_GW_maps==1) then

          GW_alldata_fn=trim(folder_input)//'GB_GW_GeoindexNEW101.asc'
          call read_ascii_grid(GW_alldata_fn, dta_index, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
          
          call subset_GWdata(dta_index, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                          hru_map, size(hru_map,2), size(hru_map,1), stats_hru(1), stats_hru(2), stats_hru(3), stats_hru(4), &
                          buffer,GEO_index, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
          

          catch_fn = trim(folder_input)//trim(cats_name)//'_GWGEOindex_clipped.asc'
          call write_ascii_grid(catch_fn, GEO_index, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, nodata, 3)
        
        elseif (save_GW_maps==0) then
           !Read in clipped Geo index data for faster runs
            GW_alldata_fn=trim(folder_input)//trim(fnames_baspred(26))//".asc"
            call read_ascii_grid(GW_alldata_fn, GEO_index, GW_ncols, GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, GW_nodata)
        end if

        
        !Feb 12th 2024 YZ
        !Step1: Read in GW OBS gauges XY locations
          filename = trim(folder_input)//'GWgauge_meta.txt'
          open (27, file = filename , status = 'unknown')
          ! First read in the meta file
          read (27, *) tmp
          read (27, *) tmp, GW_print_control
          read (27, *) tmp, ngauge
          read (27, *) tmp, num_cols
          read (27, *) tmp

          allocate(GWgauge_meta_cols(num_cols))
          allocate(GWgauge_meta(ngauge, num_cols))
          read (27, *) GWgauge_meta_cols

          ! Read in GW gauge Meta table
          !All avaliable GW gauges 
          do i = 1, ngauge
              read(27, *) GWgauge_meta(i,:)
          end do


          !Step2: Mapping GW gauges corresponding GW grids
          if (GW_print_control==1) then
                !only calcualtes the mapping when the GW resutls need to print out
                print*,'Total number of GW gauges:',size(GWgauge_meta,1)
                print*,'Print GW results:',GW_print_control
                allocate(GWtable_idxy(ngauge, num_cols))
                
                
                do i=1,size(GWgauge_meta,1)
                  !northing/easting
                  call NorthingEastingToRowCol(GWgauge_meta(i,3), GWgauge_meta(i,2), &
                          GW_nrows, GW_xllcorner, GW_yllcorner, GW_cellsize, row, col)
                    GWtable_idxy(i,1)=row
                    GWtable_idxy(i,2)=col
                    if (row .le. GW_nrows .AND. row .gt.0 .AND.col .le. GW_ncols.AND.col .gt.0) then
                        GWtable_idxy(i,3)=1
                        
                        !Only records the GW gauges that within the Catchments boundary
                        if (GWid_mask(row,col)>0) then
                            !print*,i
                            !indicator of GW gauge location
                            GWid(int(GWid_mask(row,col)),7)= 1
                            !record row in GW_mask matrix
                            GWid(int(GWid_mask(row,col)),8)=row
                            !record col in GW_mask matrix
                            GWid(int(GWid_mask(row,col)),9)=col
                            !index of GW gauge location
                            GWid(int(GWid_mask(row,col)),10)= GWgauge_meta(i,1)
                        end if
                      else
                        GWtable_idxy(i,3)=0
                    end if
                    
                end do

            !If there is no overlay GW gauge in the catchments
            !set GW_print_control to 0
            if (sum(GWid(:,7))==0) then         
              GW_print_control=0
            end if

            !print out GW obs gauge index
            if (GW_print_control==1) then

                filename = trim(folder_input)//trim(cats_name)//'_GWobs_gaugeID.txt'
                open (28, file = filename , status = 'unknown')
                write(28, '(A)') '! GW obs Gauge ID index'
                
                do i=1,size(GWid,1)
                    if (GWid(i,10)>0) then
                        write(28, '(I6)') int(GWid(i,10))
                    end if
                  end do
                  close(unit=28)

            end if

        end if

       


     end subroutine GW_read_data





end module GW_read