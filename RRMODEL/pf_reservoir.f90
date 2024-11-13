module pf_reservoir
    use dyna_common_types
    implicit none

contains

    function pf_reservoir_param_init(n_cols, column_names, column_values) result(reservoir_param)
        implicit none
        integer n_cols
        character(len=20), dimension(n_cols) :: column_names
        character(len=64), dimension(n_cols) :: column_values
        type(pf_reservoir_param_type) reservoir_param
        integer col_index
        double precision tmp

        !reservoir_param%capacity = 1
        !reservoir_param%flow_release = 0

        do col_index=1,n_cols
            select case (trim(column_names(col_index)))
                case ('s_max')
                    read(column_values(col_index),*) tmp
                    reservoir_param%s_max = tmp
                case ('s_0')
                    read(column_values(col_index),*) tmp
                    reservoir_param%s_0 = tmp
                    reservoir_param%storage = tmp
                case ('demand')
                    read(column_values(col_index),*) tmp
                    reservoir_param%demand = tmp
                case ('evap')
                    read(column_values(col_index),*) tmp
                    reservoir_param%evap = tmp
                case ('env_min')
                    read(column_values(col_index),*) tmp
                    reservoir_param%env_min = tmp
                case ('catch_area')
                    read(column_values(col_index),*) tmp
                    reservoir_param%catch_area = tmp
                case ('frac_sub_area')
                    read(column_values(col_index),*) tmp
                    reservoir_param%frac_sub_area = tmp
                case default
                    print *, 'reservoir unknown column: ', trim(column_names(col_index))
            end select
        end do
    end function pf_reservoir_param_init

    function pf_reservoir_state_init(n_cols, param_names, param_values) result(reservoir_state)
        implicit none
        integer :: n_cols
        character(len=20), dimension(n_cols) :: param_names
        character(len=64), dimension(n_cols) :: param_values
        type(pf_reservoir_state_type) :: reservoir_state

        ! todo read columns

        reservoir_state%volume = 0

    end function pf_reservoir_state_init


!Simple alternative to test whether pf_functions are properly linking into routing 

    ! function pf_reservoir_process(q_value) result (out_flow)
    !     implicit none
    !     double precision :: q_value
    !     double precision :: out_flow

    !     q_value = q_value / 2
    !     out_flow = q_value 

    ! end function pf_reservoir_process

    function pf_reservoir_process(pf_data, param_index, q_value) result (out_flow)
        implicit none

        type(pf_data_type) :: pf_data
        type(pf_policy_param_type) :: param
        integer :: param_index, n_steps, i
        double precision :: q_value, q_value_convert
        double precision :: out_flow, s
        ! Need to declare and allocate a storage array somewhere! 
        double precision :: Qreg_max, Qreg_min, Qreg_mean, u_ref, s_ref_1, s_ref_2, u_0, u_3, s_step
        double precision :: s_max, s_0, demand, evap, env_min, catch_area, frac_sub_area
        double precision :: env, spill
        double precision, allocatable, dimension(:) :: s_frac, policy_rel

        s_max = pf_data%reservoir_param(param_index)%s_max
        s_0 = pf_data%reservoir_param(param_index)%s_0
        demand = pf_data%reservoir_param(param_index)%demand
        evap = pf_data%reservoir_param(param_index)%evap
        env_min = pf_data%reservoir_param(param_index)%env_min
        catch_area = pf_data%reservoir_param(param_index)%catch_area
        frac_sub_area = pf_data%reservoir_param(param_index)%frac_sub_area

        q_value_convert = q_value * catch_area * 0.001d0 / frac_sub_area
        ! Add fractional sub area to the input file perhaps? Or call it from elsewhere in the code? 
        
        s = pf_data%reservoir_param(param_index)%storage  
        !print *, s 

        if (env_min >= s + q_value_convert - evap) then
            env = env_min
        else 
            env = q_value_convert
        endif

        if (demand >= s  + q_value_convert - evap - env) then
            s = s  + q_value_convert - evap - env - demand
        else 
        ! For now this ignores demand but maybe it should just be reduced 
            s = s  + q_value_convert - evap - env 
        endif

        spill = MAXVAL((/0.d0 ,s + q_value_convert- env - evap - s_max/))

        pf_data%reservoir_param(param_index)%storage = s - spill
        out_flow = (env + spill) / 0.001d0 / catch_area * frac_sub_area  
            
            
    end function pf_reservoir_process

!     function pf_policy_function(param, args) result (u)
!         implicit none 
!         double precision :: args, u 
!         double precision, allocatable, dimension(:,:) :: x
!         double precision, allocatable, dimension(:) :: si, ui
!         integer :: i, points_dim
!         type(pf_policy_param_type) :: param

!         allocate(x(4,2))

!         !print *, 'Using policy function'

!         x(1,:) = param%x0
!         x(2,:) = param%x1
!         x(3,:) = param%x2
!         x(4,:) = param%x3

!         ! Should i not hard code this? Won't it always be 4?
!         points_dim  = SIZE(x,1)
        
!         ! ------------------------
!         ! Variables initialization
!         ! ------------------------
!         allocate(si(points_dim))
!         allocate(ui(points_dim))

!         do i=1, points_dim
!             si(i) = MINVAL(x(i:,1))
!             ui(i) = MINVAL(x(i:,2))
!         end do 
!         si(1) = 0.d0
!         si(4) = 1.d0
        
!         ! ------------------------
!         ! Reservoir release policy
!         ! -----------------------
!         u = interp(args, si, ui)

!     end function pf_policy_function

!     function interp(xi, X, Y) result (yi)
!     implicit none 
!     ! data has to be in a certain form, should add in break clauses for if si or ui has more than one dimension etc etc 
!     integer :: i, k
!     double precision, allocatable, dimension(:) :: X, Y
!     double precision :: xi, yi 
!     double precision :: Dy, Dx, m

!      do i = 1, SIZE(X) 
!         if (X(i) == xi) then 
!             yi = Y(i)
!             !print *, yi
!             exit
!         else if (X(i) > xi) then 
!             k = i - 1
!             Dy = Y(k+1) - Y(k)
!             Dx = X(k+1) - X(k)
!             m = Dy / Dx
!             yi = Y(k) + m * (xi - X(k)) 
!             !print *, yi
!             exit
!         endif
!     end do 

!     end function interp 

end module pf_reservoir
