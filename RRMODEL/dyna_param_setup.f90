module dyna_param_setup
contains
    !
    !===============================================================
    !  ROUTINES FOR SETTING UP THE PARAMETER VALUES FOR THE DIFFERENT
    !  HRUs
    !===============================================================
    !
    subroutine param_setup (chv, &
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

        use dyna_common_types

        implicit none

        ! Argument Declares

        double precision, dimension(:,:) :: mcpar

        double precision, allocatable, dimension(:) :: chv
        double precision, allocatable, dimension(:) :: chvdt
        double precision, allocatable, dimension(:):: lnto
        double precision, allocatable, dimension(:) :: smax
        double precision, allocatable, dimension(:) :: srinit
        double precision, allocatable, dimension(:) :: srmax
        double precision, allocatable, dimension(:) :: szm
        double precision, allocatable, dimension(:) :: t0dt
        double precision, allocatable, dimension(:) :: td
        !YZ 2023
        double precision, allocatable, dimension(:) :: Drz
        double precision, allocatable, dimension(:) :: B
        double precision, allocatable, dimension(:) :: Ks

        double precision :: dt
        integer :: num_par_types

        ! Local Declares
        integer :: i

        ! End declares

        allocate(chv(num_par_types))
        allocate(chvdt(num_par_types))
        allocate(lnto(num_par_types))
        allocate(smax(num_par_types))
        allocate(srinit(num_par_types))
        allocate(srmax(num_par_types))
        allocate(szm(num_par_types))
        allocate(t0dt(num_par_types))
        allocate(td(num_par_types))

        allocate(Drz(num_par_types))
        allocate(B(num_par_types))
        allocate(Ks(num_par_types))

        do i = 1, num_par_types

            szm (i) = mcpar(i, 1)
            lnto (i) = mcpar (i, 2)
            srmax (i) = mcpar (i, 3)
            !SRinit sets to the half of the SRmax
            srinit (i) = mcpar (i, 3)*0.5 !mcpar(i, 4)
            chv (i) = mcpar (i, 5)
            td (i) = mcpar (i, 6)
            smax (i) = mcpar (i, 7)
            chvdt (i) = CHV (i) * dt
            t0dt (i) = LnTo (i) + Log (dt)
            
            B (i) = mcpar (i, 8)
            Ks (i) = mcpar (i, 9)
            Drz (i) = mcpar (i, 10)

        end do

        end subroutine

    end module dyna_param_setup
