module GW_modelNEW

contains 
      
SUBROUTINE GW_modelNEW_implicit(GWtopo,NR,NC,dt,dx,dy,S,R,h,T,hNew,runoff,&
                             max_iteration,tolerance,timestep,Solverresidual)

  IMPLICIT NONE
  
  ! Declaration

  INTEGER, INTENT(IN)::                                         &
  ! Number of rows size(GWtopo,2)
  NR                                                            &
  ! Number of columns size(GWtopo,1)
  ,NC                                                           &
  ! Time increment (d)
  ,dt                                                           &
  ! topography resolution (m)
  ,dx,dy                                                           &
  ! GW solver maximum iteration number, user define
  ,max_iteration,timestep !Model time step, number of days

  ! Loop counters (work)
  INTEGER :: i,j,c,cc,row_idx,col_idx,left_row_idx,left_col_idx &
            ,right_row_idx , right_col_idx,up_row_idx,up_col_idx &
            ,down_row_idx,down_col_idx,x_dim,y_dim

  double precision, INTENT(IN)::                                            &
  ! Topography (m)
  GWtopo(NR,NC)                                                 &
  ! Specific yield  (-)
  ,S(NR,NC)                                                     &
  ! Recharge rate (ms-1)
  ,R(NR,NC)                                                     &
  ! Groundwater head above datum (m)
  ,h(NR,NC)                                                     &
  ! Transmissivity
  ,T(NR,NC)                               

  double precision,INTENT(OUT)::                                            &
  ! New head after solving 
  hNew(NR,NC) &   
  ! runoff 
  ,runoff(NR,NC)
  
  !define grid cells
  !xdim=size(GWtopo,2)*dx;
  !ydim=size(GWtopo,1)*dx;

  !% indexing for sparse matrix
  !Declare 

  INTEGER :: LHS_sparse_row(size(GWtopo,2)*size(GWtopo,1)*5),&
  LHS_sparse_col(size(GWtopo,2)*size(GWtopo,1)*5),&
  RHS_sparse_row(size(GWtopo,2)*size(GWtopo,1)),&
  RHS_sparse_col(size(GWtopo,2)*size(GWtopo,1))

  double precision :: LHS_sparse_val(size(GWtopo,2)*size(GWtopo,1)*5),&
  RHS_sparse_val(size(GWtopo,2)*size(GWtopo,1)),&
  rr(size(GWtopo,2)*size(GWtopo,1))

  double precision r_norm,tolerance,Solverresidual
  double precision x2(size(GWtopo,2)*size(GWtopo,1))
  double precision T_left,T_right,T_up,T_down
  double precision :: hNewN(NR,NC)
 
  LHS_sparse_row=0
  LHS_sparse_col=0
  RHS_sparse_row=0
  RHS_sparse_col=0

  LHS_sparse_val=0.0
  RHS_sparse_val=0.0

  T_left =0.0
  T_right =0.0
  T_up =0.0
  T_down =0.0

  x_dim = size(GWtopo,2);
  y_dim = size(GWtopo,1);

  c = 1;
  cc = 1;

  do i=1,size(GWtopo,2) !length(x) column
    do j=1,size(GWtopo,1) !length(y)
        
        
          !% set the cell indecies
            row_idx = j; col_idx = i;
            left_col_idx = i-1; right_col_idx = i+1;
            up_row_idx = j-1; down_row_idx = j+1;
           ! % handle boundary cells
            if (i==1) then
                left_col_idx = i;
            end if
            if (i==x_dim) then
                right_col_idx = x_dim;
            end if
            if (j==1) then
                up_row_idx = 1;
            end if 
            if (j==y_dim) then
                down_row_idx = y_dim;
            end if
        
        
        !  row_idx = j 
        !  col_idx = i
         
        ! ! assign boundaries to no flow
        ! if (i == 1) then
        !     left_col_idx = i;
        ! else
        !     left_col_idx = i - 1;
        ! end if
        
        ! if (i == x_dim) then
        !     right_col_idx = i;
        ! else
        !     right_col_idx = i + 1;
        ! end if
        
        ! if (j == 1) then
        !     up_row_idx = j;
        ! else
        !     up_row_idx = j - 1;
        ! end if 
        
        ! if (j == y_dim) then
        !     down_row_idx = j;
        ! else
        !     down_row_idx = j + 1;
        ! end if
        
        !calculate inter-cell transmissivity
        T_left = harmmean([T(row_idx,col_idx),T(row_idx,left_col_idx)])
        T_right = harmmean([T(row_idx,col_idx),T(row_idx,right_col_idx)])
        T_up = harmmean([T(row_idx,col_idx),T(up_row_idx,col_idx)])
        T_down = harmmean([T(row_idx,col_idx),T(down_row_idx,col_idx)])


       !% create sparse matrices with the coefficient
        !% LHS
        !% coefficient of j,i
        LHS_sparse_row(cc) = c
        LHS_sparse_col(cc) = (col_idx-1)*y_dim+row_idx
        !LHS_sparse_val(cc) = (T_left+T_right)/(dx*dx) + (T_up+T_down)/(dy*dy) + S(row_idx,col_idx)/dt
        LHS_sparse_val(cc) = (T_left+T_right) + (T_up+T_down) + (dx*dx)*S(row_idx,col_idx)/dt
        cc = cc + 1

        !% coefficient of j-1,i (up)
        LHS_sparse_row(cc) = c
        LHS_sparse_col(cc) = (col_idx-1)*y_dim+up_row_idx
        LHS_sparse_val(cc) = -T_up
        cc = cc + 1

        !% coefficient of j+1,i (down)
        LHS_sparse_row(cc) = c
        LHS_sparse_col(cc) = (col_idx-1)*y_dim+down_row_idx
        LHS_sparse_val(cc) = -T_down
        cc = cc + 1

        !% coefficient of j,i-1 (left)
        LHS_sparse_row(cc) = c
        LHS_sparse_col(cc) = (left_col_idx-1)*y_dim+row_idx
        LHS_sparse_val(cc) = -T_left
        cc = cc + 1

        !% coefficient of j,i+1 (right)
        LHS_sparse_row(cc) = c
        LHS_sparse_col(cc) = (right_col_idx-1)*y_dim+row_idx
        LHS_sparse_val(cc) = -T_right
        cc = cc + 1        

        !% RHS
        RHS_sparse_row(c) = c
        RHS_sparse_col(c) = 1
        RHS_sparse_val(c) = (dx*dx)*R(row_idx,col_idx)+((dx*dx)*S(row_idx,col_idx)/dt)*h(row_idx,col_idx)
                              
        c = c + 1

        END DO
    END DO
    
!  Call the CG routine.
      do i = 1, size(GWtopo,2)*size(GWtopo,1)
        x2(i) = 0.0
      end do

      call r8sp_cg (size(GWtopo,2)*size(GWtopo,1),size(GWtopo,2)*size(GWtopo,1)*5,&
      LHS_sparse_row,LHS_sparse_col,LHS_sparse_val,RHS_sparse_val,x2,&
      max_iteration,tolerance,timestep,Solverresidual)

!  Compute the residual.
      call r8sp_resid (size(GWtopo,2)*size(GWtopo,1),size(GWtopo,2)*size(GWtopo,1),size(GWtopo,2)*size(GWtopo,1)*5,&
      LHS_sparse_row,LHS_sparse_col,LHS_sparse_val,x2,RHS_sparse_val,rr)
      r_norm = r8vec_norm (size(GWtopo,2)*size(GWtopo,1), rr )

!  Report.
      !write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      
     hnew=0.0
     hnewN=0.0
     c=1
     do i=1,size(GWtopo,2)
       do j=1,size(GWtopo,1)

          hnew(j,i)=x2(c)
          hnewN(j,i)=x2(c)
          c=c+1

           if ((hNew(j,i)-GWtopo(j,i))*S(j,i)> 0 ) THEN
               runoff(j,i)=(hNew(j,i)-GWtopo(j,i))*S(j,i)
            else
               runoff(j,i)=0.0
            end if

            hNew(j,i)=hNew(j,i)-(runoff(j,i)/S(j,i))

         end do
      end do

END SUBROUTINE GW_modelNEW_implicit



 double precision function harmmean(X)

      IMPLICIT NONE
      ! Declaration
      double precision ::  X(2),Sum,InverseSum,TotalValid
      INTEGER :: Count,IO

      Sum        = 0.0
      InverseSum = 0.0
      TotalValid = 0
      Count      = size(X)

      do IO=1,Count
         if (X(IO)>0) THEN
            TotalValid = TotalValid + 1
            Sum        = Sum + X(IO)
            InverseSum = InverseSum + 1.0/X(IO)
         ELSE 
            print*, 'ERROR: Input has negative value, cannot give harmmean'
            EXIT
         end IF
      end do
      harmmean   = TotalValid / InverseSum

   END function harmmean


 subroutine r8sp_cg ( n, nz_num, row, col, a, b, x,&
                   max_iteration, tolerance,timestep,Solverresidual)

    ! cc R8SP_CG uses the conjugate gradient method on an R8SP system.
    ! c  Parameters:
    ! c    Input, integer N, the order of the matrix.
    ! c    N must be positive.
    ! c    Input, integer NZ_NUM, the number of nonzero elements in
    ! c    the matrix.
    ! c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
    ! c    column indices of the nonzero elements.
    ! c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
    ! c    Input, double precision B(N), the right hand side vector.
    ! c    Input/output, double precision X(N).
    ! c    On input, an estimate for the solution, which may be 0.
    ! c    On output, the approximate solution vector.
    !  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      implicit none

      integer n
      integer nz_num
      integer max_iteration,timestep

      double precision a(nz_num)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      double precision tolerance,Solverresidual

      integer col(nz_num)
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      integer row(nz_num)
      double precision rap
      double precision x(n)
        ! c  Initialize
        ! c    AP = A * x,
        ! c    R  = b - A * x,
        ! c    P  = b - A * x.
      call r8sp_mv ( n, n, nz_num, row, col, a, x, ap )

      do i = 1, n 
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
      !c  Do the N steps of the conjugate gradient method.
      !max_iteration maximum iteration
      do it = 1, max_iteration

        !c  Compute the matrix*vector product AP=A*P.
        call r8sp_mv ( n, n, nz_num, row, col, a, p, ap )
        ! c  Compute the dot products
        ! c    PAP = P*AP,
        ! c    PR  = P*R
        ! c  Set
        ! c    ALPHA = PR / PAP.
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
        ! c  Set
        ! c    X = X + ALPHA * P
        ! c    R = R - ALPHA * AP.
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
        ! Convergence
        if (dot_product(r,r)<tolerance) then
          !write(*,'(A,I4)') "Converged after iteration number ", it
          Solverresidual=dot_product(r,r)
          exit
        end if
        ! c  Compute the vector dot product
        ! c    RAP = R*AP
        ! c  Set
        ! c    BETA = - RAP / PAP.
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
        !  Update the perturbation vector
        !    P = R + BETA * P.
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do
        
         !13th MAY 2024 YZ
         !dot_product(r,r) can be NAN values
         !Save the residual of the last iteration
        if (it==max_iteration) then
            Solverresidual=dot_product(r,r)
        end if
        ! if (it==max_iteration .AND. IEEE_IS_NAN(dot_product(r,r)) ) then
        !    print*,"At time step:",timestep
        !    print*,"Can not Converge after iteration number",it," residual:", dot_product(r,r)
        !    Solverresidual=dot_product(r,r)
        ! elseif (it==max_iteration .AND. dot_product(r,r)>=tolerance) then
        !   print*,"At time step:",timestep
        !   print*,"Can not Converge after iteration number",it," residual:", dot_product(r,r)
        !   Solverresidual=dot_product(r,r)
        ! end if


      end do
      
     


      return
 end subroutine r8sp_cg

       
  subroutine r8sp_mv ( m, n, nz_num, row, col, a, x, b )

    ! cc R8SP_MV multiplies an R8SP matrix by an R8VEC.
    ! c  Parameters:
    ! c    Input, integer M, N, the number of rows and columns of the matrix.
    ! c    Input, integer NZ_NUM, the number of nonzero elements in the matrix.
    ! c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column indices of the nonzero elements.
    ! c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
    ! c    Input, double precision X(N), the vector to be multiplied by A.
    ! c    Output, double precision B(M), the product vector A*X.

      implicit none

      integer m
      integer n
      integer nz_num
      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
      end do

      do k = 1, nz_num

        i = row(k)
        j = col(k)
        b(i) = b(i) + a(k) * x(j)

      end do

      return

  end subroutine r8sp_mv

  subroutine r8sp_resid ( m, n, nz_num, row, col, a, x, b, r )

    ! cc R8SP_RESID computes the residual R = B-A*X for R8SP matrices.
    ! c  Parameters:
    ! c    Input, integer M, the number of rows of the matrix.
    ! c    M must be positive.
    ! c    Input, integer N, the number of columns of the matrix.
    ! c    N must be positive.
    ! c    Input, integer NZ_NUM, the number of nonzero elements in
    ! c    the matrix.
    ! c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
    ! c    column indices of the nonzero elements.
    ! c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
    ! c    Input, double precision X(N), the vector to be multiplied by A.
    ! c    Input, double precision B(M), the desired result A * x.
    ! c    Output, double precision R(M), the residual R = B - A * X.

      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      double precision r(m)
      integer row(nz_num)
      double precision x(n)

      call r8sp_mv ( m, n, nz_num, row, col, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
  end subroutine r8sp_resid
     

      double precision function r8vec_norm ( n, a )
         ! cc R8VEC_NORM returns the L2 norm of an R8VEC.
         ! c    An R8VEC is a vector of R8 values.
         ! c    The vector L2 norm is defined as:
         ! c      R8VEC_NORM = sqrt ( sum ( 1 .le. I .le. N ) A(I)^2 ).
         ! c  Parameters:
         ! c    Input, integer N, the number of entries in A.
         ! c    Input, double precision A(N), the vector whose L2 norm is desired.
         ! c    Output, double precision R8VEC_NORM, the L2 norm of A.

         implicit none

         integer n
         double precision a(n)
         integer i
         double precision value

         value = 0.0D+00
         do i = 1, n
         value = value + a(i) * a(i)
         end do
         value = sqrt ( value )

         r8vec_norm = value

         return
      end function r8vec_norm


      function r8vec_dot_product ( n, v1, v2 )

         ! cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
         ! c    An R8VEC is a vector of R8 values.
         ! c    In FORTRAN90, the system routine DOT_PRODUCT should be called
         ! c    directly.
         ! c  Parameters:
         ! c    Input, integer N, the dimension of the vectors.
         ! c    Input, double precision V1(N), V2(N), the vectors.
         ! c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
         implicit none

         integer n
         integer i
         double precision r8vec_dot_product
         double precision v1(n)
         double precision v2(n)
         double precision value

         value = 0.0D+00
         do i = 1, n
         value = value + v1(i) * v2(i)
         end do

         r8vec_dot_product = value

         return
      end function r8vec_dot_product

        




end module GW_modelNEW