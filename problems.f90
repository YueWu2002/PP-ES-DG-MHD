module problems
	implicit none
	
	abstract interface
		function fun_type(x, y) result(u_pri)
			implicit none
			real, intent(in) :: x, y
			real, dimension(8) :: u_pri
		end function fun_type
	end interface
	
	! data used by 'init_by_fun'
	procedure(fun_type), pointer, public :: init_pri_fun => null()
	procedure(fun_type), pointer, public :: init_con_fun => null()
	
	!data used by 'BC_fixed'
	real, dimension(4), public :: u_pri_L, u_pri_R
	
	! data used by 'BC_by_fun'
	procedure(fun_type), pointer, public :: BC_pri_fun => null()
contains
	! implementation of initial and boundary conditions

	! initialize by the primitive variable function of space
	subroutine init_by_fun_pri(u_con, xx, yy, nx, ny, px, py)
		!! assign value to interior cells
		use MHD_PDE_m, only: MHD_pri2con
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l
		if (.not. associated(init_pri_fun)) error stop "init_by_fun_pri: required function pointer 'init_pri_fun' not initialized!"
		!$omp parallel do collapse(2) schedule(static)
		do j = 1, ny
			do i = 1, nx
				do l = 1, py+1
					do k = 1, px+1
						block
							! avoid generation of temporary array
							real :: u_pri_(8), u_con_(8)
							u_pri_ = init_pri_fun(xx(k,i), yy(l,j))
							u_con_ = MHD_pri2con(u_pri_)
							u_con(:,k,l,i,j) = u_con_
						end block
					end do
				end do
			end do
		end do
		!$omp end parallel do
	end subroutine init_by_fun_pri
	
	subroutine init_by_fun_con(u_con, xx, yy, nx, ny, px, py)
		!! assign value to interior cells
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l
		if (.not. associated(init_con_fun)) error stop "init_by_fun_con: required function pointer 'init_con_fun' not initialized!"
		!$omp parallel do collapse(2) schedule(static)
		do j = 1, ny
			do i = 1, nx
				do l = 1, py+1
					do k = 1, px+1
						u_con(:,k,l,i,j) = init_con_fun(xx(k,i), yy(l,j))
					end do
				end do
			end do
		end do
		!$omp end parallel do
	end subroutine init_by_fun_con
	
	pure function smooth_alfven_init_fun(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		
		real :: alpha, v_perp, B_perp, B_para, x_para, Bz, uz
		
		alpha = pi / 6.0
		x_para = x * cos(alpha) + y * sin(alpha)
		v_perp = 0.1 * sin(2*pi*x_para)
		B_perp = v_perp
		B_para = 1.0
		uz = 0.1 * cos(2*pi*x_para)
		Bz = uz
		
		u_pri(1) = 1.0
		u_pri(2) = v_perp * (-sin(alpha))
		u_pri(3) = v_perp * cos(alpha)
		u_pri(4) = uz
		u_pri(5) = B_para * cos(alpha) - B_perp * sin(alpha)
		u_pri(6) = B_para * sin(alpha) + B_perp * cos(alpha)
		u_pri(7) = Bz
		u_pri(8) = 0.1
	end function smooth_alfven_init_fun
	
	pure function KHI_init_pri(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		u_pri(1) = 1.0
		u_pri(2) = 5.0*( tanh(20*(y + 0.5)) - tanh(20*(y - 0.5)) - 1.0 )
		u_pri(3) = 0.25*sin(2*pi*x)*(exp(-100*(y+0.5)**2) - exp(-100*(y-0.5)**2))
		u_pri(4) = 0.0
		u_pri(5) = 1.0
		u_pri(6) = 0.0
		u_pri(7) = 0.0
		u_pri(8) = 50.0
	end function KHI_init_pri
	
	pure function Brio_Wu_init_fun(x, y) result(u_pri)
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		if (x<0.0) then
			u_pri(1) = 1.0
			u_pri(2) = 0.0
			u_pri(3) = 0.0
			u_pri(4) = 0.0
			u_pri(5) = 0.75
			u_pri(6) = 1.0
			u_pri(7) = 0.0
			u_pri(8) = 1.0
		else
			u_pri(1) = 0.125
			u_pri(2) = 0.0
			u_pri(3) = 0.0
			u_pri(4) = 0.0
			u_pri(5) = 0.75
			u_pri(6) = -1.0
			u_pri(7) = 0.0
			u_pri(8) = 0.1
		end if
	end function Brio_Wu_init_fun
	
	pure function smooth_1Dwave_init_fun(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		u_pri(1) = 1.0 + 0.2*sin(x)
		u_pri(2) = 1.0
		u_pri(3) = 0.0
		u_pri(4) = 0.0
		u_pri(5) = 0.5
		u_pri(6) = 1.0
		u_pri(7) = 1.5
		u_pri(8) = 2.0
	end function smooth_1Dwave_init_fun
	
	pure function OrszagTang_init_pri(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		u_pri(1) = 25.0/(36.0*pi)
		u_pri(2) = -sin(2*pi*y)
		u_pri(3) = sin(2*pi*x)
		u_pri(4) = 0.0
		u_pri(5) = -sin(2*pi*y) / sqrt(4.0*pi)
		u_pri(6) = sin(4*pi*x) / sqrt(4.0*pi)
		u_pri(7) = 0.0
		u_pri(8) = 5.0/(12.0*pi)
	end function OrszagTang_init_pri
	
	pure function Rotor_init_pri(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		real, parameter :: r0 = 0.1
		real, parameter :: r1 = 0.115
		real, parameter :: u0 = 2.0
		real :: r, f
		r = sqrt((x-0.5)**2 + (y-0.5)**2)
		f = (r1-r)/(r1-r0)
		if (r < r0) then
			u_pri(1) = 10.0
			u_pri(2) = -(u0/r0)*(y-0.5)
			u_pri(3) = (u0/r0)*(x-0.5)
		elseif (r < r1) then
			u_pri(1) = 1.0 + 9*f
			u_pri(2) = -(u0/r0)*(y-0.5) * f
			u_pri(3) = (u0/r0)*(x-0.5) * f
		else
			u_pri(1) = 1.0
			u_pri(2) = 0.0
			u_pri(3) = 0.0
		end if
		u_pri(4) = 0.0
		u_pri(5) = 5.0/sqrt(4*pi)
		u_pri(6) = 0.0
		u_pri(7) = 0.0
		u_pri(8) = 1.0
	end function Rotor_init_pri
	
	pure function Blast_init_pri(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		u_pri(1) = 1.0
		u_pri(2) = 0.0
		u_pri(3) = 0.0
		u_pri(4) = 0.0
		u_pri(5) = 100.0 / sqrt(4*pi)
		u_pri(6) = 0.0
		u_pri(7) = 0.0
		if ((x-0.5)**2 + (y-0.5)**2 < 0.01) then
			u_pri(8) = 1000.0
		else
			u_pri(8) = 0.1
		end if
	end function Blast_init_pri
	
	pure function cloud_shock_init_pri(x, y) result(u_pri)
		use utility, only: pi=>PI_math
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		real, parameter :: x0 = 0.8, y0 = 0.5, r0 = 0.15
		
		real, parameter :: rr = 3.868605300827511504664324644266228 ! compression ratio, -(5/6) (123 + 103 \[Pi]) + 5/6 Sqrt[16281 + 26298 \[Pi] + 10609 \[Pi]^2]
		real, parameter :: cf = 15.17658208897592333562014288162906 ! fast magnetosonic speed on right side, gamma must be 5.0/3.0
		real, parameter :: vl = 11.25357084620851400613837284504711 ! left speed
		real, parameter :: byl = 2.182626813584516145993203496356033 ! Byl
		real, parameter :: pl = 167.3451918203954512529698897850511 ! left pressure
		
		if (x<0.6) then
			! post-shock region
			u_pri = [3.86859, 0.0, 0.0, 0.0, 0.0, 2.1826182, -2.1826182, 167.345] ! original formulation
			! u_pri = [rr, 0.0, 0.0, 0.0, 0.0, byl, -byl, pl] ! refined
		else
			! pre-shock region
			u_pri = [1.0, -11.2536, 0.0, 0.0, 0.0, 0.56418958, 0.56418958, 1.0] ! original formulation
			! u_pri = [1.0, -vl, 0.0, 0.0, 0.0, 1.0/sqrt(pi), 1.0/sqrt(pi), 1.0] ! refined
			
			! cloud region
			if ((x-x0)**2 + (y-y0)**2 <= r0**2) then
				u_pri(1) = 10.0
			end if
		end if
	end function cloud_shock_init_pri
	
	pure function jet1_init_pri(x, y) result(u_pri)
		use MHD_PDE_m, only: gamma
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_pri
		real, parameter :: B2 = sqrt(200.0)
		real, parameter :: v2 = 800.0
		u_pri = [0.1*gamma, 0.0, 0.0, 0.0, 0.0, B2, 0.0, 1.0] ! ambient region
	end function jet1_init_pri
	
	pure function YeeSjogreenInitCon(x, y) result(u_con)
		implicit none
		real, intent(in) :: x, y
		real, dimension(8) :: u_con
		if (y>=0.0) then
			if (x>=0.0) then
				! 1st Quadrant
				u_con(1) = 0.9308
				u_con(2) = 1.4557
				u_con(3) = -0.4633
				u_con(4) = 0.0575
				u_con(5) = 0.3501
				u_con(6) = 0.9830
				u_con(7) = 0.3050
				u_con(8) = 5.0838
			else
				! 2st Quadrant
				u_con(1) = 1.0304
				u_con(2) = 1.5774
				u_con(3) = -1.0455
				u_con(4) = -0.1016
				u_con(5) = 0.3501
				u_con(6) = 0.5078
				u_con(7) = 0.1576
				u_con(8) = 5.7813
			end if
		else
			if (x<=0.0) then
				! 3rd Quadrant
				u_con(1) = 1.0000
				u_con(2) = 1.7500
				u_con(3) = -1.0000
				u_con(4) = 0.0000
				u_con(5) = 0.5642
				u_con(6) = 0.5078
				u_con(7) = 0.2539
				u_con(8) = 6.0000
			else
				! 4th Quadrant
				u_con(1) = 1.8887
				u_con(2) = 0.2334
				u_con(3) = -1.7422
				u_con(4) = 0.0733
				u_con(5) = 0.5642
				u_con(6) = 0.9830
				u_con(7) = 0.4915
				u_con(8) = 12.999
			end if
		end if
	end function YeeSjogreenInitCon
	
	!-------------------------------------------------------------------------!
	
	subroutine BC_periodic(u_con, xx, yy, nx, ny, px, py)
		!! assign value to ghost cells
		use MHD_PDE_m, only: MHD_con2pri
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j
		! !$omp parallel workshare ! excecute one-by-one
		! u_con(:,:,:,0,1:ny) = u_con(:,:,:,nx,1:ny) ! assign left ghost cells
		! u_con(:,:,:,nx+1,1:ny) = u_con(:,:,:,1,1:ny) ! assign right ghost cells
		! u_con(:,:,:,1:nx,0) = u_con(:,:,:,1:nx,ny) ! assign bottom ghost cells
		! u_con(:,:,:,1:nx,ny+1) = u_con(:,:,:,1:nx,1) ! assign top ghost cells
		! !$omp end parallel workshare
		!$omp parallel
		!$omp do schedule(static)
		do j = 1, ny
			u_con(:,:,:,0,j) = u_con(:,:,:,nx,j) ! assign left ghost cells
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do j = 1, ny
			u_con(:,:,:,nx+1,j) = u_con(:,:,:,1,j) ! assign right ghost cells
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			u_con(:,:,:,i,0) = u_con(:,:,:,i,ny) ! assign bottom ghost cells
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			u_con(:,:,:,i,ny+1) = u_con(:,:,:,i,1) ! assign top ghost cells
		end do
		!$omp end do nowait
		!$omp end parallel
	end subroutine BC_periodic
	
	subroutine BC_BrioWu(u_con, xx, yy, nx, ny, px, py)
		!! statically extrapolate in x and periodic in y
		!! assign value to ghost cells
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l
		
		!$omp parallel
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					! assign left ghost cells
					u_con(:,k,l,0,j) = u_con(:,1,l,1,j)
					
					! assign right ghost cells
					u_con(:,k,l,nx+1,j) = u_con(:,px+1,l,nx,j)
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			u_con(:,:,:,i,0) = u_con(:,:,:,i,ny) ! assign bottom ghost cells
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			u_con(:,:,:,i,ny+1) = u_con(:,:,:,i,1) ! assign top ghost cells
		end do
		!$omp end do nowait
		!$omp end parallel
	end subroutine BC_BrioWu

	subroutine BC_YeeSjogreen(u_con, xx, yy, nx, ny, px, py)
		!! statically extrapolate
		!! assign value to ghost cells
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l
		
		!$omp parallel
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					! assign left ghost cells
					u_con(:,k,l,0,j) = u_con(:,1,l,1,j)
					
					! assign right ghost cells
					u_con(:,k,l,nx+1,j) = u_con(:,px+1,l,nx,j)
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			do l = 1, py+1
				do k = 1, px+1
					! assign bottom ghost cells
					u_con(:,k,l,i,0) = u_con(:,k,1,i,1)
					
					! assign top ghost cells
					u_con(:,k,l,i,ny+1) = u_con(:,k,py+1,i,ny)
				end do
			end do
		end do
		!$omp end do nowait
		!$omp end parallel
	end subroutine BC_YeeSjogreen

	subroutine BC_cloud_shock(u_con, xx, yy, nx, ny, px, py)
		!! assign value to ghost cells
		!! new version: copy the whole boundary element to outside???
		use MHD_PDE_m, only: MHD_pri2con
		use utility, only: pi => PI_math
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l, s
		real, dimension(8) :: u_con_inflow, u_pri_inflow
		
		real, parameter :: rr = 3.868605300827511504664324644266228 ! compression ratio, -(5/6) (123 + 103 \[Pi]) + 5/6 Sqrt[16281 + 26298 \[Pi] + 10609 \[Pi]^2]
		real, parameter :: cf = 15.17658208897592333562014288162906 ! fast magnetosonic speed on right side, gamma must be 5.0/3.0
		real, parameter :: vl = 11.25357084620851400613837284504711 ! left speed
		real, parameter :: byl = 2.182626813584516145993203496356033 ! Byl
		real, parameter :: pl = 167.3451918203954512529698897850511 ! left pressure
		
		! u_pri_inflow = [1.0, -vl, 0.0, 0.0, 0.0, 1.0/sqrt(pi), 1.0/sqrt(pi), 1.0]
		u_pri_inflow = [1.0, -11.2536, 0.0, 0.0, 0.0, 0.56418958, 0.56418958, 1.0]
		
		u_con_inflow = MHD_pri2con(u_pri_inflow)
		
		!$omp parallel
		!$omp do schedule(static)
		do j = 1, ny
			! assign left ghost cells: copy insider element directly
			u_con(:,:,:,0,j) = u_con(:,:,:,1,j)
			
			! assign right ghost cells (inflow)
			do l = 1, py+1
				do k = 1, px+1
					u_con(:,k,l,nx+1,j) = u_con_inflow(:)
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			! assign bottom ghost cells: copy insider element directly
			u_con(:,:,:,i,0) = u_con(:,:,:,i,1)
			
			! assign top ghost cells: copy insider element directly
			u_con(:,:,:,i,ny+1) = u_con(:,:,:,i,ny)
		end do
		!$omp end do nowait
		!$omp end parallel
	end subroutine BC_cloud_shock
	
	subroutine BC_jet1_cell_ext(u_con, xx, yy, nx, ny, px, py)
		!! assign cell values directly to ghost elements
		!! according to private communication with Prof. Yong Liu
		!! assign value to ghost cells
		use MHD_PDE_m, only: gamma, MHD_pri2con
		use gauss_m
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		integer :: i, j, k, l, s
		real, parameter :: B2 = sqrt(200.0)
		real, parameter :: v2 = 800.0
		real, dimension(8) :: u_con_nozzle, u_pri_nozzle, u_con_ambient, u_pri_ambient
		
		real, dimension(px+1) :: qx_x, qw_x
		real, dimension(py+1) :: qx_y, qw_y
		
		call gauss_lobatto(px+1, qx_x, qw_x)
		call gauss_lobatto(py+1, qx_y, qw_y)
		
		u_pri_nozzle = [gamma, 0.0, v2, 0.0, 0.0, B2, 0.0, 1.0]
		u_pri_ambient = [0.1*gamma, 0.0, 0.0, 0.0, 0.0, B2, 0.0, 1.0]
		
		u_con_nozzle = MHD_pri2con(u_pri_nozzle)
		u_con_ambient = MHD_pri2con(u_pri_ambient)
		
		!$omp parallel
		
		! assign left ghost cells (by mirror symmetry) (only compute right part)
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					! mirror symmetry boundary condition in x-direction
					! (not compatible with the magnetic vector potential formulation)
					! MUST BE EXAMINED CAREFULLY TO ENSURE CONSISTENCY WITH THE PROBLEM SETTING
					
					! mirror of scalar field (symmetric)
					u_con(1,k,l,0,j) = u_con(1,px+2-k,l,1,j)
					
					! mirror of normal velocity (skew-symmetric)
					u_con(2,k,l,0,j) = -u_con(2,px+2-k,l,1,j)
					
					! mirror of tangential velocity (symmetric)
					u_con(3,k,l,0,j) = u_con(3,px+2-k,l,1,j)
					u_con(4,k,l,0,j) = u_con(4,px+2-k,l,1,j)
					
					! mirror of normal magnetic field (skew-symmetric)
					u_con(5,k,l,0,j) = -u_con(5,px+2-k,l,1,j)
					
					! mirror of tangential magnetic field (symmetric)
					u_con(6,k,l,0,j) = u_con(6,px+2-k,l,1,j)
					u_con(7,k,l,0,j) = u_con(7,px+2-k,l,1,j)
					
					! mirror of scalar field (symmetric)
					u_con(8,k,l,0,j) = u_con(8,px+2-k,l,1,j)
				end do
			end do
		end do
		!$omp end do nowait
		
		! assign right ghost cells
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					! copy insider element directly according to Prof. Yong Liu
					u_con(:,k,l,nx+1,j) = u_con(:,k,l,nx,j)
				end do
			end do
		end do
		!$omp end do nowait
		
		! assign bottom ghost cells
		!$omp do schedule(static)
		do i = 1, nx
			if (i <= nx/10) then ! WARNING: require nx to be multiple of 10!
				do l = 1, py+1
					do k = 1, px+1
						! bottom inflow nozzle region (should be an entire element to ensure LDF is maintained)
						u_con(:,k,l,i,0) = u_con_nozzle(:)
					end do
				end do
			else
				block
				! real, dimension(8) :: u_avg_
				! u_avg_ = 0.0
				! do l = 1, py+1
				! 	do k = 1, px+1
				! 		u_avg_ = u_avg_ + (0.25*qw_x(k)*qw_y(l)) * u_con(:,k,l,i,1)
				! 	end do
				! end do
				
				do l = 1, py+1
					do k = 1, px+1
						! copy insider element directly according to Prof. Yong Liu
						u_con(:,k,l,i,0) = u_con(:,k,l,i,1)
						
						! use cell average (very dissipative, but no help to the LDF phenomenon)
						! u_con(:,k,l,i,0) = u_avg_
					end do
				end do
				end block
			end if
		end do
		!$omp end do nowait
		
		! assign top ghost cells
		!$omp do schedule(static)
		do i = 1, nx
			do l = 1, py+1
				do k = 1, px+1
					! copy insider element directly according to Prof. Yong Liu
					u_con(:,k,l,i,ny+1) = u_con(:,k,l,i,ny)
				end do
			end do
		end do
		!$omp end do nowait
		
		!$omp end parallel
	end subroutine BC_jet1_cell_ext
	
	subroutine BC_rotor(u_con, xx, yy, nx, ny, px, py)
		!! assign value to ghost cells
		use MHD_PDE_m, only: MHD_pri2con
		use utility, only: pi=>PI_math
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
		real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		real, dimension(8) :: u_pri, u_con_
		integer :: i, j, k, l
		u_pri(1) = 1.0
		u_pri(2) = 0.0
		u_pri(3) = 0.0
		u_pri(4) = 0.0
		u_pri(5) = 5.0/sqrt(4*pi)
		u_pri(6) = 0.0
		u_pri(7) = 0.0
		u_pri(8) = 1.0
		u_con_ = MHD_pri2con(u_pri)
		!$omp parallel
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					u_con(:,k,l,0,j) = u_con_ ! assign left ghost cells
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do j = 1, ny
			do l = 1, py+1
				do k = 1, px+1
					u_con(:,k,l,nx+1,j) = u_con_ ! assign right ghost cells
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			do l = 1, py+1
				do k = 1, px+1
					u_con(:,k,l,i,0) = u_con_ ! assign bottom ghost cells
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do schedule(static)
		do i = 1, nx
			do l = 1, py+1
				do k = 1, px+1
					u_con(:,k,l,i,ny+1) = u_con_ ! assign top ghost cells
				end do
			end do
		end do
		!$omp end do nowait
		!$omp end parallel
	end subroutine BC_rotor
	
	!-------------------------------------------------------------------------!
	! compute error
	
	subroutine compute_error_smooth_alfven(u_con, hx, hy, xx, yy, nx, ny, px, py, err_l2, err_max)
		! error in B_perp
		use utility, only: pi=>PI_math
		use gauss_m, only: gauss_lobatto
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: hx, hy
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny)
		real, dimension(8,px+1,py+1,0:nx+1,0:ny+1), intent(in) :: u_con ! must be contiguous
		real, intent(out) :: err_l2, err_max
		
		real, dimension(px+1) :: qx_x, qw_x
		real, dimension(py+1) :: qx_y, qw_y
		real, parameter :: alpha = pi/6.0
		integer :: i, j, k, l
		
		call gauss_lobatto(px+1, qx_x, qw_x)
		call gauss_lobatto(py+1, qx_y, qw_y)
		err_l2 = 0.0
		err_max = 0.0
		! disable OMP to ensure reproducibility
		do j = 1, ny
			do i = 1, nx
				do l = 1, py+1
					do k = 1, px+1
						block
							real :: x_para, B_perp, Bh_perp
							x_para = xx(k,i) * cos(alpha) + yy(l,j) * sin(alpha)
							B_perp = 0.1 * sin(2*pi*x_para)
							Bh_perp = u_con(5,k,l,i,j) * (-sin(alpha)) + u_con(6,k,l,i,j) * cos(alpha)
							err_l2 = err_l2 + (0.25 * hx * hy * qw_x(k) * qw_y(l)) * (Bh_perp - B_perp)**2
							err_max = max(err_max, abs(Bh_perp - B_perp))
						end block
					end do
				end do
			end do
		end do
		err_l2 = sqrt(err_l2)
	end subroutine compute_error_smooth_alfven
	
	subroutine compute_error_smooth1D(u_con, hx, hy, xx, yy, nx, ny, px, py, T, err_l2, err_max)
		! error in B_perp
		use gauss_m, only: gauss_lobatto
		use MHD_PDE_m, only: MHD_pri2con
		implicit none
		integer, intent(in) :: nx, ny, px, py
		real, intent(in) :: hx, hy, T
		real, intent(in) :: xx(px+1, nx), yy(py+1, ny)
		real, dimension(8,px+1,py+1,0:nx+1,0:ny+1), intent(in) :: u_con ! must be contiguous
		real, intent(out) :: err_l2, err_max
		
		real, dimension(px+1) :: qx_x, qw_x
		real, dimension(py+1) :: qx_y, qw_y
		integer :: i, j, k, l
		
		call gauss_lobatto(px+1, qx_x, qw_x)
		call gauss_lobatto(py+1, qx_y, qw_y)
		err_l2 = 0.0
		err_max = 0.0
		! disable OMP to ensure reproducibility
		do j = 1, ny
			do i = 1, nx
				do l = 1, py+1
					do k = 1, px+1
						block
							real :: exact(8), u_pri(8)
							u_pri(1) = 1.0 + 0.2*sin(xx(k,i)-T)
							u_pri(2) = 1.0
							u_pri(3) = 0.0
							u_pri(4) = 0.0
							u_pri(5) = 0.5
							u_pri(6) = 1.0
							u_pri(7) = 1.5
							u_pri(8) = 2.0
							exact = MHD_pri2con(u_pri)
							err_l2 = err_l2 + (0.25*hx*hy*qw_x(k)*qw_y(l)) * sum((u_con(:,k,l,i,j) - exact)**2)
							err_max = max(err_max, maxval(abs(u_con(:,k,l,i,j) - exact)))
						end block
					end do
				end do
			end do
		end do
		err_l2 = sqrt(err_l2)
	end subroutine compute_error_smooth1D
	
end module problems