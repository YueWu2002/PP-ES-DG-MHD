program main
	use gauss_m
	use utility
	use drive_m
	use DG_solver_m
	use problems
	use io
	use MHD_PDE_m
	implicit none
	character(len=100) :: problem_name
	! select the problem model
	integer :: Problem_ID ! 0: smooth Alfven wave
	real :: left, right, bottom, top
	real :: final_T, CFL, err_l2, err_max
	!-------------------------------------------------------------------------!
	program_para_output: block
		use utility
		integer :: file_unit
		open(newunit=file_unit, file="args.log", status="replace", action="write")
		call compiler_para_output(file_unit)
		call omp_para_output(file_unit)
		close(file_unit)
	end block program_para_output
	!-------------------------------------------------------------------------!
	write(stdout, *) gauss_lobatto_test_legendre()
	block
		integer :: p
		do p = 0, 6
			block
			real, dimension(0:p) :: x, w
			real, dimension(0:p, 0:p) :: M, Minv, D, Pmn, Pnm
			call gauss_lobatto(p+1, x, w)
			call gauss_lobatto_matrix(p, M, Minv, D, Pmn, Pnm)
			call test_lobatto_mat(p, x, w, M, Minv, D, Pmn, Pnm)
			end block
		end do
		print *, "Gauss-Lobatto verification finished."
	end block
	block
		integer :: p
		do p = 0, 6
			block
			real, dimension(0:p) :: m, minv, tl, tr
			real, dimension(0:p, 0:p) :: D
			call gauss_legendre_matrix(p, m, minv, D, tl, tr)
			call test_legendre_mat(p, m, D, tl, tr)
			end block
		end do
		print *, "Gauss-Lobatto verification finished."
	end block
	!-------------------------------------------------------------------------!
	
	CFL = 0.6
	
	Problem_ID = 0
	select case (Problem_ID)
	case (0)
		problem_name = "1D_Brio_Wu_shock_tube"
		left = -1.0; right = 1.0; bottom = -0.5; top = 0.5
		final_T = 0.2
		Init_handle => init_by_fun_pri
		init_pri_fun => Brio_Wu_init_fun
		BC_handle => BC_BrioWu
		call init_DG_solver(480, 1, 3, 0, 5, left, right, bottom, top) ! 1D problem
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (1)
		problem_name = "smooth_1D_density_wave"
		left = 0.0; right = 2*PI_math; bottom = -0.5; top = 0.5
		final_T = 1.3
		Init_handle => init_by_fun_pri
		init_pri_fun => smooth_1Dwave_init_fun
		BC_handle => BC_periodic
		call init_DG_solver(480, 1, 3, 0, 5, left, right, bottom, top) ! 1D problem
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .false., .false.)
		call compute_error_smooth1D(u_con, hx, hy, xx, yy, nx, ny, px, py, final_T, err_l2, err_max)
		print *, 'l2_err = ', err_l2
		print *, 'max_err = ', err_max
	case (2)
		problem_name = "smooth_Alfven_wave"
		left = 0.0; right = 1.0/cos(PI_math/6.0); bottom = 0.0; top = 1.0/sin(PI_math/6.0)
		final_T = 5.0
		Init_handle => init_by_fun_pri
		init_pri_fun => smooth_alfven_init_fun
		BC_handle => BC_periodic
		call init_DG_solver(128, 128, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .false., .false.)
		call compute_error_smooth_alfven(u_con, hx, hy, xx, yy, nx, ny, px, py, err_l2, err_max)
		print *, 'l2_err = ', err_l2
		print *, 'max_err = ', err_max
	case (3)
		problem_name = "2D Riemann Yee Sjogreen"
		left = -1.0; right = 1.0; bottom = -1.0; top = 1.0
		final_T = 0.2
		Init_handle => init_by_fun_con
		init_con_fun => YeeSjogreenInitCon
		BC_handle => BC_YeeSjogreen
		call init_DG_solver(512, 512, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (4)
		problem_name = "Orszag-Tang vortex"
		left = 0.0; right = 1.0; bottom = 0.0; top = 1.0
		final_T = 0.5
		Init_handle => init_by_fun_pri
		init_pri_fun => OrszagTang_init_pri
		BC_handle => BC_periodic
		call init_DG_solver(128, 128, 3, 3, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (5)
		problem_name = "rotor"
		left = 0.0; right = 1.0; bottom = 0.0; top = 1.0
		final_T = 0.15
		call set_gamma(1.4)
		Init_handle => init_by_fun_pri
		init_pri_fun => Rotor_init_pri
		BC_handle => BC_rotor
		call init_DG_solver(256, 256, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (6)
		problem_name = "blast_wave"
		left = 0.0; right = 1.0; bottom = 0.0; top = 1.0
		final_T = 0.01
		call set_gamma(1.4)
		Init_handle => init_by_fun_pri
		init_pri_fun => Blast_init_pri
		BC_handle => BC_periodic ! according to Balsara's GDF-DG-MHD paper
		call init_DG_solver(256, 256, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (7)
		problem_name = "cloud_shock" ! original version (bubble hit the flow from the right)
		! WARNING: cannot use mirror symmetry to half the domain, because By is not skew symmetric.
		left = 0.0; right = 1.0; bottom = 0.0; top = 1.0
		final_T = 0.06
		Init_handle => init_by_fun_pri
		init_pri_fun => cloud_shock_init_pri
		BC_handle => BC_cloud_shock
		! Lowering CFL will not affect oscillation.
		! Must use Prof. Kailiang Wu's original version of OEDG to suppress oscillation. 
		! Simplified version that only consider jump in value and derivative does not suffice. 
		! Jump filter-type convex blending is too dissipative.
		call init_DG_solver(512, 512, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (8)
		! We compute only the rigth half. (using mirror symmetry)
		problem_name = "Astrophysical-jet"
		left = 0.0; right = 0.5; bottom = 0.0; top = 1.5
		final_T = 2.0e-3
		call set_gamma(1.4)
		Init_handle => init_by_fun_pri
		init_pri_fun => jet1_init_pri
		BC_handle => BC_jet1_cell_ext
		! We have verified that the mirror symmetry BC is implemented correctly.
		call init_DG_solver(200, 600, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case (9)
		problem_name = "Kelvin-Helmholtz instability"
		left = 0.0; right = 1.0; bottom = -1.0; top = 1.0
		final_T = 0.5
		call set_gamma(1.4)
		Init_handle => init_by_fun_pri
		init_pri_fun => KHI_init_pri
		BC_handle => BC_periodic
		call init_DG_solver(128, 256, 2, 2, 3, left, right, bottom, top)
		print *, 'Problem ID=', Problem_ID, ', nx=', nx, ', ny=', ny, ', px=', px, ', py=', py
		call solve_execute(final_T, CFL, .true., .true., .false.)
		call save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, 0, final_T)
	case default
		error stop 'invalid Problem_ID!'
	end select
	call finalize_DG_solver
	print '(A)', 'Program MAIN terminated successfully.'
end program main