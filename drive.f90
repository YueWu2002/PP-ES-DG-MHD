module drive_m
	! total control
	implicit none
	
	abstract interface
		subroutine Init_type(u_con, xx, yy, nx, ny, px, py)
			!! assign value to interior cells
			implicit none
			integer, intent(in) :: nx, ny, px, py
			real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
			real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		end subroutine Init_type
		
		subroutine BC_type(u_con, xx, yy, nx, ny, px, py)
			!! assign value to ghost cells (not only the ghost nodes)
			implicit none
			integer, intent(in) :: nx, ny, px, py
			real, intent(in) :: xx(px+1, nx), yy(py+1, ny) ! must be contiguous
			real, intent(inout) :: u_con(8, px+1, py+1, 0:nx+1, 0:ny+1) ! must be contiguous & with one layer of ghost cells
		end subroutine BC_type
	end interface
	
	! pointers to functions of the declared types
	procedure(Init_type),	pointer, public :: Init_handle	=> null()
	procedure(BC_type),		pointer, public :: BC_handle	=> null()
	
	public	:: Init_type, BC_type ! expose the interface to any place that uses this module
	public	:: solve_execute
	
contains

	subroutine solve_execute(final_T, CFL, use_oedg, use_LDF, ent_offset)
		use utility
		use DG_solver_m
		implicit none
		real, intent(in) :: final_T, CFL
		logical, intent(in) :: use_oedg, use_LDF, ent_offset
		
		! time marching variables
		logical :: flag
		real :: t, dt, dt_lim, CFL_RK
		real, dimension(rknstg) :: t_RK ! pseudo-time in the Runge-Kutta procedure
		integer :: i, j, k, l, s
		
		! diagnostic data
		real :: tic_cpu_time, toc_cpu_time, tot_ent
		integer :: time_step_num, retry_num, num_negcell
		logical :: neg_avg
		
		! data logging
		integer :: log_count
		real, dimension(:,:), allocatable :: my_log ! history of (time, entropy)
		
		if (.not. associated(Init_handle)) error stop "drive: solve_execute: procedure pointer 'Init_handle' not initialized!"
		if (.not. associated(BC_handle)) error stop "drive: solve_execute: procedure pointer 'BC_handle' not initialized!"
		
		if (rknstg <= 3) then
			CFL_RK = 1.0
		else if (rknstg == 5) then
			CFL_RK = 1.508 ! SSP-RK(5,4)
		else
			error stop 'unsupported yet'
		end if
		
		! assign initial value for both the conservative and primitive variables
		call Init_handle(u_con, xx, yy, nx, ny, px, py)
		t = 0.0
		call BC_handle(u_con, xx, yy, nx, ny, px, py)
		call post_process(u_con, u_avg, neg_avg, num_negcell, apply_oedg=.false., apply_LDF=use_LDF, oedg_dt=0.0, ent_offset=.false., comput_ent=.true., tot_ent=tot_ent)
		if (neg_avg) error stop 'Error: Initial value has negative cell average!'
		call BC_handle(u_con, xx, yy, nx, ny, px, py) ! because post_process is not completely cell-wise
		! compute the time derivative for the current stage u_con
		call DG_semidiscretize(ut_con(:,:,:,:,:,1), dt_lim, u_con)
		write(stdout, "(A, ES14.7, A, ES14.7, A, ES14.7, A, ES14.7)") 't=', t, ', ent=', tot_ent, ', dt_lim=', dt_lim, ', time per it:', 0.0
		call log_append(my_log, log_count, [0.0, tot_ent])
		
		! Remark: We evolve the stage together with its time derivative and dt_lim. (In this way, the residue is also computed.)
		
		flag = .true.
		time_step_num = 0
		retry_num = 0
		call cpu_time(tic_cpu_time)
		solving_loop: do while(flag)
			time_step_num = time_step_num + 1
			
			try_RK_block: block
				integer :: try_num
				integer, parameter :: max_try_num = 10 ! must be bigger than 1
				logical :: try_again
				
				try_num = 0
				try_again = .true.
				try_RK_loop: do while(try_again)
					try_again = .false.
					try_num = try_num + 1
					
					if (try_num > max_try_num .or. dt_lim < 1e-14) then
						write(stderr, *)"Computation failed: time step converges to 0! dt_lim=", dt_lim
						error stop 1
					else if (try_num > 1) then
						! This is a re-try. We need to re-compute the conservative variables at ghost points and the primitive variables and re-compute the time derivative for u_con.
						retry_num = retry_num + 1 ! increase the global retry number
						! re-compute the time derivative for the current stage u_con
						call DG_semidiscretize(ut_con(:,:,:,:,:,1), dt_lim, u_con)
					end if
					
					if (try_num <= 1) then
						! only use this CFL condition for the first try (linear stability)
						dt = CFL_RK * CFL * dt_lim
					else
						! retry with a stricter time step (shrink with golden ratio)
						dt = ((sqrt(1.25)-0.5)**(try_num-1)) * dt
					end if
					
					if (t + dt >= final_T) then
						dt = final_T - t
						flag = .false. ! get ready for termination
					else
						flag = .true. ! reserved for cases where dt is reduced by a re-try loop
					end if
					
					RK_block: block
						integer :: rks_id ! ID of the current to-compute Runge-Kutta stage
						logical :: rk_final
						real :: oedg_beta ! coefficient used for dt in the OEDG limiter
						
						! For the coefficient beta, we use the following rule to determine it, which agrees with "An entropy stable essentially oscillation-free discontinuous Galerkin method for solving ideal magnetohydrodynamic equations":
						! For a stage, if it can be rewritten as convex combination of past stages plus some time derivatives, then beta will be the dt coefficient of that derivative. 
						
						RK_loop : do rks_id = 1, rknstg
							select case (rknstg)
							case (1) ! SSP-RK(1,1) / Euler forward
								select case (rks_id)
								case (1)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,1) = u_con(:,:,:,i,j) + dt*ut_con(:,:,:,i,j,1)
										end do
									end do
									!$omp end parallel do
									t_RK(1) = t + dt*1.0
									oedg_beta = 1.0
									rk_final = .true.
								case default
									error stop 'Invalid stage for SSPRK(1,1)!'
								end select
							case (2) ! SSP-RK(2,2) / Heun's method
								select case (rks_id)
								case (1)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,1) = u_con(:,:,:,i,j) + dt*ut_con(:,:,:,i,j,1)
										end do
									end do
									!$omp end parallel do
									t_RK(1) = t + dt*1.0
									oedg_beta = 1.0
									rk_final = .false.
								case (2)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,2) = 0.5*u_con(:,:,:,i,j) + 0.5*(uRK_con(:,:,:,i,j,1) + dt*ut_con(:,:,:,i,j,2))
										end do
									end do
									!$omp end parallel do
									t_RK(2) = 0.5*t + 0.5*(t_RK(1) + dt*1.0)
									oedg_beta = 0.5
									rk_final = .true.
								case default
									error stop 'Invalid stage for SSPRK(2,2)!'
								end select
							case (3) ! SSP-RK(3,3)
								select case (rks_id)
								case (1)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,1) = u_con(:,:,:,i,j) + dt*ut_con(:,:,:,i,j,1)
										end do
									end do
									!$omp end parallel do
									t_RK(1) = t + dt*1.0
									oedg_beta = 1.0
									rk_final = .false.
								case (2)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,2) = 0.75*u_con(:,:,:,i,j) + 0.25*(uRK_con(:,:,:,i,j,1) + dt*ut_con(:,:,:,i,j,2))
										end do
									end do
									!$omp end parallel do
									t_RK(2) = 0.75*t + 0.25*(t_RK(1) + dt*1.0)
									oedg_beta = 0.25
									rk_final = .false.
								case (3)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,3) = (1.0/3.0)*u_con(:,:,:,i,j) + (2.0/3.0)*(uRK_con(:,:,:,i,j,2) + dt*ut_con(:,:,:,i,j,3))
										end do
									end do
									!$omp end parallel do
									t_RK(3) = (1.0/3.0)*t + (2.0/3.0)*(t_RK(2) + dt*1.0)
									oedg_beta = 2.0/3.0
									rk_final = .true.
								case default
									error stop 'Invalid stage for SSPRK(3,3)!'
								end select
							case (5) ! SSP-RK(5,4) (Ref: s10915-008-9239-z)
								select case (rks_id)
								case (1)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,1) = u_con(:,:,:,i,j) + (0.391752226571890*dt)*ut_con(:,:,:,i,j,1)
										end do
									end do
									!$omp end parallel do
									t_RK(1) = t + (0.391752226571890*dt)*1.0
									oedg_beta = 0.391752226571890
									rk_final = .false.
								case (2)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,2) = 0.444370493651235*u_con(:,:,:,i,j) + 0.555629506348765*uRK_con(:,:,:,i,j,1) + (0.368410593050371*dt)*ut_con(:,:,:,i,j,2)
										end do
									end do
									!$omp end parallel do
									t_RK(2) = 0.444370493651235*t + 0.555629506348765*t_RK(1) + (0.368410593050371*dt)*1.0
									oedg_beta = 0.368410593050371
									rk_final = .false.
								case (3)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,3) = 0.620101851488403*u_con(:,:,:,i,j) + 0.379898148511597*uRK_con(:,:,:,i,j,2) + (0.251891774271694*dt)*ut_con(:,:,:,i,j,3)
										end do
									end do
									!$omp end parallel do
									t_RK(3) = 0.620101851488403*t + 0.379898148511597*t_RK(2) + (0.251891774271694*dt)*1.0
									oedg_beta = 0.251891774271694
									rk_final = .false.
								case (4)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,4) = 0.178079954393132*u_con(:,:,:,i,j) + 0.821920045606868*uRK_con(:,:,:,i,j,3) + (0.544974750228521*dt)*ut_con(:,:,:,i,j,4)
										end do
									end do
									!$omp end parallel do
									t_RK(4) = 0.178079954393132*t + 0.821920045606868*t_RK(3) + (0.544974750228521*dt)*1.0
									oedg_beta = 0.544974750228521
									rk_final = .false.
								case (5)
									!$omp parallel do collapse(2) schedule(static)
									do j = 1, ny
										do i = 1, nx
											uRK_con(:,:,:,i,j,5) = 0.517231671970585*uRK_con(:,:,:,i,j,2) + 0.096059710526147*uRK_con(:,:,:,i,j,3) + (0.063692468666290*dt)*ut_con(:,:,:,i,j,4) + 0.386708617503269*uRK_con(:,:,:,i,j,4) + (0.226007483236906*dt)*ut_con(:,:,:,i,j,5)
										end do
									end do
									!$omp end parallel do
									t_RK(5) = 0.517231671970585*t_RK(2) + 0.096059710526147*t_RK(3) + (0.063692468666290*dt)*1.0 + 0.386708617503269*t_RK(4) + (0.226007483236906*dt)*1.0
									oedg_beta = 0.063692468666290 + 0.226007483236906
									rk_final = .true.
								case default
									error stop 'Invalid stage for SSPRK(5,4)!'
								end select
							case default
								error stop 'Unsupported RK num stages!'
							end select
							
							call BC_handle(uRK_con(:,:,:,:,:,rks_id), xx, yy, nx, ny, px, py)
							
							! apply OEDG at time step (more efficient)
							call post_process(uRK_con(:,:,:,:,:,rks_id), u_avg, neg_avg, num_negcell, apply_oedg=(use_oedg .and. rk_final), apply_LDF=use_LDF, oedg_dt=dt, ent_offset=ent_offset, comput_ent=rk_final, tot_ent=tot_ent)
							
							! apply OEDG at each stage (very costly) (not much difference)
							! call post_process(uRK_con(:,:,:,:,:,rks_id), u_avg, neg_avg, num_negcell, apply_oedg=use_oedg, apply_LDF=use_LDF, oedg_dt=(oedg_beta*dt), ent_offset=ent_offset, comput_ent=rk_final, tot_ent=tot_ent)
							
							if (neg_avg) then
								! There is negativity in the cell average!
								write(stderr, '(A, g0, A, g0)') "RKstg ", rks_id, ": Inadmissible average! Restart RK with reduced dt. cur_try_num=", try_num
								try_again = .true.
								exit RK_loop
							end if
							call BC_handle(uRK_con(:,:,:,:,:,rks_id), xx, yy, nx, ny, px, py) ! because post_process is not completely cell-wise
							
							if (num_negcell > 0) then
								write(stderr, '(A, g0, A, g0, A, g0)') "RKstg ", rks_id, ": PP limited: ", num_negcell, " cells, cur_try_num=", try_num
							end if
							! If the limiter is correct, there will not be any negative nodal values.
							! compute the time derivative for the current stage uRK_con
							call DG_semidiscretize(ut_con(:,:,:,:,:,mod(rks_id,rknstg)+1), dt_lim, uRK_con(:,:,:,:,:,rks_id))
						end do RK_loop
					end block RK_block
				end do try_RK_loop
			end block try_RK_block ! RK success here
			
			!$omp parallel do collapse(2) schedule(static)
			do j = 1, ny
				do i = 1, nx
					! update u with the latest RK stage (and also accept the dt_lim, dt_lim_PP, tot_ent and ut_con implicitly)
					u_con(:,:,:,i,j) = uRK_con(:,:,:,i,j,rknstg)
				end do
			end do
			!$omp end parallel do
			t = t_RK(rknstg)
			
			call log_append(my_log, log_count, [t, tot_ent])
			
			if (mod(time_step_num, 256)==0 .or. (flag.eqv..false.)) then
				call cpu_time(toc_cpu_time)
				write(stdout, "(A, ES14.7, A, ES14.7, A, ES14.7, A, ES14.7)") 't=', t, ', ent=', tot_ent, ', dt_lim=', dt_lim, ', time per it:', (toc_cpu_time - tic_cpu_time) / time_step_num
			end if
		end do solving_loop
		! call cpu_time(toc_cpu_time) ! because the CPU time is already recorded above
		
		! These are not supported by GCC version 11, but supported by newer versions. We can enable them when we switch to a newer compiler.
		!Init_handle => null()
		!BC_handle => null()
		
		! final outputs
		block
			use io
			real :: cpu_time_duration
			character(len=128) :: filename, zonename
			cpu_time_duration = toc_cpu_time - tic_cpu_time
			write(stdout, '(A)') "drive: solve_execute: Solving completed."
			write(stdout, '(A, g0)') "# time steps = ", time_step_num
			write(stdout, '(A, g0)') "# re-trys = ", retry_num
			write(stdout, '(A, ES11.4, A)') "total CPU time = ", cpu_time_duration, " s"
			write(stdout, '(A, ES11.4, A)') "CPU time per RK = ", cpu_time_duration / time_step_num, " s"
			
			write(filename, '(A)') 'entropy.plt'
			write(zonename, '(I0, 1X, I0, 1X, I0, 1X, I0)') nx, ny, px, py
			call save_tecbin_entlog(my_log, log_count, filename, zonename)
		end block
	end subroutine solve_execute

	pure subroutine log_append(arr, used_size, new_vec)
		!! append data at the end of arr(:,1,...,used_size), where the array is actually arr(:,1,...,max_size)
		implicit none
		real, dimension(:,:), allocatable, intent(inout) :: arr
		integer, intent(inout) :: used_size
		real, dimension(:), intent(in) :: new_vec
		real, dimension(:,:), allocatable :: temp
		integer :: max_size, vec_size
		
		vec_size = size(new_vec)
		
		if (.not.allocated(arr)) then
			allocate(arr(vec_size, 8))
			used_size = 1
			arr(:,1) = new_vec(:)
			return
		end if
		
		used_size = used_size + 1
		max_size = size(arr, 2)
		if (used_size > max_size) then
			allocate(temp(vec_size, 2*max_size))
			temp(:, 1:used_size-1) = arr(:, 1:used_size-1)
			call move_alloc(temp, arr) ! move temp to arr and deallocate the original arr
		end if
		arr(:, used_size) = new_vec(:)
	end subroutine log_append
	
end module drive_m