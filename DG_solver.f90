module DG_solver_m
	! spatial discretization (no BC control) (no actual storage of solution)
	use LDF_m
	implicit none
	
	! public data: for computation (use AoS storage as FLEXI does)
	! u_avg has size (8, 0:nx+1, 0:ny+1) (one layer of ghost cells)
	! u_con has size (8, px+1, py+1, 0:nx+1, 0:ny+1) (one layer of ghost cells)
	! ut_con has size (8, px+1, py+1, nx, ny, rknstg)
	! uRK_con has size (8, px+1, py+1, 0:nx+1, 0:ny+1, rknstg) (one layer of ghost cells)
	real, dimension(:,:,:), allocatable, public :: u_avg				! cell average data
	real, dimension(:,:,:,:,:), allocatable, public :: u_con			! conservatie variables / numerical solution
	real, dimension(:,:,:,:,:,:), allocatable, public :: ut_con			! time derivative
	real, dimension(:,:,:,:,:,:), allocatable, public :: uRK_con		! Runge Kutta stages
	
	! private data: internal working storage
	real, parameter, public :: PP_eps = 1.0e-11								! for PP limiter (too small: bad for FP stability)
	
	real, protected :: domain_left, domain_right, domain_top, domain_bottom	! boundary of domain
	integer, protected :: nx, ny, px, py									! number of element and polynomial degree
	integer, protected :: rknstg											! number of RK stages
	real, protected :: hx, hy
	real, dimension(:), allocatable, protected :: grid_x, grid_y			! boundary of cells
	real, dimension(:,:), allocatable, protected :: xx, yy					! collocation nodes
	
	! interface data
	real, dimension(:,:,:),		allocatable, private :: sxl, sxr, syd, syu		! signal speed estimates at interface nodes
	real, dimension(:,:,:,:),	allocatable, private :: nflux_x, nflux_y		! numerical flux at interface nodes
	real, dimension(:,:,:),		allocatable, private :: nB_x, nB_y				! numerical flux at interface nodes
	real, dimension(:,:,:),		allocatable, private :: OEDG_sigma_x, OEDG_sigma_y		! Yong Liu's OEDG indicators (sum) at interfaces (non affine-invariant)
	real, dimension(:,:,:),		allocatable, private :: OEDG_sigsum_x, OEDG_sigsum_y	! Kailiang Wu's OEDG indicators (summed) at interfaces (non affine-invariant) for jump of equal max degree
	
	! volume data
	logical, dimension(:,:),	allocatable, private :: cell_ent_nan			! whether the cell math entropy is NaN
	real, dimension(:,:),		allocatable, private :: cell_ent_DG				! cell integrated entropy
	real, dimension(:,:),		allocatable, private :: cell_ent_FV				! cell integrated entropy of cell average
	
	real, dimension(:), save,	allocatable, protected :: qx_x, qx_y, qw_x, qw_y		! Gauss Lobatto nodes and weights
	real, dimension(:,:), save,	allocatable, protected :: diff_x, diff_y				! Gauss Lobatto differentiation matrix
	real, dimension(:,:), save,	allocatable, protected :: M_x, Minv_x, M_y, Minv_y		! Gauss Lobatto mass matrix
	real, dimension(:,:), save,	allocatable, protected :: Pmn_x, Pmn_y, Pnm_x, Pnm_y	! Gauss Lobatto mapping matrix
	
	type(LDF_projector), save,	allocatable, private :: LDFproj_			! LDF projector
	
contains

	subroutine init_DG_solver(nx_in, ny_in, px_in, py_in, rknstg_in, domain_left_in, domain_right_in, domain_bottom_in, domain_top_in)
		use gauss_m
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer, intent(in) :: nx_in, ny_in, px_in, py_in, rknstg_in
		real, intent(in) :: domain_left_in, domain_right_in, domain_bottom_in, domain_top_in
		
		integer :: i, j
		real :: nan_ 
		nan_ = ieee_value(1.0, ieee_signaling_nan)
		
		nx = nx_in
		ny = ny_in
		px = px_in
		py = py_in
		rknstg = rknstg_in
		domain_left = domain_left_in
		domain_right = domain_right_in
		domain_bottom = domain_bottom_in
		domain_top = domain_top_in

		allocate(u_avg(8, 0:nx+1, 0:ny+1), source = nan_)
		allocate(u_con(8, px+1, py+1, 0:nx+1, 0:ny+1), source = nan_)
		allocate(ut_con(8, px+1, py+1, nx, ny, rknstg), source = nan_)
		allocate(uRK_con(8, px+1, py+1, 0:nx+1, 0:ny+1, rknstg), source = nan_)
		
		allocate(grid_x(nx+1), grid_y(ny+1), source = nan_)
		allocate(xx(px+1, nx), yy(py+1, ny), source = nan_)
		allocate(sxl(py+1,nx+1,ny), sxr(py+1,nx+1,ny), syd(px+1,nx,ny+1), syu(px+1,nx,ny+1), source = nan_)
		allocate(nflux_x(8,py+1,nx+1,ny), nflux_y(8,px+1,nx,ny+1), nB_x(py+1,nx+1,ny), nB_y(px+1,nx,ny+1), source = nan_)
		allocate(OEDG_sigma_x(8, nx+1, ny), OEDG_sigma_y(8, nx, ny+1), source = nan_)
		allocate(OEDG_sigsum_x(8, nx+1, ny), OEDG_sigsum_y(8, nx, ny+1), source = nan_)
		allocate(cell_ent_DG(0:nx+1, 0:ny+1), cell_ent_FV(0:nx+1, 0:ny+1), source = nan_)
		allocate(cell_ent_nan(0:nx+1, 0:ny+1), source = .true.)

		allocate(qx_x(px+1), qw_x(px+1), qx_y(py+1), qw_y(py+1), source = nan_)
		allocate(diff_x(px+1,px+1), diff_y(py+1,py+1), source = nan_)
		allocate(M_x(px+1,px+1), Minv_x(px+1,px+1), M_y(py+1,py+1), Minv_y(py+1,py+1), source = nan_)
		allocate(Pmn_x(px+1,px+1), Pnm_x(px+1,px+1), Pmn_y(py+1,py+1), Pnm_y(py+1,py+1), source = nan_)
		call gauss_lobatto(px+1, qx_x, qw_x)
		call gauss_lobatto(py+1, qx_y, qw_y)
		call gauss_lobatto_matrix(px, M_x, Minv_x, diff_x, Pmn_x, Pnm_x)
		call gauss_lobatto_matrix(py, M_y, Minv_y, diff_y, Pmn_y, Pnm_y)
		allocate(LDFproj_)
		call LDFproj_%init(px, py)
		
		hx = (domain_right - domain_left) / nx
		hy = (domain_top - domain_bottom) / ny
		do i = 1, nx+1
			grid_x(i) = domain_left + (i-1) * hx
		end do
		do j = 1, ny+1
			grid_y(j) = domain_bottom + (j-1) * hy
		end do
		do i = 1, nx
			xx(:,i) = (0.5*hx)*qx_x(:) + (domain_left + (i-0.5)*hx)
		end do
		do j = 1, ny
			yy(:,j) = (0.5*hy)*qx_y(:) + (domain_bottom + (j-0.5)*hy)
		end do
	end subroutine init_DG_solver
	
	subroutine finalize_DG_solver
		implicit none
		deallocate(u_avg, u_con, ut_con, uRK_con)
		deallocate(cell_ent_DG, cell_ent_FV, cell_ent_nan)
		
		deallocate(grid_x, grid_y)
		deallocate(xx, yy)
		deallocate(sxl, sxr, syd, syu)
		deallocate(nflux_x, nflux_y, nB_x, nB_y)
		deallocate(OEDG_sigma_x, OEDG_sigma_y)
		deallocate(OEDG_sigsum_x, OEDG_sigsum_y)
		
		deallocate(LDFproj_)
		deallocate(qx_x, qw_x, qx_y, qw_y, diff_x, diff_y, M_x, Minv_x, M_y, Minv_y, Pmn_x, Pnm_x, Pmn_y, Pnm_y)
	end subroutine finalize_DG_solver
	
#define _macro_prs_(rho, mx, my, mz, bx, by, bz, e_tot, gamma) (((gamma)-1.0)*((e_tot) - 0.5*(((mx)**2 + (my)**2 + (mz)**2)/(rho) + ((bx)**2 + (by)**2 + (bz)**2))))
#define _macro_ent_(rho, prs, gamma) (-(rho)*(log((prs)) - (gamma)*log((rho)))/((gamma)-1.0))

	subroutine post_process(u_con, u_avg, neg_avg, num_negcell, apply_oedg, apply_LDF, oedg_dt, ent_offset, comput_ent, tot_ent)
		!! find cell average, apply LDF proj and apply convex limiter (to only interior cells)
		!! Warning: need to handle BC again after calling this subroutine, since this procedure is not cell-wise.
		use MHD_PDE_m, only: gamma
		use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
		use utility
		implicit none
		real, dimension(8, px+1, py+1, 0:nx+1, 0:ny+1), intent(inout) :: u_con ! with one layer of ghost cells
		real, dimension(8, 0:nx+1, 0:ny+1), intent(inout) :: u_avg ! with one layer of ghost cells
		logical, intent(out) :: neg_avg ! found inadmissible cell average
		integer, intent(out) :: num_negcell ! number of negative cells / cells that PP limiter is activated
		logical, intent(in) :: apply_oedg ! whether should apply OEDG / OFDG limiter
		logical, intent(in) :: apply_LDF ! whether should apply LDF projection
		real, intent(in) :: oedg_dt ! dt used in OEDG convex limiter
		logical, intent(in) :: ent_offset ! whether should apply convex limiting to offset potential cell entropy increase due to LDF projection
		logical, intent(in) :: comput_ent ! whether should compute the total entropy or not
		real, intent(out) :: tot_ent ! total (math) entropy (by quadrature)
		
		real, dimension(8) :: u_avg_global, u_vari_global
		
		real, parameter :: OEDG_eps = 1.0e-8
		integer :: i,j,k,l,s
		
		neg_avg = .false.
		num_negcell = 0
		tot_ent = 0.0
		u_avg_global(:) = 0.0
		u_vari_global(:) = 0.0
		
		!$omp parallel
		
		! compute cell average (including ghost cells)
		! Reduction for Fortran arrays is supported on OpenMP 4.5 standard.
		!$omp do collapse(2) schedule(static) reduction(+: u_avg_global) reduction(.or.: neg_avg)
		do j = 0, ny+1
			do i = 0, nx+1
				if (neg_avg) cycle ! The rest part is meaningless. Restart the RK stage as soon as possible.
				if ((i>=1.and.i<=nx) .or. (j>=1.and.j<=ny)) then
				block
					real :: u_avg_(8)
					real, dimension(8, px+1, py+1) :: u_con_
					real :: den_avg, prs_avg, ent_avg
					real :: cell_ent_
					logical :: cell_ent_isnan_ ! whether the original cell total entropy is NaN (can be caused by negative density or pressure at some nodes)
					
					! read
					u_con_ = u_con(:,:,:,i,j)
					
					! compute the cell average
					u_avg_(:) = 0.0
					do l = 1, py+1
						do k = 1, px+1
							do s = 1, 8 ! completely unrolled by 8 by IFX
								u_avg_(s) = u_avg_(s) + (0.25*qw_x(k)*qw_y(l)) * u_con_(s,k,l)
							end do
						end do
					end do
					
					! variables evaluated at the cell average
					den_avg = u_avg_(1)
					prs_avg = _macro_prs_(u_avg_(1), u_avg_(2), u_avg_(3), u_avg_(4), u_avg_(5), u_avg_(6), u_avg_(7), u_avg_(8), gamma)
					if (den_avg < PP_eps .or. prs_avg < PP_eps) then
						neg_avg = .true.
						ent_avg = ieee_value(0.0, ieee_signaling_nan)
						cycle ! The rest part is meaningless. Restart the RK stage as soon as possible.
					end if
					
					! compute the cell total math entropy (cell_ent_old and cell_ent_nan_)
					cell_ent_isnan_ = .true.
					cell_ent_ = ieee_value(0.0, ieee_signaling_nan)
					if (ent_offset) then
						! guaranteed to be legal, because (den_avg >= PP_eps .and. prs_avg >= PP_eps) is just ensured above
						ent_avg = _macro_ent_(den_avg, prs_avg, gamma)
						block
							logical :: flag ! whether all nodes are admissible (positive density and pressure)
							real :: rho(px+1, py+1), prs(px+1, py+1)
							flag = .true.
							do l = 1, py+1
								do k = 1, px+1
									! (will trigger VLS-optimized vector load by IFX compiler)
									rho(k,l) = u_con_(1,k,l)
									prs(k,l) = _macro_prs_(u_con_(1,k,l), u_con_(2,k,l), u_con_(3,k,l), u_con_(4,k,l), u_con_(5,k,l), u_con_(6,k,l), u_con_(7,k,l), u_con_(8,k,l), gamma)
									flag = flag .and. (rho(k,l)>=PP_eps) .and. (prs(k,l)>=PP_eps)
								end do
							end do
							if (flag) then
								cell_ent_isnan_ = .false.
								cell_ent_ = 0.0
								do l = 1, py+1
									do k = 1, px+1
										! guaranteed to be legal, because (rho(k,l)>=PP_eps) .and. (prs(k,l)>=PP_eps) is just checked above
										cell_ent_ = cell_ent_ + (0.25*hx*hy*qw_x(k)*qw_y(l)) * _macro_ent_(rho(k,l), prs(k,l), gamma)
									end do
								end do
							end if
						end block
					end if
					
					! update global average
					if ((i>=1.and.i<=nx) .and. (j>=1.and.j<=ny)) then
						u_avg_global(:) = u_avg_global(:) + (hx*hy) * u_avg_(:)
					end if
					
					if (apply_LDF) then
					! if (apply_LDF .and. (grid_y(j)>=0.1)) then !DEBUG: not applied to bottom cells of a certain level
					block
						real :: mag_val(px+1, py+1, 2), mag_coef(px+1, py+1, 2)
						! LDF projection (structure also applicable to curvilinear elements)
						
						! step 1: use contravariant Piola mapping to transform to the reference element: uhat(xhat) = det(J) * J^{-1} * u(x)
						do l = 1, py+1
							do k = 1, px+1
								! det(J) * J^{-1} = diag(hy/2, hx/2)
								mag_val(k,l,1) = (0.5*hy) * (u_con_(5,k,l) - u_avg_(5))
								mag_val(k,l,2) = (0.5*hx) * (u_con_(6,k,l) - u_avg_(6))
							end do
						end do
						
						! step 2: transform to Legendre modal coefficients
						mag_coef(:,:,1) = matmul(matmul(Pnm_x, mag_val(:,:,1)), transpose(Pnm_y))
						mag_coef(:,:,2) = matmul(matmul(Pnm_x, mag_val(:,:,2)), transpose(Pnm_y))
						
						! step 3: LDF projection
						call LDFproj_%apply(mag_coef)
						
						! step 4: transform back to Lobatto nodal values
						mag_val(:,:,1) = matmul(matmul(Pmn_x, mag_coef(:,:,1)), transpose(Pmn_y))
						mag_val(:,:,2) = matmul(matmul(Pmn_x, mag_coef(:,:,2)), transpose(Pmn_y))
						
						! step 5: use contravariant Piola mapping to transform to the physical element: u(x) = 1/det(J) * J * uhat(xhat)
						do l = 1, py+1
							do k = 1, px+1
								! 1/det(J) * J = diag(2/hy, 2/hx)
								u_con_(5,k,l) = (2.0/hy) * mag_val(k,l,1) + u_avg_(5)
								u_con_(6,k,l) = (2.0/hx) * mag_val(k,l,2) + u_avg_(6)
							end do
						end do
					end block
					end if
					
					! write
					cell_ent_FV(i,j) = ent_avg
					cell_ent_DG(i,j) = cell_ent_
					cell_ent_nan(i,j) = cell_ent_isnan_
					u_con(:,:,:,i,j) = u_con_
					u_avg(:,i,j) = u_avg_
				end block
				end if
			end do
		end do
		!$omp end do
		! implicit barrier here
		
		if (apply_oedg) then
			!$omp do collapse(2) schedule(static) reduction(max: u_vari_global) ! Reduction for Fortran arrays is supported on OpenMP 4.5 standard.
			do j = 1, ny
				do i = 1, nx
					do l = 1, py+1
						do k = 1, px+1
							do s = 1, 8
								u_vari_global(s) = max(u_vari_global(s), abs(u_con(s,k,l,i,j) - u_avg_global(s)))
							end do
						end do
					end do
				end do
			end do
			!$omp end do nowait
			!$omp do collapse(2) schedule(static)
			do j = 1, ny
				do i = 1, nx+1
					! x-face between cells (i-1,j) and (i,j)
					block
					real, dimension(8,px+1,py+1) :: u_con_this, u_con_left
					integer :: kk, ll, mm, dd
					integer :: fac_kk, fac_ll ! factorial number
					real :: coef
					real :: hx_pow_kk, hx_pow ! hx to some powers
					real, dimension(px+1, py+1) :: dudx_ll, dudx_rr
					real, dimension(py+1) :: dudxdy_jump
					real :: total_jump_L1, dudxdy_jump_L1
					
					u_con_this = u_con(:,:,:,i  ,j)
					u_con_left = u_con(:,:,:,i-1,j)
					
					do s = 1, 8
						!<-------------------------------------------------------------------------------------------------------->!
						! full version / more dissipative than Kailiang Wu's version
						total_jump_L1 = 0.0
						
						dudx_ll = u_con_left(s,:,:)
						dudx_rr = u_con_this(s,:,:)
						fac_kk = 1
						hx_pow_kk = 1.0
						do kk = 0, px
							dudxdy_jump = dudx_rr(1   ,:) - dudx_ll(px+1,:)
							fac_ll = 1
							hx_pow = hx_pow_kk
							do ll = 0, max(px,py)-kk
								mm = kk + ll ! total order of derivative
								dudxdy_jump_L1 = dot_product(0.5*qw_y, abs(dudxdy_jump))
								coef = (((2*mm+1)*hx_pow)/(2*max(2*px-1,1)*fac_kk*fac_ll))
								total_jump_L1 = total_jump_L1 + coef * dudxdy_jump_L1
								if (ll < py) then
									dudxdy_jump = matmul((2.0/hy)*diff_y, dudxdy_jump)
									fac_ll = fac_ll * (ll+1)
									hx_pow = hx_pow * hx
								end if
							end do
							if (kk < px) then
								dudx_ll = matmul((2.0/hx)*diff_x, dudx_ll)
								dudx_rr = matmul((2.0/hx)*diff_x, dudx_rr)
								fac_kk = fac_kk * (kk+1)
								hx_pow_kk = hx_pow_kk * hx
							end if
						end do
						OEDG_sigsum_x(s,i,j) = total_jump_L1
						
						!<-------------------------------------------------------------------------------------------------------->!
						
						! simplified version / Yong Liu's version (less dissipative) (only jump in value and gradient)
						OEDG_sigma_x(s,i,j) = dot_product(0.5*qw_y, &
							&   1.0*abs(u_con_this(s,1,:) - u_con_left(s,px+1,:)) &
							& + (3*hx)*abs((2.0/hx)*matmul(diff_x(1,:), u_con_this(s,:,:)) - (2.0/hx)*matmul(diff_x(px+1,:), u_con_left(s,:,:))) &
							& + (3*hx)*abs((2.0/hy)*matmul(diff_y(:,:), u_con_this(s,1,:) - u_con_left(s,px+1,:))) &
							& ) / (2*max(2*px-1,1))
						!<-------------------------------------------------------------------------------------------------------->!
					end do
					end block
				end do
			end do
			!$omp end do nowait
			!$omp do collapse(2) schedule(static)
			do i = 1, nx
				do j = 1, ny+1
					! y-face between cells (i,j-1) and (i,j)
					block
					real, dimension(8,px+1,py+1) :: u_con_this, u_con_bott
					integer :: kk, ll, mm, dd
					integer :: fac_kk, fac_ll ! factorial number
					real :: coef
					real :: hy_pow_ll, hy_pow ! hy to some powers
					real, dimension(px+1, py+1) :: dudy_dd, dudy_uu
					real, dimension(px+1) :: dudxdy_jump
					real :: total_jump_L1, dudxdy_jump_L1
					
					u_con_this = u_con(:,:,:,i,j  )
					u_con_bott = u_con(:,:,:,i,j-1)
					
					do s = 1, 8
						!<-------------------------------------------------------------------------------------------------------->!
						! full version / more dissipative than Kailiang Wu's version
						total_jump_L1 = 0.0
						
						dudy_dd = u_con_bott(s,:,:)
						dudy_uu = u_con_this(s,:,:)
						fac_ll = 1
						hy_pow_ll = 1.0
						do ll = 0, py
							dudxdy_jump = dudy_uu(:,1   ) - dudy_dd(:,py+1)
							fac_kk = 1
							hy_pow = hy_pow_ll
							do kk = 0, max(px,py)-ll
								mm = kk + ll ! total order of derivative
								dudxdy_jump_L1 = dot_product(0.5*qw_x, abs(dudxdy_jump))
								coef = (((2*mm+1)*hy_pow)/(2*max(2*py-1,1)*fac_kk*fac_ll))
								total_jump_L1 = total_jump_L1 + coef * dudxdy_jump_L1
								if (kk < px) then
									dudxdy_jump = matmul((2.0/hx)*diff_x, dudxdy_jump)
									fac_kk = fac_kk * (kk+1)
									hy_pow = hy_pow * hy
								end if
							end do
							if (ll < py) then
								dudy_dd = matmul(dudy_dd, transpose((2.0/hy)*diff_y))
								dudy_uu = matmul(dudy_uu, transpose((2.0/hy)*diff_y))
								fac_ll = fac_ll * (ll+1)
								hy_pow_ll = hy_pow_ll * hy
							end if
						end do
						OEDG_sigsum_y(s,i,j) = total_jump_L1
						
						!<-------------------------------------------------------------------------------------------------------->!
						
						! simplified version / Yong Liu's version (less dissipative) (only jump in value and gradient)
						OEDG_sigma_y(s,i,j) = dot_product(0.5*qw_x, &
							&   1.0*abs(u_con_this(s,:,1) - u_con_bott(s,:,py+1)) &
							& + (3*hy)*abs((2.0/hx)*matmul(diff_x(:,:), u_con_this(s,:,1) - u_con_bott(s,:,py+1))) &
							& + (3*hy)*abs((2.0/hy)*matmul(u_con_this(s,:,:), diff_y(1,:)) - (2.0/hy)*matmul(u_con_bott(s,:,:), diff_y(py+1,:))) &
							& ) / (2*max(2*py-1,1))
						!<-------------------------------------------------------------------------------------------------------->!
					end do
					end block
				end do
			end do
			!$omp end do nowait
			
			!$omp barrier
		end if
		
		! apply limiters and compute total entropy (only interior cells)
		!$omp do collapse(2) schedule(static) reduction(+: num_negcell, tot_ent)
		do j = 1, ny
			do i = 1, nx
				if (neg_avg) cycle ! The rest part is meaningless. Restart the RK stage as soon as possible.
				limit_this_cell: block
					real :: u_avg_(8), u_con_(8, px+1, py+1)
					real :: den_avg, prs_avg, ent_avg
					
					! fetch vars
					u_con_ = u_con(:,:,:,i,j)
					u_avg_ = u_avg(:,i,j)
					
					! first apply OEDG limiter (simplified version so that it is a convex limiter)
					! (entropy stable; positivity maintaining)
					! (This one should be applied first, because it depends on neighbour information, which may change after LDF and PP limiters.)
					! Remark: Integration factor RK methods cannot be rewritten as applying convex limiters applied to each stage.
					! Ref: OEDG: Oscillation-eliminating discontinuous Galerkin method for hyperbolic conservation laws
					! Ref: Robust DG schemes on unstructured triangular meshes: Oscillation elimination and bound preservation via optimal convex decomposition
					! Ref: An Entropy Stable Essentially Oscillation-Free Discontinuous Galerkin Method for Hyperbolic Conservation Laws
					! Ref: An entropy stable essentially oscillation-free discontinuous Galerkin method for solving ideal magnetohydrodynamic equations
					if (apply_oedg) then
					block
						integer :: kk, ll, dd
						real :: maxeig_x, maxeig_y, theta
						real, dimension(8) :: vec_delta
						real :: rho, inv_rho, mx, my, mz, vx, vy, vz, bx, by, bz, e_tot, e_knt, b_sq, p_fld, p_mag, snc_sq, alfven_sq, sum_sq, coef, sum_sqsq
						
						! compute the maximum local characteristic speed (within the cell)
						maxeig_x = 0.0
						maxeig_y = 0.0
						do l = 1, py+1
							do k = 1, px+1
								rho = max(u_con_(1,k,l), PP_eps)
								mx = u_con_(2,k,l)
								my = u_con_(3,k,l)
								mz = u_con_(4,k,l)
								bx = u_con_(5,k,l)
								by = u_con_(6,k,l)
								bz = u_con_(7,k,l)
								e_tot = u_con_(8,k,l)
								
								inv_rho = 1.0 / rho ! specific volume
								vx = inv_rho*mx ! x velocity
								vy = inv_rho*my ! y velocity
								vz = inv_rho*mz ! z velocity
								e_knt = 0.5*inv_rho*(mx**2 + my**2 + mz**2) ! kinetic energy
								b_sq = bx**2 + by**2 + bz**2 ! squared norm of magnetic field
								p_mag = 0.5*b_sq ! magnetic pressure
								p_fld = max((gamma-1.0)*(e_tot - (e_knt + p_mag)), PP_eps) ! fluid pressure (might be negative due to the previous LDF mapping)
								snc_sq = gamma * p_fld * inv_rho ! squared fluid sonic speed
								alfven_sq = b_sq * inv_rho ! squared Alfven speed
								sum_sq = snc_sq + alfven_sq
								sum_sqsq = sum_sq**2
								coef = 4*snc_sq*inv_rho
								
								maxeig_x = max(maxeig_x, abs(vx) + sqrt(0.5*max(sum_sq + sqrt(max(sum_sqsq - coef * bx**2, 0.0)), 0.0)))
								maxeig_y = max(maxeig_y, abs(vy) + sqrt(0.5*max(sum_sq + sqrt(max(sum_sqsq - coef * by**2, 0.0)), 0.0)))
							end do
						end do
						
						!<-------------------------------------------------------------------------------------------------------->!
						! simplified OEDG (only lowest-order dissipation) (globally affine-invariant) (seems not stable enough for k>=3)
						! Ref: Structure-preserving oscillation-eliminating discontinuous Galerkin schemes for ideal MHD equations Locally divergence-free and positivity-preserving
						! do s = 1, 8
						! 	if (u_vari_global(s) > OEDG_eps) then
						! 		vec_delta(s) = ( &
						! 			&   (maxeig_x / hx) * (OEDG_sigma_x(s,i,j) + OEDG_sigma_x(s,i+1,j)) &
						! 			& + (maxeig_y / hy) * (OEDG_sigma_y(s,i,j) + OEDG_sigma_y(s,i,j+1)) &
						! 			& ) / u_vari_global(s)
						! 	else
						! 		vec_delta(s) = 0.0
						! 	end if
						! end do
						
						!<-------------------------------------------------------------------------------------------------------->!
						! full OEDG (convex version) (globally affine-invariant)
						! Ref: Structure-preserving oscillation-eliminating discontinuous Galerkin schemes for ideal MHD equations Locally divergence-free and positivity-preserving
						do s = 1, 8
							if (u_vari_global(s) > OEDG_eps) then
								vec_delta(s) = ( &
									&   (maxeig_x / hx) * (OEDG_sigsum_x(s,i,j) + OEDG_sigsum_x(s,i+1,j)) &
									& + (maxeig_y / hy) * (OEDG_sigsum_y(s,i,j) + OEDG_sigsum_y(s,i,j+1)) &
									& ) / u_vari_global(s)
							else
								vec_delta(s) = 0.0
							end if
						end do
						
						!<-------------------------------------------------------------------------------------------------------->!
						theta = exp(-oedg_dt * maxval(vec_delta(:)))
						! convex blending (maintain divergence for x-y magnetic field)
						do s = 1, 8
							u_con_(s,:,:) = u_avg_(s) + theta * (u_con_(s,:,:) - u_avg_(s))
						end do
						
					end block
					end if
					
					! finally apply PP limiter (not the optimal convex decomposition (OCD))
					! (entropy stable; maintain LDF)
					block
						real :: theta, den_min, prs_min
						logical :: neg_node
						integer :: it_count
						neg_node = .false.
						
						den_avg = u_avg_(1)
						prs_avg = _macro_prs_(u_avg_(1), u_avg_(2), u_avg_(3), u_avg_(4), u_avg_(5), u_avg_(6), u_avg_(7), u_avg_(8), gamma)
						
						! limit the density (only using the density) (once is enough)
						den_min = huge(1.0)
						do l = 1, py+1
							do k = 1, px+1
								den_min = min(den_min, u_con_(1,k,l))
							end do
						end do
						if (den_min < PP_eps) then
							neg_node = .true.
							theta = max(min((den_avg - PP_eps)/(den_avg - den_min), 1.0), 0.0)
							do l = 1, py+1
								do k = 1, px+1
									u_con_(1,k,l) = u_avg_(1) + theta * (u_con_(1,k,l) - u_avg_(1))
								end do
							end do
							u_con_(1,:,:) = max(u_con_(1,:,:), PP_eps) ! can use max() because this rho is linear in conservative vars
						end if
						
						! limit the pressure (using the whole variable vector)
						prs_min = huge(1.0)
						do l = 1, py+1
							do k = 1, px+1
								prs_min = min(prs_min, _macro_prs_(u_con_(1,k,l), u_con_(2,k,l), u_con_(3,k,l), u_con_(4,k,l), u_con_(5,k,l), u_con_(6,k,l), u_con_(7,k,l), u_con_(8,k,l), gamma))
							end do
						end do
						if (prs_min < PP_eps) then
							neg_node = .true.
							theta = max(min((prs_avg - PP_eps)/(prs_avg - prs_min), 1.0), 0.0)
							do l = 1, py+1
								do k = 1, px+1
									do s = 1, 8
										u_con_(s,k,l) = u_avg_(s) + theta * (u_con_(s,k,l) - u_avg_(s))
									end do
								end do
							end do
							
							! (recurring in case of FP error) (note that the theta formula is not sharp) (maybe necessary for extreme problems like MHD jet)
							it_count = 0
							do
								prs_min = huge(1.0)
								do l = 1, py+1
									do k = 1, px+1
										prs_min = min(prs_min, _macro_prs_(u_con_(1,k,l), u_con_(2,k,l), u_con_(3,k,l), u_con_(4,k,l), u_con_(5,k,l), u_con_(6,k,l), u_con_(7,k,l), u_con_(8,k,l), gamma))
									end do
								end do
								if (prs_min < 0.99*PP_eps) then
									if(it_count>1) then
										write(stderr, *) 're-enter limiting due to FP error in pressure, min_prs=', prs_min
										error stop 123456789
									end if
									theta = 0.6*theta ! hard shrink (only when extreme FP errors occur)
									do l = 1, py+1
										do k = 1, px+1
											do s = 1, 8
												u_con_(s,k,l) = u_avg_(s) + theta * (u_con_(s,k,l) - u_avg_(s))
											end do
										end do
									end do
								else
									exit
								end if
								it_count = it_count + 1
							end do
						end if
						
						if (neg_node) num_negcell = num_negcell + 1
					end block
					
					! convex limiting based to offset potential increase in cell total entropy
					if (ent_offset .and. .not.cell_ent_nan(i,j)) then
					block
						real :: cell_ent_old, cell_ent_new, theta
						cell_ent_old = cell_ent_DG(i,j)
						ent_avg = cell_ent_FV(i,j)
						
						cell_ent_new = 0.0
						do l = 1, py+1
							do k = 1, px+1
								! Warning: not necessarily legal here, because floating point error might make density or pressure negative, even if PP limiter is applied
								! (robustness enhancement applied)
								cell_ent_new = cell_ent_new + (0.25*hx*hy*qw_x(k)*qw_y(l)) * _macro_ent_(u_con_(1,k,l), max(_macro_prs_(u_con_(1,k,l), u_con_(2,k,l), u_con_(3,k,l), u_con_(4,k,l), u_con_(5,k,l), u_con_(6,k,l), u_con_(7,k,l), u_con_(8,k,l), gamma), PP_eps), gamma)
							end do
						end do
						if (cell_ent_new > cell_ent_old) then
							! convex limiting based on cell total entropy
							! Applicability is by the convexity of (math) entropy w.r.t. conservative variables
							theta = max(min((cell_ent_old - ent_avg)/(cell_ent_new - ent_avg), 1.0), 0.0)
							do l = 1, py+1
								do k = 1, px+1
									do s = 1, 8
										u_con_(s,k,l) = u_avg_(s) + theta * (u_con_(s,k,l) - u_avg_(s))
									end do
								end do
							end do
						end if
					end block
					end if
					
					if (comput_ent) then
					block
						real :: cell_ent_final
						cell_ent_final = 0.0
						do l = 1, py+1
							do k = 1, px+1
								! (will trigger VLS-optimized vector load by IFX compiler)
								! Warning: not necessarily legal here, because floating point error might make density or pressure negative, even if PP limiter is applied
								! (robustness enhancement applied)
								cell_ent_final = cell_ent_final + (0.25*hx*hy*qw_x(k)*qw_y(l)) * _macro_ent_(u_con_(1,k,l), max(_macro_prs_(u_con_(1,k,l), u_con_(2,k,l), u_con_(3,k,l), u_con_(4,k,l), u_con_(5,k,l), u_con_(6,k,l), u_con_(7,k,l), u_con_(8,k,l), gamma), PP_eps), gamma)
							end do
						end do
						tot_ent = tot_ent + cell_ent_final
					end block
					end if
					
					! write to memory
					u_con(:,:,:,i,j) = u_con_
				end block limit_this_cell
			end do
		end do
		!$omp end do nowait
		
		!$omp end parallel
		
		if (neg_avg .or. .not.comput_ent) tot_ent = ieee_value(0.0, ieee_signaling_nan)
	end subroutine post_process
	
	subroutine DG_semidiscretize(ut_con, dt_lim, u_con)
		!! compute the time derivative
		use MHD_PDE_m, only: gamma
		implicit none
		
		real, dimension(8, px+1, py+1, 0:nx+1, 0:ny+1), intent(in) :: u_con ! with one layer of ghost cells
		
		real, dimension(8, px+1, py+1, nx, ny), intent(out) :: ut_con ! only interior cells
		real, intent(out) :: dt_lim
		
		integer :: i,j,k,l,r,s
		
		dt_lim = huge(1.0) ! dt for linear stability (using CFL=1.0 outside of this subroutine)
		
		!$omp parallel
		
		! compute numerical flux (after assigning values to ghost cells)
		!$omp do collapse(2) schedule(static)
		do j = 1, ny
			do i = 1, nx+1
				! x-face between cells (i-1,j) and (i,j)
				do l = 1, py+1
					call ESPP_HLL_flux_x(u_con(:,px+1,l,i-1,j), u_con(:,1,l,i,j), gamma, nflux_x(:,l,i,j), nB_x(l,i,j), sxl(l,i,j), sxr(l,i,j))
					
					! call PP_HLL_flux_x(u_con(:,px+1,l,i-1,j), u_con(:,1,l,i,j), gamma, nflux_x(:,l,i,j), nB_x(l,i,j), sxl(l,i,j), sxr(l,i,j))
					! call ESPP_LF_flux_x(u_con(:,px+1,l,i-1,j), u_con(:,1,l,i,j), gamma, nflux_x(:,l,i,j), nB_x(l,i,j), sxl(l,i,j), sxr(l,i,j))
				end do
			end do
		end do
		!$omp end do nowait
		!$omp do collapse(2) schedule(static)
		do j = 1, ny+1
			do i = 1, nx
				! y-face between cells (i,j-1) and (i,j)
				do k = 1, px+1
					block
					real, dimension(8) :: temp1, temp2, temp_flux
					temp1 = swap_xy(u_con(:,k,py+1,i,j-1))
					temp2 = swap_xy(u_con(:,k,1,i,j))
					
					call ESPP_HLL_flux_x(temp1, temp2, gamma, temp_flux, nB_y(k,i,j), syd(k,i,j), syu(k,i,j))
					
					! call PP_HLL_flux_x(temp1, temp2, gamma, temp_flux, nB_y(k,i,j), syd(k,i,j), syu(k,i,j))
					! call ESPP_LF_flux_x(temp1, temp2, gamma, temp_flux, nB_y(k,i,j), syd(k,i,j), syu(k,i,j))
					!
					nflux_y(:,k,i,j) = swap_xy(temp_flux)
					end block
				end do
			end do
		end do
		!$omp end do nowait
		
		!$omp barrier
		
		!$omp do collapse(2) schedule(static) reduction(min: dt_lim)
		do j = 1, ny
			do i = 1, nx
				! for (i,j)-element
				block
					real :: max_sigma_x, max_sigma_y
					real, dimension(px+1, py+1) :: rho, bx, by, bz, vx, vy, vz, beta, divB
					real, dimension(8, px+1, py+1) :: ut_con_, GPcoefs
					real, dimension(8, px+1, py+1, 2) :: flux
					!---------------------------------------------------------------------------------------------------!
					real :: mx, my, mz, e_tot, ex, ey, ez, inv_rho, e_knt, b_sq, p_mag, p_fld, p_tot, b_dot_u, snc_sq, alfven_sq, sum_sq, coef
					
					! manually conduct loop fusion and loop-carried scalar replacement optimizations: separate target and temporary variables
					! Target variables (to feed the flux difference scheme) remain vectors, while temporary variables are replaced with scalars to fit into register.
					! 10 unmasked unaligned unit stride store: rho, bx, by, bz, vx, vy, vz, beta, sigma_x, sigma_y
					! 8 unmasked VLS-optimized loads
					! 24 unmasked VLS-optimized stores
					max_sigma_x = 0.0
					max_sigma_y = 0.0
					do l = 1, py+1
						do k = 1, px+1
							! fetch conservative variables (will enable VLS-optimized vector load by IFX compiler)
							rho(k,l) = u_con(1,k,l,i,j) ! density (target)
							mx = u_con(2,k,l,i,j) ! x momentum (aux)
							my = u_con(3,k,l,i,j) ! y momentum (aux)
							mz = u_con(4,k,l,i,j) ! z momentum (aux)
							bx(k,l) = u_con(5,k,l,i,j) ! x magnetic field (target)
							by(k,l) = u_con(6,k,l,i,j) ! y magnetic field (target)
							bz(k,l) = u_con(7,k,l,i,j) ! z magnetic field (target)
							e_tot = u_con(8,k,l,i,j) ! total energy (aux)
							
							! intermediate variables (will be transformed into fused loops by IFX)
							inv_rho = 1.0 / rho(k,l) ! specific volume (aux)
							vx(k,l) = inv_rho*mx ! x velocity (target)
							vy(k,l) = inv_rho*my ! y velocity (target)
							vz(k,l) = inv_rho*mz ! z velocity (target)
							e_knt = 0.5*inv_rho*(mx**2 + my**2 + mz**2) ! kinetic energy (aux)
							ex = by(k,l)*vz(k,l) - bz(k,l)*vy(k,l) ! x electric field (aux)
							ey = bz(k,l)*vx(k,l) - bx(k,l)*vz(k,l) ! y electric field (aux)
							ez = bx(k,l)*vy(k,l) - by(k,l)*vx(k,l) ! z electric field (aux)
							b_sq = bx(k,l)**2 + by(k,l)**2 + bz(k,l)**2 ! squared norm of magnetic field (aux)
							p_mag = 0.5*b_sq ! magnetic pressure (aux)
							p_fld = max((gamma-1.0)*(e_tot - (e_knt + p_mag)), PP_eps) ! fluid pressure (aux) (robust implementation with PP limiters)
							beta(k,l) = 0.5 * rho(k,l) / p_fld ! (target)
							p_tot = p_fld + p_mag ! total pressure (aux)
							b_dot_u = bx(k,l)*vx(k,l) + by(k,l)*vy(k,l) + bz(k,l)*vz(k,l) ! inner product of magnetic field and velocity (aux)
							snc_sq = gamma * p_fld * inv_rho ! squared fluid sonic speed (aux)
							alfven_sq = b_sq * inv_rho ! squared Alfven speed (aux)
							sum_sq = snc_sq + alfven_sq ! (aux)
							coef = 4*snc_sq*inv_rho ! (aux)
							
							! maximum spectral radius of the Jacobian
							max_sigma_x = max(max_sigma_x, abs(vx(k,l)) + sqrt(0.5*max(sum_sq + sqrt(max(sum_sq**2 - coef * bx(k,l)**2, 0.0)), 0.0)))
							max_sigma_y = max(max_sigma_y, abs(vy(k,l)) + sqrt(0.5*max(sum_sq + sqrt(max(sum_sq**2 - coef * by(k,l)**2, 0.0)), 0.0)))
							
							! compute flux function (will enable 2 VLS-optimized vector stores by IFX compiler)
							flux(1,k,l,1) = mx ! density flux
							flux(2,k,l,1) = vx(k,l)*mx - bx(k,l)*bx(k,l) + p_tot ! mx flux
							flux(3,k,l,1) = vx(k,l)*my - by(k,l)*bx(k,l) + 0.0 ! my flux
							flux(4,k,l,1) = vx(k,l)*mz - bz(k,l)*bx(k,l) + 0.0 ! mz flux
							flux(5,k,l,1) = 0.0 ! bx flux
							flux(6,k,l,1) = -ez ! by flux
							flux(7,k,l,1) = ey ! bz flux
							flux(8,k,l,1) = vx(k,l)*(e_tot + p_tot) - bx(k,l)*b_dot_u ! e_tot flux
							!
							flux(1,k,l,2) = my ! density flux
							flux(2,k,l,2) = vy(k,l)*mx - bx(k,l)*by(k,l) + 0.0 ! mx flux
							flux(3,k,l,2) = vy(k,l)*my - by(k,l)*by(k,l) + p_tot ! my flux
							flux(4,k,l,2) = vy(k,l)*mz - bz(k,l)*by(k,l) + 0.0 ! mz flux
							flux(5,k,l,2) = ez ! bx flux
							flux(6,k,l,2) = 0.0 ! by flux
							flux(7,k,l,2) = -ex ! bz flux
							flux(8,k,l,2) = vy(k,l)*(e_tot + p_tot) - by(k,l)*b_dot_u ! e_tot flux
							
							! Godunov--Powell coefficients (will enable VLS-optimized vector store by IFX compiler)
							GPcoefs(1,k,l) = 0.0
							GPcoefs(2,k,l) = bx(k,l)
							GPcoefs(3,k,l) = by(k,l)
							GPcoefs(4,k,l) = bz(k,l)
							GPcoefs(5,k,l) = vx(k,l)
							GPcoefs(6,k,l) = vy(k,l)
							GPcoefs(7,k,l) = vz(k,l)
							GPcoefs(8,k,l) = b_dot_u
						end do
					end do
					
					!---------------------------------------------------------------------------------------------------!
					
					! dt for DG linear stability (only considering cell-wise nodal contribution; no inter-cell effects; order of polynomial considered)
					dt_lim = min(dt_lim, 1.0/( (2*px+1)*max(max_sigma_x, maxval(sxr(:,i,j)), -minval(sxl(:,i+1,j)))/hx + (2*py+1)*max(max_sigma_y, maxval(syu(:,i,j)), -minval(syd(:,i,j+1)))/hy ))
					
					! initialize local time derivative (will be converted to memset call by IFX compiler)
					ut_con_ = 0.0
					
					! compute local divergence of magnetic field, and add the Godunov-Powell source term to the time derivative
					! (can be avoided if LDF projection is applied before this step)
					divB = (2.0/hx)*matmul(diff_x, bx) + (2.0/hy)*matmul(by, transpose(diff_y))
					
					!<-------------------------------------------------------------------------------------------------------->!
					! Legendre-Gauss-Lobatto flux difference
					
					! formulation in x: ut(k,l) -= 2*(2/hx)*D(k,r)*Fc((k,l),(r,l))
					! find Fc((k,l),(k,l))=F(k,l), then broadcast to ut(k,l)
					! find Fc((k,l),(r,l)), then broadcast to ut(k,l) and ut(r,l) (for k neq r)
					
					! formulation in y: ut(k,l) -= 2*(2/hy)*D(l,s)*Fc((k,l),(k,s))
					! find F(k,l), then broadcast to ut(k,l)
					! find Fc((k,l),(k,s)), then broadcast to ut(k,l) and ut(k,s) (for l neq s)
					
					! plus: Godunov-Powell volume source term (in flux differencing form; belong to the ESDG formulation)
					
					block
					real, dimension(8) :: ecflx_buf
					do l = 1, py+1
						do k = 1, px+1
							ut_con_(:,k,l) = ut_con_(:,k,l) - 2.0*(2.0/hx)*diff_x(k,k)*flux(:,k,l,1) - 2.0*(2.0/hy)*diff_y(l,l)*flux(:,k,l,2) - divB(k,l)*GPcoefs(:,k,l)
							do r = 1, k-1 ! (avoid redundancy)
								ecflx_buf = ec_flux_x( &
									& rho(k,l), vx(k,l), vy(k,l), vz(k,l), bx(k,l), by(k,l), bz(k,l), beta(k,l), &
									& rho(r,l), vx(r,l), vy(r,l), vz(r,l), bx(r,l), by(r,l), bz(r,l), beta(r,l), &
									& gamma)
								ut_con_(:,k,l) = ut_con_(:,k,l) - 2.0*(2.0/hx)*diff_x(k,r)*ecflx_buf ! unrolled by IFX
								ut_con_(:,r,l) = ut_con_(:,r,l) - 2.0*(2.0/hx)*diff_x(r,k)*ecflx_buf ! unrolled by IFX
							end do
							do s = 1, l-1 ! (avoid redundancy)
								ecflx_buf = ec_flux_y( &
									& rho(k,l), vx(k,l), vy(k,l), vz(k,l), bx(k,l), by(k,l), bz(k,l), beta(k,l), &
									& rho(k,s), vx(k,s), vy(k,s), vz(k,s), bx(k,s), by(k,s), bz(k,s), beta(k,s), &
									& gamma)
								ut_con_(:,k,l) = ut_con_(:,k,l) - 2.0*(2.0/hy)*diff_y(l,s)*ecflx_buf ! unrolled by IFX
								ut_con_(:,k,s) = ut_con_(:,k,s) - 2.0*(2.0/hy)*diff_y(s,l)*ecflx_buf ! unrolled by IFX
							end do
						end do
					end do
					end block
					
					!<-------------------------------------------------------------------------------------------------------->!
					! weak-form volume terms
					
					! do l = 1, py+1
					! 	do k = 1, px+1
					! 		do r = 1, px+1
					! 			ut_con_(:,k,l) = ut_con_(:,k,l) - (2.0/hx)*diff_x(k,r)*flux(:,r,l,1)
					! 		end do
					! 		do s = 1, py+1
					! 			ut_con_(:,k,l) = ut_con_(:,k,l) - (2.0/hy)*diff_y(l,s)*flux(:,k,s,2)
					! 		end do
					! 		ut_con_(:,k,l) = ut_con_(:,k,l) - divB(k,l)*GPcoefs(:,k,l)
					! 	end do
					! end do
					
					!<-------------------------------------------------------------------------------------------------------->!
					
					! x-face term using numerical flux and normal magnetic field jump (will be loop fused and vectorized by IFX compiler)
					! Warning: may have race condition when px=0!
					do l = 1, py+1
						ut_con_(:,   1,l) = ut_con_(:,   1,l) - ((2.0/hx)/qw_x(   1)) * (flux(:,   1,l,1) - nflux_x(:,l,  i,j) + GPcoefs(:,   1,l)*(bx(   1,l) - nB_x(l,  i,j)))
						ut_con_(:,px+1,l) = ut_con_(:,px+1,l) + ((2.0/hx)/qw_x(px+1)) * (flux(:,px+1,l,1) - nflux_x(:,l,i+1,j) + GPcoefs(:,px+1,l)*(bx(px+1,l) - nB_x(l,i+1,j)))
					end do
					
					! y-face term using numerical flux and normal magnetic field jump (will be loop fused and vectorized by IFX compiler)
					! Warning: may have race condition when py=0!
					do k = 1, px+1
						ut_con_(:,k,   1) = ut_con_(:,k,   1) - ((2.0/hy)/qw_y(   1)) * (flux(:,k,   1,2) - nflux_y(:,k,i,  j) + GPcoefs(:,k,   1)*(by(k,   1) - nB_y(k,i,  j)))
						ut_con_(:,k,py+1) = ut_con_(:,k,py+1) + ((2.0/hy)/qw_y(py+1)) * (flux(:,k,py+1,2) - nflux_y(:,k,i,j+1) + GPcoefs(:,k,py+1)*(by(k,py+1) - nB_y(k,i,j+1)))
					end do
					
					! write
					ut_con(:,:,:,i,j) = ut_con_
				end block
			end do
		end do
		!$omp end do nowait
		
		!$omp end parallel
	end subroutine DG_semidiscretize
	
	subroutine ESPP_HLL_flux_x(w_l, w_r, gamma, nflux, nBn, sl, sr)
		implicit none
		real, dimension(8), intent(in) :: w_l, w_r ! conservative variables
		real, intent(in) :: gamma
		real, dimension(8), intent(out) :: nflux ! numerical flux in x-direction
		real, intent(out) :: nBn ! numerical normal magnetic field in x-direction
		real, intent(out) :: sl, sr ! signal speed estimates (sl <= 0.0 <= sr)
		
		real :: phi_l, phi_r, ent_flux_l, ent_flux_r, ent_fun_l, ent_fun_r, cfx_l, cfx_r
		real, dimension(8) :: entv_l, entv_r, flux_l, flux_r
		
		real :: rho_l, rho_r, mx_l, mx_r, my_l, my_r, mz_l, mz_r, bx_l, bx_r, by_l, by_r, bz_l, bz_r, e_tot_l, e_tot_r
		real :: inv_rho_l, inv_rho_r, vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, ex_l, ex_r, ey_l, ey_r, ez_l, ez_r, e_knt_l, e_knt_r, b_sq_l, b_sq_r, p_mag_l, p_mag_r, p_fld_l, p_fld_r, p_tot_l, p_tot_r, b_dot_u_l, b_dot_u_r, snc_sq_l, snc_sq_r, alfven_sq_l, alfven_sq_r, sum_sq_l, sum_sq_r, inv_theta_l, inv_theta_r, sp_ent_l, sp_ent_r
		real :: sqrt_rho_l, sqrt_rho_r, sum_sqrt_rho, db_aux, vx_roe
		real :: a_form, bl_form, br_form, al_form, ar_form, aa_form
		real :: sl_pp, sr_pp, sl_es, sr_es, dslr, wei_l, wei_r
		
		! fetch vars
		rho_l = w_l(1); mx_l = w_l(2); my_l = w_l(3); mz_l = w_l(4); bx_l = w_l(5); by_l = w_l(6); bz_l = w_l(7); e_tot_l = w_l(8)
		rho_r = w_r(1); mx_r = w_r(2); my_r = w_r(3); mz_r = w_r(4); bx_r = w_r(5); by_r = w_r(6); bz_r = w_r(7); e_tot_r = w_r(8)
				
		! specific volume
		inv_rho_l = 1.0 / rho_l
		inv_rho_r = 1.0 / rho_r
		
		! x velocity
		vx_l = inv_rho_l * mx_l
		vx_r = inv_rho_r * mx_r
		
		! y velocity
		vy_l = inv_rho_l * my_l
		vy_r = inv_rho_r * my_r
		
		! z velocity
		vz_l = inv_rho_l * mz_l
		vz_r = inv_rho_r * mz_r
		
		! kinetic energy
		e_knt_l = 0.5*inv_rho_l*(mx_l**2 + my_l**2 + mz_l**2)
		e_knt_r = 0.5*inv_rho_r*(mx_r**2 + my_r**2 + mz_r**2)
		
		! x electric field
		ex_l = by_l*vz_l - bz_l*vy_l
		ex_r = by_r*vz_r - bz_r*vy_r
		
		! y electric field
		ey_l = bz_l*vx_l - bx_l*vz_l
		ey_r = bz_r*vx_r - bx_r*vz_r
		
		! z electric field
		ez_l = bx_l*vy_l - by_l*vx_l
		ez_r = bx_r*vy_r - by_r*vx_r
		
		! squared norm of magnetic field
		b_sq_l = bx_l**2 + by_l**2 + bz_l**2
		b_sq_r = bx_r**2 + by_r**2 + bz_r**2
		
		! magnetic pressure
		p_mag_l = 0.5*b_sq_l
		p_mag_r = 0.5*b_sq_r
		
		! fluid pressure (robust implementation)
		p_fld_l = max((gamma-1.0)*(e_tot_l - (e_knt_l + p_mag_l)), PP_eps)
		p_fld_r = max((gamma-1.0)*(e_tot_r - (e_knt_r + p_mag_r)), PP_eps)
		
		! total pressure
		p_tot_l = p_fld_l + p_mag_l
		p_tot_r = p_fld_r + p_mag_r
		
		! dot product of magnetic field and velocity
		b_dot_u_l = bx_l*vx_l + by_l*vy_l + bz_l*vz_l
		b_dot_u_r = bx_r*vx_r + by_r*vy_r + bz_r*vz_r
		
		! fluid sonic speed squared
		snc_sq_l = gamma * p_fld_l * inv_rho_l
		snc_sq_r = gamma * p_fld_r * inv_rho_r
		
		! Alfven speed squared
		alfven_sq_l = b_sq_l * inv_rho_l
		alfven_sq_r = b_sq_r * inv_rho_r
		
		! sum of squares
		sum_sq_l = snc_sq_l + alfven_sq_l
		sum_sq_r = snc_sq_r + alfven_sq_r
		
		! fast magneto-acoustic wave speed
		cfx_l = sqrt(0.5*(sum_sq_l + sqrt(max(sum_sq_l**2 - 4*snc_sq_l*inv_rho_l * bx_l**2, 0.0))))
		cfx_r = sqrt(0.5*(sum_sq_r + sqrt(max(sum_sq_r**2 - 4*snc_sq_r*inv_rho_r * bx_r**2, 0.0))))
		
		! inverse of temperature
		inv_theta_l = rho_l/p_fld_l
		inv_theta_r = rho_r/p_fld_r
		
		phi_l = inv_theta_l*b_dot_u_l
		phi_r = inv_theta_r*b_dot_u_r
		
		! (physical) specific entropy (costly) (robustness enhancement)
		sp_ent_l = log(p_fld_l) - gamma*log(rho_l)
		sp_ent_r = log(p_fld_r) - gamma*log(rho_r)
		
		! (math) entropy function
		ent_fun_l = -rho_l * sp_ent_l / (gamma - 1.0)
		ent_fun_r = -rho_r * sp_ent_r / (gamma - 1.0)
		
		! entropy flux function
		ent_flux_l = ent_fun_l * vx_l
		ent_flux_r = ent_fun_r * vx_r
		
		! compute flux: left
		flux_l(1) = mx_l
		flux_l(2) = vx_l*mx_l - bx_l*bx_l + p_tot_l
		flux_l(3) = vx_l*my_l - by_l*bx_l
		flux_l(4) = vx_l*mz_l - bz_l*bx_l
		flux_l(5) = 0.0
		flux_l(6) = -ez_l
		flux_l(7) = ey_l
		flux_l(8) = vx_l*(e_tot_l + p_tot_l) - bx_l*b_dot_u_l
		
		! compute flux: right
		flux_r(1) = mx_r
		flux_r(2) = vx_r*mx_r - bx_r*bx_r + p_tot_r
		flux_r(3) = vx_r*my_r - by_r*bx_r
		flux_r(4) = vx_r*mz_r - bz_r*bx_r
		flux_r(5) = 0.0
		flux_r(6) = -ez_r
		flux_r(7) = ey_r
		flux_r(8) = vx_r*(e_tot_r + p_tot_r) - bx_r*b_dot_u_r
		
		! entropy variables: left
		entv_l(1) = (gamma - sp_ent_l) / (gamma - 1.0) - e_knt_l / p_fld_l
		entv_l(2) = mx_l/p_fld_l
		entv_l(3) = my_l/p_fld_l
		entv_l(4) = mz_l/p_fld_l
		entv_l(5) = inv_theta_l*bx_l
		entv_l(6) = inv_theta_l*by_l
		entv_l(7) = inv_theta_l*bz_l
		entv_l(8) = -inv_theta_l
		
		! entropy variables: right
		entv_r(1) = (gamma - sp_ent_r) / (gamma - 1.0) - e_knt_r / p_fld_r
		entv_r(2) = mx_r/p_fld_r
		entv_r(3) = my_r/p_fld_r
		entv_r(4) = mz_r/p_fld_r
		entv_r(5) = inv_theta_r*bx_r
		entv_r(6) = inv_theta_r*by_r
		entv_r(7) = inv_theta_r*bz_r
		entv_r(8) = -inv_theta_r
		
		! PP signal speed estimates
		sqrt_rho_l = sqrt(rho_l)
		sqrt_rho_r = sqrt(rho_r)
		sum_sqrt_rho = sqrt_rho_l + sqrt_rho_r
		vx_roe = (vx_l*sqrt_rho_l + vx_r*sqrt_rho_r) / sum_sqrt_rho
		db_aux = sqrt((bx_l - bx_r)**2 + (by_l - by_r)**2 + (bz_l - bz_r)**2) / sum_sqrt_rho
		sl_pp = min(vx_l, vx_roe) - cfx_l - db_aux
		sr_pp = max(vx_r, vx_roe) + cfx_r + db_aux
		
		! ES signal speed estimates
		a_form = max(dot_product(entv_r - entv_l, w_r - w_l), 0.0) + 1e-8 ! A simple cure is applied (depending on floating point precision).
		bl_form = (ent_flux_r - ent_flux_l) - dot_product(entv_l, flux_r - flux_l) - phi_l*(bx_r - bx_l)
		br_form = (ent_flux_r - ent_flux_l) - dot_product(entv_r, flux_r - flux_l) - phi_r*(bx_r - bx_l)
		al_form = max(bl_form / a_form, 0.0)
		ar_form = max(br_form / a_form, 0.0)
		aa_form = sqrt(al_form * ar_form)
		sl_es = -(ar_form + aa_form) ! checked.
		sr_es = +(al_form + aa_form) ! checked.
		
		! PP and semidiscrete ES estimate
		sl = min(sl_es, sl_pp, 0.0)
		sr = max(sr_es, sr_pp, 0.0)
		
		! HLL rule
		dslr = sr - sl
		wei_l = sr/dslr
		wei_r = -sl/dslr
		nflux = wei_l * flux_l + wei_r * flux_r + ((sl*sr)/dslr) * (w_r - w_l)
		nBn = wei_l * bx_l + wei_r * bx_r
	end subroutine ESPP_HLL_flux_x
	
	subroutine PP_HLL_flux_x(w_l, w_r, gamma, nflux, nBn, sl, sr)
		!! Kailiang Wu's PP-HLL GP-MHD flux
		implicit none
		real, dimension(8), intent(in) :: w_l, w_r ! conservative variables
		real, intent(in) :: gamma
		real, dimension(8), intent(out) :: nflux ! numerical flux in x-direction
		real, intent(out) :: nBn ! numerical normal magnetic field in x-direction
		real, intent(out) :: sl, sr ! signal speed estimates (sl <= 0.0 <= sr)
		
		real :: phi_l, phi_r, ent_flux_l, ent_flux_r, ent_fun_l, ent_fun_r, cfx_l, cfx_r
		real, dimension(8) :: entv_l, entv_r, flux_l, flux_r
		
		real :: rho_l, rho_r, mx_l, mx_r, my_l, my_r, mz_l, mz_r, bx_l, bx_r, by_l, by_r, bz_l, bz_r, e_tot_l, e_tot_r
		real :: inv_rho_l, inv_rho_r, vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, ex_l, ex_r, ey_l, ey_r, ez_l, ez_r, e_knt_l, e_knt_r, b_sq_l, b_sq_r, p_mag_l, p_mag_r, p_fld_l, p_fld_r, p_tot_l, p_tot_r, b_dot_u_l, b_dot_u_r, snc_sq_l, snc_sq_r, alfven_sq_l, alfven_sq_r, sum_sq_l, sum_sq_r, inv_theta_l, inv_theta_r, sp_ent_l, sp_ent_r
		real :: sqrt_rho_l, sqrt_rho_r, sum_sqrt_rho, db_aux, vx_roe
		real :: a_form, bl_form, br_form, al_form, ar_form, aa_form
		real :: sl_pp, sr_pp, sl_es, sr_es, dslr, wei_l, wei_r
		
		! fetch vars
		rho_l = w_l(1); mx_l = w_l(2); my_l = w_l(3); mz_l = w_l(4); bx_l = w_l(5); by_l = w_l(6); bz_l = w_l(7); e_tot_l = w_l(8)
		rho_r = w_r(1); mx_r = w_r(2); my_r = w_r(3); mz_r = w_r(4); bx_r = w_r(5); by_r = w_r(6); bz_r = w_r(7); e_tot_r = w_r(8)
				
		! specific volume
		inv_rho_l = 1.0 / rho_l
		inv_rho_r = 1.0 / rho_r
		
		! x velocity
		vx_l = inv_rho_l * mx_l
		vx_r = inv_rho_r * mx_r
		
		! y velocity
		vy_l = inv_rho_l * my_l
		vy_r = inv_rho_r * my_r
		
		! z velocity
		vz_l = inv_rho_l * mz_l
		vz_r = inv_rho_r * mz_r
		
		! kinetic energy
		e_knt_l = 0.5*inv_rho_l*(mx_l**2 + my_l**2 + mz_l**2)
		e_knt_r = 0.5*inv_rho_r*(mx_r**2 + my_r**2 + mz_r**2)
		
		! x electric field
		ex_l = by_l*vz_l - bz_l*vy_l
		ex_r = by_r*vz_r - bz_r*vy_r
		
		! y electric field
		ey_l = bz_l*vx_l - bx_l*vz_l
		ey_r = bz_r*vx_r - bx_r*vz_r
		
		! z electric field
		ez_l = bx_l*vy_l - by_l*vx_l
		ez_r = bx_r*vy_r - by_r*vx_r
		
		! squared norm of magnetic field
		b_sq_l = bx_l**2 + by_l**2 + bz_l**2
		b_sq_r = bx_r**2 + by_r**2 + bz_r**2
		
		! magnetic pressure
		p_mag_l = 0.5*b_sq_l
		p_mag_r = 0.5*b_sq_r
		
		! fluid pressure (robust implementation)
		p_fld_l = max((gamma-1.0)*(e_tot_l - (e_knt_l + p_mag_l)), PP_eps)
		p_fld_r = max((gamma-1.0)*(e_tot_r - (e_knt_r + p_mag_r)), PP_eps)
		
		! total pressure
		p_tot_l = p_fld_l + p_mag_l
		p_tot_r = p_fld_r + p_mag_r
		
		! dot product of magnetic field and velocity
		b_dot_u_l = bx_l*vx_l + by_l*vy_l + bz_l*vz_l
		b_dot_u_r = bx_r*vx_r + by_r*vy_r + bz_r*vz_r
		
		! fluid sonic speed squared
		snc_sq_l = gamma * p_fld_l * inv_rho_l
		snc_sq_r = gamma * p_fld_r * inv_rho_r
		
		! Alfven speed squared
		alfven_sq_l = b_sq_l * inv_rho_l
		alfven_sq_r = b_sq_r * inv_rho_r
		
		! sum of squares
		sum_sq_l = snc_sq_l + alfven_sq_l
		sum_sq_r = snc_sq_r + alfven_sq_r
		
		! fast magneto-acoustic wave speed
		cfx_l = sqrt(0.5*(sum_sq_l + sqrt(max(sum_sq_l**2 - 4*snc_sq_l*inv_rho_l * bx_l**2, 0.0))))
		cfx_r = sqrt(0.5*(sum_sq_r + sqrt(max(sum_sq_r**2 - 4*snc_sq_r*inv_rho_r * bx_r**2, 0.0))))
		
		! inverse of temperature
		inv_theta_l = rho_l/p_fld_l
		inv_theta_r = rho_r/p_fld_r
		
		phi_l = inv_theta_l*b_dot_u_l
		phi_r = inv_theta_r*b_dot_u_r
		
		! (physical) specific entropy (costly) (robustness enhancement)
		sp_ent_l = log(p_fld_l) - gamma*log(rho_l)
		sp_ent_r = log(p_fld_r) - gamma*log(rho_r)
		
		! (math) entropy function
		ent_fun_l = -rho_l * sp_ent_l / (gamma - 1.0)
		ent_fun_r = -rho_r * sp_ent_r / (gamma - 1.0)
		
		! entropy flux function
		ent_flux_l = ent_fun_l * vx_l
		ent_flux_r = ent_fun_r * vx_r
		
		! compute flux: left
		flux_l(1) = mx_l
		flux_l(2) = vx_l*mx_l - bx_l*bx_l + p_tot_l
		flux_l(3) = vx_l*my_l - by_l*bx_l
		flux_l(4) = vx_l*mz_l - bz_l*bx_l
		flux_l(5) = 0.0
		flux_l(6) = -ez_l
		flux_l(7) = ey_l
		flux_l(8) = vx_l*(e_tot_l + p_tot_l) - bx_l*b_dot_u_l
		
		! compute flux: right
		flux_r(1) = mx_r
		flux_r(2) = vx_r*mx_r - bx_r*bx_r + p_tot_r
		flux_r(3) = vx_r*my_r - by_r*bx_r
		flux_r(4) = vx_r*mz_r - bz_r*bx_r
		flux_r(5) = 0.0
		flux_r(6) = -ez_r
		flux_r(7) = ey_r
		flux_r(8) = vx_r*(e_tot_r + p_tot_r) - bx_r*b_dot_u_r
		
		! PP signal speed estimates
		sqrt_rho_l = sqrt(rho_l)
		sqrt_rho_r = sqrt(rho_r)
		sum_sqrt_rho = sqrt_rho_l + sqrt_rho_r
		vx_roe = (vx_l*sqrt_rho_l + vx_r*sqrt_rho_r) / sum_sqrt_rho
		db_aux = sqrt((bx_l - bx_r)**2 + (by_l - by_r)**2 + (bz_l - bz_r)**2) / sum_sqrt_rho
		sl_pp = min(vx_l, vx_roe) - cfx_l - db_aux
		sr_pp = max(vx_r, vx_roe) + cfx_r + db_aux
		
		! PP and semidiscrete ES estimate
		sl = min(sl_pp, 0.0)
		sr = max(sr_pp, 0.0)
		
		! HLL rule
		dslr = sr - sl
		wei_l = sr/dslr
		wei_r = -sl/dslr
		nflux = wei_l * flux_l + wei_r * flux_r + ((sl*sr)/dslr) * (w_r - w_l)
		nBn = wei_l * bx_l + wei_r * bx_r
	end subroutine PP_HLL_flux_x
	
	subroutine ESPP_LF_flux_x(w_l, w_r, gamma, nflux, nBn, sl, sr)
		implicit none
		real, dimension(8), intent(in) :: w_l, w_r ! conservative variables
		real, intent(in) :: gamma
		real, dimension(8), intent(out) :: nflux ! numerical flux in x-direction
		real, intent(out) :: nBn ! numerical normal magnetic field in x-direction
		real, intent(out) :: sl, sr ! signal speed estimates (sl <= 0.0 <= sr)
		
		real :: phi_l, phi_r, ent_flux_l, ent_flux_r, ent_fun_l, ent_fun_r, cfx_l, cfx_r
		real, dimension(8) :: entv_l, entv_r, flux_l, flux_r
		
		real :: rho_l, rho_r, mx_l, mx_r, my_l, my_r, mz_l, mz_r, bx_l, bx_r, by_l, by_r, bz_l, bz_r, e_tot_l, e_tot_r
		real :: inv_rho_l, inv_rho_r, vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, ex_l, ex_r, ey_l, ey_r, ez_l, ez_r, e_knt_l, e_knt_r, b_sq_l, b_sq_r, p_mag_l, p_mag_r, p_fld_l, p_fld_r, p_tot_l, p_tot_r, b_dot_u_l, b_dot_u_r, snc_sq_l, snc_sq_r, alfven_sq_l, alfven_sq_r, sum_sq_l, sum_sq_r, inv_theta_l, inv_theta_r, sp_ent_l, sp_ent_r
		real :: sqrt_rho_l, sqrt_rho_r, sum_sqrt_rho, db_aux, vx_roe
		real :: a_form, bb_form
		real :: sl_pp, sr_pp, sl_es, sr_es, dslr, wei_l, wei_r
		
		! fetch vars
		rho_l = w_l(1); mx_l = w_l(2); my_l = w_l(3); mz_l = w_l(4); bx_l = w_l(5); by_l = w_l(6); bz_l = w_l(7); e_tot_l = w_l(8)
		rho_r = w_r(1); mx_r = w_r(2); my_r = w_r(3); mz_r = w_r(4); bx_r = w_r(5); by_r = w_r(6); bz_r = w_r(7); e_tot_r = w_r(8)
				
		! specific volume
		inv_rho_l = 1.0 / rho_l
		inv_rho_r = 1.0 / rho_r
		
		! x velocity
		vx_l = inv_rho_l * mx_l
		vx_r = inv_rho_r * mx_r
		
		! y velocity
		vy_l = inv_rho_l * my_l
		vy_r = inv_rho_r * my_r
		
		! z velocity
		vz_l = inv_rho_l * mz_l
		vz_r = inv_rho_r * mz_r
		
		! kinetic energy
		e_knt_l = 0.5*inv_rho_l*(mx_l**2 + my_l**2 + mz_l**2)
		e_knt_r = 0.5*inv_rho_r*(mx_r**2 + my_r**2 + mz_r**2)
		
		! x electric field
		ex_l = by_l*vz_l - bz_l*vy_l
		ex_r = by_r*vz_r - bz_r*vy_r
		
		! y electric field
		ey_l = bz_l*vx_l - bx_l*vz_l
		ey_r = bz_r*vx_r - bx_r*vz_r
		
		! z electric field
		ez_l = bx_l*vy_l - by_l*vx_l
		ez_r = bx_r*vy_r - by_r*vx_r
		
		! squared norm of magnetic field
		b_sq_l = bx_l**2 + by_l**2 + bz_l**2
		b_sq_r = bx_r**2 + by_r**2 + bz_r**2
		
		! magnetic pressure
		p_mag_l = 0.5*b_sq_l
		p_mag_r = 0.5*b_sq_r
		
		! fluid pressure (robust implementation)
		p_fld_l = max((gamma-1.0)*(e_tot_l - (e_knt_l + p_mag_l)), PP_eps)
		p_fld_r = max((gamma-1.0)*(e_tot_r - (e_knt_r + p_mag_r)), PP_eps)
		
		! total pressure
		p_tot_l = p_fld_l + p_mag_l
		p_tot_r = p_fld_r + p_mag_r
		
		! dot product of magnetic field and velocity
		b_dot_u_l = bx_l*vx_l + by_l*vy_l + bz_l*vz_l
		b_dot_u_r = bx_r*vx_r + by_r*vy_r + bz_r*vz_r
		
		! fluid sonic speed squared
		snc_sq_l = gamma * p_fld_l * inv_rho_l
		snc_sq_r = gamma * p_fld_r * inv_rho_r
		
		! Alfven speed squared
		alfven_sq_l = b_sq_l * inv_rho_l
		alfven_sq_r = b_sq_r * inv_rho_r
		
		! sum of squares
		sum_sq_l = snc_sq_l + alfven_sq_l
		sum_sq_r = snc_sq_r + alfven_sq_r
		
		! fast magneto-acoustic wave speed
		cfx_l = sqrt(0.5*(sum_sq_l + sqrt(max(sum_sq_l**2 - 4*snc_sq_l*inv_rho_l * bx_l**2, 0.0))))
		cfx_r = sqrt(0.5*(sum_sq_r + sqrt(max(sum_sq_r**2 - 4*snc_sq_r*inv_rho_r * bx_r**2, 0.0))))
		
		! inverse of temperature
		inv_theta_l = rho_l/p_fld_l
		inv_theta_r = rho_r/p_fld_r
		
		phi_l = inv_theta_l*b_dot_u_l
		phi_r = inv_theta_r*b_dot_u_r
		
		! (physical) specific entropy (costly) (robustness enhancement)
		sp_ent_l = log(p_fld_l) - gamma*log(rho_l)
		sp_ent_r = log(p_fld_r) - gamma*log(rho_r)
		
		! (math) entropy function
		ent_fun_l = -rho_l * sp_ent_l / (gamma - 1.0)
		ent_fun_r = -rho_r * sp_ent_r / (gamma - 1.0)
		
		! entropy flux function
		ent_flux_l = ent_fun_l * vx_l
		ent_flux_r = ent_fun_r * vx_r
		
		! compute flux: left
		flux_l(1) = mx_l
		flux_l(2) = vx_l*mx_l - bx_l*bx_l + p_tot_l
		flux_l(3) = vx_l*my_l - by_l*bx_l
		flux_l(4) = vx_l*mz_l - bz_l*bx_l
		flux_l(5) = 0.0
		flux_l(6) = -ez_l
		flux_l(7) = ey_l
		flux_l(8) = vx_l*(e_tot_l + p_tot_l) - bx_l*b_dot_u_l
		
		! compute flux: right
		flux_r(1) = mx_r
		flux_r(2) = vx_r*mx_r - bx_r*bx_r + p_tot_r
		flux_r(3) = vx_r*my_r - by_r*bx_r
		flux_r(4) = vx_r*mz_r - bz_r*bx_r
		flux_r(5) = 0.0
		flux_r(6) = -ez_r
		flux_r(7) = ey_r
		flux_r(8) = vx_r*(e_tot_r + p_tot_r) - bx_r*b_dot_u_r
		
		! entropy variables: left
		entv_l(1) = (gamma - sp_ent_l) / (gamma - 1.0) - e_knt_l / p_fld_l
		entv_l(2) = mx_l/p_fld_l
		entv_l(3) = my_l/p_fld_l
		entv_l(4) = mz_l/p_fld_l
		entv_l(5) = inv_theta_l*bx_l
		entv_l(6) = inv_theta_l*by_l
		entv_l(7) = inv_theta_l*bz_l
		entv_l(8) = -inv_theta_l
		
		! entropy variables: right
		entv_r(1) = (gamma - sp_ent_r) / (gamma - 1.0) - e_knt_r / p_fld_r
		entv_r(2) = mx_r/p_fld_r
		entv_r(3) = my_r/p_fld_r
		entv_r(4) = mz_r/p_fld_r
		entv_r(5) = inv_theta_r*bx_r
		entv_r(6) = inv_theta_r*by_r
		entv_r(7) = inv_theta_r*bz_r
		entv_r(8) = -inv_theta_r
		
		! PP signal speed estimates
		sqrt_rho_l = sqrt(rho_l)
		sqrt_rho_r = sqrt(rho_r)
		sum_sqrt_rho = sqrt_rho_l + sqrt_rho_r
		vx_roe = (vx_l*sqrt_rho_l + vx_r*sqrt_rho_r) / sum_sqrt_rho
		db_aux = sqrt((bx_l - bx_r)**2 + (by_l - by_r)**2 + (bz_l - bz_r)**2) / sum_sqrt_rho
		sl_pp = min(vx_l, vx_roe) - cfx_l - db_aux
		sr_pp = max(vx_r, vx_roe) + cfx_r + db_aux
		
		! LF version (K. Wu's paper)
		sr_pp = max(abs(vx_l) + cfx_l, abs(vx_r) + cfx_r, abs(vx_roe) + max(cfx_l, cfx_r)) + db_aux
		
		! ES signal speed estimates (Lax-Friedrichs type: same for left and right)
		a_form = max(dot_product(entv_r - entv_l, w_r - w_l), 0.0) + 1e-9 ! A simple cure is applied (depending on floating point precision).
		bb_form = 2.0*(ent_flux_r - ent_flux_l) - dot_product(entv_l + entv_r, flux_r - flux_l) - (phi_l + phi_r)*(bx_r - bx_l)
		sr_es = max(bb_form / a_form, 0.0)
		
		! PP and semidiscrete ES estimate (Lax-Friedrichs type: same for left and right)
		sr = max(sr_es, -sl_pp, sr_pp, 0.0)
		sl = -sr
		
		! LF rule
		nflux = 0.5 * flux_l + 0.5 * flux_r - 0.5 * sr * (w_r - w_l)
		nBn = 0.5 * bx_l + 0.5 * bx_r
	end subroutine ESPP_LF_flux_x
	
	subroutine EC_LF_flux_x(w_l, w_r, gamma, nflux, nBn, sl, sr)
		!! Yong Liu's ES flux: EC flux ples LF-type dissipation
		implicit none
		real, dimension(8), intent(in) :: w_l, w_r ! conservative variables
		real, intent(in) :: gamma
		real, dimension(8), intent(out) :: nflux ! numerical flux in x-direction
		real, intent(out) :: nBn ! numerical normal magnetic field in x-direction
		real, intent(out) :: sl, sr ! signal speed estimates (sl <= 0.0 <= sr)
		
		real :: phi_l, phi_r, ent_flux_l, ent_flux_r, ent_fun_l, ent_fun_r, cfx_l, cfx_r
		real, dimension(8) :: entv_l, entv_r, flux_l, flux_r
		
		real :: rho_l, rho_r, mx_l, mx_r, my_l, my_r, mz_l, mz_r, bx_l, bx_r, by_l, by_r, bz_l, bz_r, e_tot_l, e_tot_r, beta_l, beta_r
		real :: inv_rho_l, inv_rho_r, vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, ex_l, ex_r, ey_l, ey_r, ez_l, ez_r, e_knt_l, e_knt_r, b_sq_l, b_sq_r, p_mag_l, p_mag_r, p_fld_l, p_fld_r, p_tot_l, p_tot_r, b_dot_u_l, b_dot_u_r, snc_sq_l, snc_sq_r, alfven_sq_l, alfven_sq_r, sum_sq_l, sum_sq_r, inv_theta_l, inv_theta_r, sp_ent_l, sp_ent_r
		real :: sqrt_rho_l, sqrt_rho_r, sum_sqrt_rho, db_aux, vx_roe
		real :: a_form, bl_form, br_form, al_form, ar_form, aa_form
		real :: sl_pp, sr_pp, sl_es, sr_es, dslr, wei_l, wei_r
		
		! fetch vars
		rho_l = w_l(1); mx_l = w_l(2); my_l = w_l(3); mz_l = w_l(4); bx_l = w_l(5); by_l = w_l(6); bz_l = w_l(7); e_tot_l = w_l(8)
		rho_r = w_r(1); mx_r = w_r(2); my_r = w_r(3); mz_r = w_r(4); bx_r = w_r(5); by_r = w_r(6); bz_r = w_r(7); e_tot_r = w_r(8)
		
		! specific volume
		inv_rho_l = 1.0 / rho_l
		inv_rho_r = 1.0 / rho_r
		
		! x velocity
		vx_l = inv_rho_l * mx_l
		vx_r = inv_rho_r * mx_r
		
		! y velocity
		vy_l = inv_rho_l * my_l
		vy_r = inv_rho_r * my_r
		
		! z velocity
		vz_l = inv_rho_l * mz_l
		vz_r = inv_rho_r * mz_r
		
		! kinetic energy
		e_knt_l = 0.5*inv_rho_l*(mx_l**2 + my_l**2 + mz_l**2)
		e_knt_r = 0.5*inv_rho_r*(mx_r**2 + my_r**2 + mz_r**2)
		
		! x electric field
		ex_l = by_l*vz_l - bz_l*vy_l
		ex_r = by_r*vz_r - bz_r*vy_r
		
		! y electric field
		ey_l = bz_l*vx_l - bx_l*vz_l
		ey_r = bz_r*vx_r - bx_r*vz_r
		
		! z electric field
		ez_l = bx_l*vy_l - by_l*vx_l
		ez_r = bx_r*vy_r - by_r*vx_r
		
		! squared norm of magnetic field
		b_sq_l = bx_l**2 + by_l**2 + bz_l**2
		b_sq_r = bx_r**2 + by_r**2 + bz_r**2
		
		! magnetic pressure
		p_mag_l = 0.5*b_sq_l
		p_mag_r = 0.5*b_sq_r
		
		! fluid pressure
		p_fld_l = max((gamma-1.0)*(e_tot_l - (e_knt_l + p_mag_l)), PP_eps)
		p_fld_r = max((gamma-1.0)*(e_tot_r - (e_knt_r + p_mag_r)), PP_eps)
		
		! beta value
		beta_l = 0.5 * rho_l / p_fld_l
		beta_r = 0.5 * rho_r / p_fld_r
		
		! total pressure
		p_tot_l = p_fld_l + p_mag_l
		p_tot_r = p_fld_r + p_mag_r
		
		! dot product of magnetic field and velocity
		b_dot_u_l = bx_l*vx_l + by_l*vy_l + bz_l*vz_l
		b_dot_u_r = bx_r*vx_r + by_r*vy_r + bz_r*vz_r
		
		! fluid sonic speed squared
		snc_sq_l = gamma * p_fld_l * inv_rho_l
		snc_sq_r = gamma * p_fld_r * inv_rho_r
		
		! Alfven speed squared
		alfven_sq_l = b_sq_l * inv_rho_l
		alfven_sq_r = b_sq_r * inv_rho_r
		
		! sum of squares
		sum_sq_l = snc_sq_l + alfven_sq_l
		sum_sq_r = snc_sq_r + alfven_sq_r
		
		! fast magneto-acoustic wave speed
		cfx_l = sqrt(0.5*(sum_sq_l + sqrt(max(sum_sq_l**2 - 4*snc_sq_l*inv_rho_l * bx_l**2, 0.0))))
		cfx_r = sqrt(0.5*(sum_sq_r + sqrt(max(sum_sq_r**2 - 4*snc_sq_r*inv_rho_r * bx_r**2, 0.0))))
		
		! David's estimate (not PP or ES)
		sr = max(abs(vx_l) + cfx_l, abs(vx_r) + cfx_r)
		sl = -sr
		
		nflux = ec_flux_x(rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l, beta_l, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r, beta_r, gamma) - 0.5 * sr * (w_r - w_l)
		nBn = 0.5 * (bx_l + bx_r)
	end subroutine EC_LF_flux_x
	
	!DIR$ ATTRIBUTES FORCEINLINE :: aux
	pure subroutine aux(rho, mx, my, mz, bx, by, bz, e_tot, gamma, cfx, flux, entv, ent_flux, phi)
		implicit none
		real, intent(in) :: rho, mx, my, mz, bx, by, bz, e_tot, gamma
		
		real, intent(out) :: phi, ent_flux, cfx
		real, dimension(8), intent(out) :: entv, flux
		
		real :: vx, vy, vz, ex, ey, ez, inv_rho, e_knt, b_sq, p_fld, p_mag, p_tot, b_dot_u, snc_sq, alfven_sq, sum_sq, inv_theta, sp_ent, ent_fun
		
		! intermediate variables
		inv_rho = 1.0 / rho ! specific volume
		vx = inv_rho*mx ! x velocity
		vy = inv_rho*my ! y velocity
		vz = inv_rho*mz ! z velocity
		e_knt = 0.5*inv_rho*(mx**2 + my**2 + mz**2) ! kinetic energy
		ex = by*vz - bz*vy ! x electric field
		ey = bz*vx - bx*vz ! y electric field
		ez = bx*vy - by*vx ! z electric field
		b_sq = bx**2 + by**2 + bz**2 ! squared norm of magnetic field
		p_mag = 0.5*b_sq ! magnetic pressure
		p_fld = max((gamma-1.0)*(e_tot - (e_knt + p_mag)), PP_eps) ! fluid pressure (robust implementation with PP limiter)
		p_tot = p_fld + p_mag ! total pressure
		b_dot_u = bx*vx + by*vy + bz*vz ! inner product of magnetic field and velocity
		snc_sq = gamma * p_fld * inv_rho ! fluid sonic speed squared
		alfven_sq = b_sq * inv_rho ! Alfven speed squared
		sum_sq = snc_sq + alfven_sq
		cfx = sqrt(0.5*(sum_sq + sqrt(max(sum_sq**2 - 4*snc_sq*inv_rho * bx**2, 0.0))))
		
		inv_theta = rho/p_fld ! inverse of temperature
		phi = inv_theta*b_dot_u
		
		! compute flux
		flux(1) = mx ! density flux
		flux(2) = vx*mx - bx*bx + p_tot ! mx flux
		flux(3) = vx*my - by*bx ! my flux
		flux(4) = vx*mz - bz*bx ! mz flux
		flux(5) = 0.0 ! bx flux
		flux(6) = -ez ! by flux
		flux(7) = ey ! bz flux
		flux(8) = vx*(e_tot + p_tot) - bx*b_dot_u ! e_tot flux
		
		sp_ent = log(p_fld) - gamma*log(rho) ! specific entropy (costly) (robustness enhancement)
		ent_fun = -rho * sp_ent / (gamma-1.0) ! entropy function
		ent_flux = ent_fun * vx ! entropy flux function
		
		! entropy variable
		entv(1) = (gamma - sp_ent) / (gamma - 1.0) - e_knt / p_fld
		entv(2) = mx/p_fld
		entv(3) = my/p_fld
		entv(4) = mz/p_fld
		entv(5) = inv_theta*bx
		entv(6) = inv_theta*by
		entv(7) = inv_theta*bz
		entv(8) = -inv_theta
	end subroutine aux
	
	pure function swap_xy(uin) result(uout)
		implicit none
		real, dimension(8), intent(in) :: uin
		real, dimension(8) :: uout
		uout(1) = uin(1)
		uout(2) = uin(3) ! this
		uout(3) = uin(2) ! this
		uout(4) = uin(4)
		uout(5) = uin(6) ! this
		uout(6) = uin(5) ! this
		uout(7) = uin(7)
		uout(8) = uin(8)
	end function swap_xy
	
	!DIR$ ATTRIBUTES FORCEINLINE :: ec_flux_x
	pure function ec_flux_x(rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l, beta_l, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r, beta_r, gamma) result(flux)
		!! very costly procedure
		!! two-point symmetric entropy conservative flux in x-direction
		!! Ref: Entropy Stable Finite Volume Scheme for Ideal Compressible MHD on 2-D Cartesian Meshes
		implicit none
		real, intent(in) :: rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l, beta_l, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r, beta_r, gamma
		real, dimension(8) :: flux
		
		real :: rho_hat, beta_hat
		real :: v_sq_bar, b_sq_bar, beta_bar, rho_bar, vx_bar, vy_bar, vz_bar, bx_bar, by_bar, bz_bar, betavx_bar, betavy_bar, betavz_bar
		
		rho_bar = 0.5 * (rho_l + rho_r)
		vx_bar = 0.5 * (vx_l + vx_r)
		vy_bar = 0.5 * (vy_l + vy_r)
		vz_bar = 0.5 * (vz_l + vz_r)
		bx_bar = 0.5 * (bx_l + bx_r)
		by_bar = 0.5 * (by_l + by_r)
		bz_bar = 0.5 * (bz_l + bz_r)
		beta_bar = 0.5 * (beta_l + beta_r)
		
		v_sq_bar = 0.5 * ((vx_l**2 + vy_l**2 + vz_l**2) + (vx_r**2 + vy_r**2 + vz_r**2))
		b_sq_bar = 0.5 * ((bx_l**2 + by_l**2 + bz_l**2) + (bx_r**2 + by_r**2 + bz_r**2))
		
		betavx_bar = 0.5 * (beta_l*vx_l + beta_r*vx_r)
		betavy_bar = 0.5 * (beta_l*vy_l + beta_r*vy_r)
		betavz_bar = 0.5 * (beta_l*vz_l + beta_r*vz_r)
		
		rho_hat = rho_bar / log_avg_divisor(rho_l, rho_r)
		beta_hat = beta_bar / log_avg_divisor(beta_l, beta_r)
		
		flux(1) = rho_hat * vx_bar
		flux(2) = vx_bar * flux(1) - bx_bar*bx_bar + 0.5 * (b_sq_bar + rho_bar / beta_bar)
		flux(3) = vy_bar * flux(1) - bx_bar*by_bar + 0.0
		flux(4) = vz_bar * flux(1) - bx_bar*bz_bar + 0.0
		flux(5) = 0.0
		flux(6) = (betavx_bar * by_bar - betavy_bar * bx_bar) / beta_bar
		flux(7) = (betavx_bar * bz_bar - betavz_bar * bx_bar) / beta_bar
		flux(8) = 0.5 * (1.0/((gamma-1.0)*beta_hat) - v_sq_bar) * flux(1) + (vx_bar * flux(2) + vy_bar * flux(3) + vz_bar * flux(4)) + (bx_bar * flux(5) + by_bar * flux(6) + bz_bar * flux(7)) - 0.5 * vx_bar * b_sq_bar + (vx_bar * bx_bar + vy_bar * by_bar + vz_bar * bz_bar) * bx_bar
	end function ec_flux_x
	
	!DIR$ ATTRIBUTES FORCEINLINE :: ec_flux_y
	pure function ec_flux_y(rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l, beta_l, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r, beta_r, gamma) result(flux)
		!! very costly procedure
		!! two-point symmetric entropy conservative flux in y-direction
		!! Ref: Entropy Stable Finite Volume Scheme for Ideal Compressible MHD on 2-D Cartesian Meshes
		implicit none
		real, intent(in) :: rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l, beta_l, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r, beta_r, gamma
		real, dimension(8) :: flux
		
		real :: rho_hat, beta_hat
		real :: v_sq_bar, b_sq_bar, beta_bar, rho_bar, vx_bar, vy_bar, vz_bar, bx_bar, by_bar, bz_bar, betavx_bar, betavy_bar, betavz_bar
		
		rho_bar = 0.5 * (rho_l + rho_r)
		vx_bar = 0.5 * (vx_l + vx_r)
		vy_bar = 0.5 * (vy_l + vy_r)
		vz_bar = 0.5 * (vz_l + vz_r)
		bx_bar = 0.5 * (bx_l + bx_r)
		by_bar = 0.5 * (by_l + by_r)
		bz_bar = 0.5 * (bz_l + bz_r)
		beta_bar = 0.5 * (beta_l + beta_r)
		
		v_sq_bar = 0.5 * ((vx_l**2 + vy_l**2 + vz_l**2) + (vx_r**2 + vy_r**2 + vz_r**2))
		b_sq_bar = 0.5 * ((bx_l**2 + by_l**2 + bz_l**2) + (bx_r**2 + by_r**2 + bz_r**2))
		
		betavx_bar = 0.5 * (beta_l*vx_l + beta_r*vx_r)
		betavy_bar = 0.5 * (beta_l*vy_l + beta_r*vy_r)
		betavz_bar = 0.5 * (beta_l*vz_l + beta_r*vz_r)
		
		rho_hat = rho_bar / log_avg_divisor(rho_l, rho_r)
		beta_hat = beta_bar / log_avg_divisor(beta_l, beta_r)
		
		flux(1) = rho_hat * vy_bar
		flux(2) = vx_bar * flux(1) - by_bar*bx_bar + 0.0
		flux(3) = vy_bar * flux(1) - by_bar*by_bar + 0.5 * (b_sq_bar + rho_bar / beta_bar)
		flux(4) = vz_bar * flux(1) - by_bar*bz_bar + 0.0
		flux(5) = (betavy_bar * bx_bar - betavx_bar * by_bar) / beta_bar
		flux(6) = 0.0
		flux(7) = (betavy_bar * bz_bar - betavz_bar * by_bar) / beta_bar
		flux(8) = 0.5 * (1.0/((gamma-1.0)*beta_hat) - v_sq_bar) * flux(1) + (vx_bar * flux(2) + vy_bar * flux(3) + vz_bar * flux(4)) + (bx_bar * flux(5) + by_bar * flux(6) + bz_bar * flux(7)) - 0.5 * vy_bar * b_sq_bar + (vx_bar * bx_bar + vy_bar * by_bar + vz_bar * bz_bar) * by_bar
	end function ec_flux_y
	
	!DIR$ ATTRIBUTES FORCEINLINE :: log_avg_divisor
	pure elemental function log_avg_divisor(LL, RR) result(ff)
		!! very costly procedure
		!! return the divisior that makes the arithmatic average become the logarithmic average
		!! i.e., 0.5*(LL+RR)/ff will become the logarithmic average (LL-RR)/(log(LL)-log(RR))
		!! Two inputs must be of the same sign (non-zero)!!! This affects the stability of the code!!!
		!! The output is symmetric w.r.t. LL and RR.
		!! Ref: Affordable, entropy-consistent Euler flux functions II: Entropy production at shocks
		implicit none
		real, intent(in) :: LL, RR
		real, parameter :: eps = 1.4e-4 ! for single or double precision / float 32 or float64 only
		real, parameter :: a0 = 1.0, a1 = 1.0/3.0, a2 = 1.0/5.0, a3 = 1.0/7.0
		real :: zeta, f, u, ff
		
		zeta = LL / RR
		f = (zeta - 1.0) / (zeta + 1.0) ! f in (-1,1), zeta == (1+f)/(1-f)
		u = f**2 ! u in [0,1) if inputs are of the same sign and non-zero
		
		!! ff = 0.5*log(zeta)/f = arctanh(f)/(f) = \sum_{k=0}^{\infty} f^(2*k)/(2*k+1) = \sum_{k=0}^{\infty} u^k/(2*k+1)
		!! ff >= 1.0
		!! For truncated series, the error can be bounded by (for u in [0, eps)): 
		!! 0 <= ff - \sum_{k=0}^N f^(2*k)/(2*k+1) = \sum_{k=N+1}^{\infty} u^k/(2*k+1) <= u^(N+1)/((1-u)*(2*N+3))
		!! Hence, for single or double precision, eps = 1.4e-4 and N = 3 is enough for error <= 0.25*eps
		
		if (u < eps) then
			ff = a0 + u*(a1 + u*(a2 + u*a3)) ! Horner's summation
		else
			ff = 0.5 * log(zeta) / f ! warning: log is vulnerable to floating point errors
		end if
	end function log_avg_divisor
	
end module DG_solver_m