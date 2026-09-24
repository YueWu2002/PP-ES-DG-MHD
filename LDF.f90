module LDF_m
	! module for locally divergence-free mapping
	! Derivation of LDF projection:
	! We are using W = Q_{px,py}^2 space (dim = 2*(px+1)*(py+1)), whose divergence space is M = Q_{px-1,py} + Q_{px,py-1} (dim = (px+1)*(py+1)-1). 
	! Inspired by practice in mixed finite elements and the geom_free flag in NGSolve (and also the article "Arbitrarily high-order globally divergence-free DG method for compressible ideal MHD equations on unstructured meshes"), consider the following transformations that maps functions on K to those on reference element Khat, based upon the physical coordinate mapping T: Khat -> K: x = T(xhat). 
	! 1. div-conforming contravariant Piola mapping for vector-valued functions u: uhat(xhat) = det(J) * J^{-1} * u(x), 
	! 2. grad-conforming Piola mapping for scalar function p: phat(xhat) = p(x).
	! The first one satisfies divhat(uhat) = det(J) * div(u), so u is LDF i.f.f. uhat is LDF. 
	! The underlying projection scheme on the reference element is: 
	!    find (uhat, phat) in What*Mhat s.t. for any (vhat, qhat) in What*Mhat: 
	!    (uhat, vhat)_Khat + (phat, divhat(vhat))_Khat = (u0hat, vhat)_Khat,
	!    (divhat(uhat), qhat)_Khat = 0. 
	! For LDF property of uhat, we need Mhat \supset divhat(What). For inf-sup stability / well-posedness, we need Mhat \subset divhat(What). Hence we must take Mhat = divhat(What). The matrix equation is 
	!    [diag(M,M), B^T*M] [uhat] = [diag(M,M)*u0hat],
	!    [M*B,       0] [phat] = [0],
	! where B is the divergence mapping and M is the mass matrix for scalar field. The exact solution for uhat is:
	!    uhat = (I - diag(M,M)^{-1} * B^T * sqrt(M) * S^{-1} * sqrt(M) * B) * u0hat, 
	! where the Schur complement matrix (Do not eliminate sqrt(M) because it helps to control conditonal number! Numerical experience shows that sqrt(M) is better than either 1 or M.)
	!    S = sqrt(M) * B * diag(M,M)^{-1} * B^T * sqrt(M)
	! is symmetric positive definite. The associated continuous problem is:
	!    find (u, p) in W*M s.t. for any (v, q) in W*M: 
	!    (J^{-1}*u; det(J); J^{-1}*v)_K + (p, div(v))_K = (J^{-1}*u0; det(J); J^{-1}*v)_K, 
	!    (div(u), q)_K = 0. 
	! In regard of u, the above formulation finds the best LDF projector of u0 in the (J^{-1}*u; det(J); J^{-1}*v)_K inner product, which agree with "Arbitrarily high-order globally divergence-free DG method for compressible ideal MHD equations on unstructured meshes" in the constant Jacobian case. 
	! In computation, we first use nodal values of u0 to compute the Legendre modal coefficient of u0hat, then proceed as above for the Legendre modal coefficient of uhat, and finally obtain the nodal values of the LDF projected u. Using Legendre modal basis will make the Schur complement sparse.
	! With the help of Legendre basis, assembly of S is easy as follows. Suppose we store all Legendre modal coefs of Bx then By, the index vary first in x then in y. Then the matrices (with hierarchical Legendre basis; on reference element) are: (All are 1D matrices.)
	!    M = M_py otimes M_px, 
	!    B = B_full(1:end-1, :), where:
	!    B_full = [I_py otimes D_px, D_py otimes I_px] (Last row are 0's.),
	!    S = S_full(1:end-1, 1:end-1), where:
	!    S_full = sqrt(M) * B_full * diag(M,M)^{-1} * B_full^T * sqrt(M) (Last row and column are 0's.)
	!          = I_py otimes (sqrt(M_px)*D_px*M_px^{-1}*D_px^T*sqrt(M_px)) + (sqrt(M_py)*D_py*M_py^{-1}*D_py^T*sqrt(M_py)) otimes I_px
	!          = I_py otimes S_full_px + S_full_py otimes I_px 
	!           (S_full_px and S_full_py are 1D full Schur cores: S_full_px = sqrt(M_px)*D_px*M_px^{-1}*D_px^T*sqrt(M_px)). 
	! Here, S_full is symmetric positive semi-definite with only one zero eigenvalue for the last row and column, while S is s.p.d.
	! To simplify S_full and invert S, we first decompose:
	!    S_full_px = G_px * G_px^T, where G_px = sqrt(M_px)*D_px*sqrt(M_px)^{-1},
	!    S_full_py = G_py * G_py^T, where G_py = sqrt(M_py)*D_py*sqrt(M_py)^{-1}.
	! Then conduct SVD decomposition (use LAPACK routine DGESVD or more advanced DGESDD): 
	!    G_px = U_px * L_px * V_px^T, (U_px and V_px has orthonormal columns, L_px is diagonal and descending order)
	!    G_py = U_py * L_py * V_py^T, (U_py and V_py has orthonormal columns, L_py is diagonal and descending order)
	! which implies
	!    S_full_px = U_px * L_px.^2 * U_px^T,
	!    S_full_py = U_py * L_py.^2 * U_py^T.
	! Then there holds:
	!    (U_py^T otimes U_px^T) * S_full * (U_py otimes U_px) = I_py otimes L_px.^2 + L_py.^2 otimes I_px := L_full, 
	! where L_full's last element is zero,
	! or equivalently,
	!    S_full = (U_py^{-T} otimes U_px^{-T}) * L_full * (U_py^{-1} otimes U_px^{-1}),
	! thus:
	!    S^{-1} = (U_py otimes U_px) * L_full^{†} * (U_py^T otimes VU_px^T),
	! so solution is easy:
	!    uhat = u0hat - ucorr, 
	! where:
	!    ucorr = diag(M,M)^{-1} * B^T * sqrt(M) * S^{-1} * sqrt(M) * B * u0hat 
	!          = diag(M,M)^{-1} * B_full^T * sqrt(M) * (U_py otimes V_px) * L_full^{†} * (U_py^T otimes U_px^T) * sqrt(M) * B_full * u0hat
	!          = diag(M,M)^{-1} * R^T * L_full^{†} * R * u0hat,
	! where:
	!    R = (U_py^T otimes U_px^T) * sqrt(M) * B_full
	! Application of R, R^T and diag(M,M)^{-1} can be accelerated using the tensor structure. Application of L_full^{†} can be achieved through entry-wise operation due to diagonality.
	! 
! matlab test code (with my LGLDGSEM code)
! % set nx and ny
! [Mx,~,Dx,~] = Legendre_basis_mat(nx);
! [My,~,Dy,~] = Legendre_basis_mat(ny);
! Idx = speye(nx+1);
! Idy = speye(ny+1);
! B_full = [kron(Idy,Dx), kron(Dy,Idx)];
! B = B_full(1:end-1,:);
! M = kron(My, Mx);
! MM = blkdiag(M, M);
! S_full = sqrt(M)*B_full*inv(MM)*B_full'*sqrt(M); % full Schur
! S = S_full(1:end-1, 1:end-1); % core Schur
! A = MM \ (B' * sqrt(M(1:end-1,1:end-1)) * inv(S) * sqrt(M(1:end-1,1:end-1)) * B); % correction mapping
! Gx = sqrt(Mx)*Dx*inv(sqrt(Mx)); [Ux,Lx,Vx]=svd(full(Gx),"vector");
! Gy = sqrt(My)*Dy*inv(sqrt(My)); [Uy,Ly,Vy]=svd(full(Gy),"vector");
! Llist = kron(ones(ny+1,1), Lx.^2) + kron(Ly.^2, ones(nx+1,1));
! Ldaglist = 1.0 ./ Llist; Ldaglist(end) = 0.0;
! Ldag = sparse(1:(nx+1)*(ny+1), 1:(nx+1)*(ny+1), Ldaglist);
! R = kron(Uy'*sqrt(My), Ux'*sqrt(Mx)) * B_full;
! A1 = MM \ (R' * Ldag * R);
! norm(full(A)-full(A1),1)
! norm(sqrt(M)*B_full*(eye(size(A)) - A)) % check that the result is LDF (L2 norm of div)
! norm(sqrt(M)*B_full*(eye(size(A1)) - A1)) % check that the result is LDF (L2 norm of div)
! P1 = eye(size(A1)) - A1; % projection matrix
! norm(P1*P1 - P1, 1) % check that P is a projection matrix / idempotent operator
! norm((MM*P1)'-MM*P1, 1) % check that P is self-adjoint under the inner product with Gram matrix MM
! % At this point, we found P1 is better than P in that its eigenvalues do not have small imaginary error part.
	implicit none
	
	! locally divergence-free projector for Q_{p1,p2}^2 magnetic field
	! only for on the reference element [-1,1]^2 with Legendre modal coefficients as input
	! For physical fields: first use contravariant Piola mapping to transform onto Legendre coefs on [-1,1]^2
	type, public :: LDF_projector
		private ! default private
		logical :: ready = .false.
		integer :: px = -1, py = -1
		real, allocatable :: Lx(:), Ly(:)						! eigenvalues of Sx and Sy
		real, allocatable :: mx(:), my(:), invmx(:), invmy(:)	! diagonal Legendre mass matrix and inverse
		real, allocatable :: sqrtm2d(:,:), invm2d(:,:)			! square root and inverse of 2D mass matrix
		real, allocatable :: Dx(:,:), Dy(:,:)					! Legendre modal differentiation matrices
		real, allocatable :: Ux(:,:), Uy(:,:)					! orthonormal basis for SVD decompositions
		real, allocatable :: Ldagger(:,:)						! diagonal Moore-Penrose pseudo-inverse of L
	contains
		procedure, pass :: init
		procedure, pass :: apply ! apply the LDF projection (on reference element)
		procedure, pass :: div ! compute the divergence (on reference element)
		final :: destroy
	end type LDF_projector
	
	private ! default private
contains

	subroutine init(this, px, py)
#ifdef _OPENMP
		use omp_lib
#endif
		use, intrinsic :: iso_fortran_env, only: stdin  => input_unit, stdout => output_unit, stderr => error_unit ! need Fortran 2003
		use gauss_m, only: gauss_legendre_matrix
		implicit none
		class(LDF_projector), intent(inout) :: this
		integer, intent(in) :: px, py
        ! use interface rather than KWD "external" to enhance stability
        ! The following interfaces are adapted from the f77_lapack_single_double.f90 of the Netlib-LAPACK95 source code.
        interface
        PURE SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
            implicit none
            CHARACTER(LEN=1), INTENT(IN) :: JOBZ
            INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
            INTEGER, INTENT(OUT) :: INFO
            REAL(8), INTENT(OUT) :: S(*)
            REAL(8), INTENT(INOUT) :: A(LDA,*)
            REAL(8), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
            INTEGER, INTENT(OUT) :: IWORK(*)
            END SUBROUTINE DGESDD
        end interface
		integer :: i, j
		if (px<0 .or. py<0) error stop 'LDF_projector: error para in init'
		this%px = px
		this%py = py
		allocate(this%Lx(px+1), this%Ly(py+1), this%mx(px+1), this%my(py+1), this%invmx(px+1), this%invmy(py+1))
		allocate(this%sqrtm2d(px+1,py+1), this%invm2d(px+1,py+1), this%Dx(px+1,px+1), this%Dy(py+1,py+1), this%Ux(px+1,px+1), this%Uy(py+1,py+1), this%Ldagger(px+1,py+1))
		
		block ! obtain the Legendre modal matrices
		real :: tlx(0:px), trx(0:px), tly(0:py), try(0:py)
		call gauss_legendre_matrix(px, this%mx, this%invmx, this%Dx, tlx, trx)
		call gauss_legendre_matrix(py, this%my, this%invmy, this%Dy, tly, try)
		end block
		
		do j = 1, py+1
			do i = 1, px+1
				this%sqrtm2d(i,j) = sqrt(this%mx(i) * this%my(j))
			end do
		end do
		do j = 1, py+1
			do i = 1, px+1
				this%invm2d(i,j) = this%invmx(i) * this%invmy(j)
			end do
		end do
		
		block ! compute the SVD decompositions
		integer :: ierrx, ierry, lwork
		real(8) :: sqrtmx(px+1), sqrtmy(py+1), Gx(px+1,px+1), Gy(py+1,py+1), Lx(px+1), Ly(py+1), Ux(px+1,px+1), Uy(py+1,py+1), VxT(px+1,px+1), VyT(py+1,py+1)
		real(8), allocatable :: work(:)
		integer :: iwork(8*max(px+1,py+1))
		sqrtmx = sqrt(this%mx)
		sqrtmy = sqrt(this%my)
		do j = 1, px+1
			do i = 1, px+1
				Gx(i,j) = sqrtmx(i) * this%Dx(i,j) / sqrtmx(j)
			end do
		end do
		do j = 1, py+1
			do i = 1, py+1
				Gy(i,j) = sqrtmy(i) * this%Dy(i,j) / sqrtmy(j)
			end do
		end do
		
		! use the original Fortran77-style API
		allocate(work(1)) ! pre-allocate
		!<---begin for Gx--->
		call DGESDD('A', px+1, px+1, Gx, px+1, Lx, Ux, px+1, VxT, px+1, work, -1, iwork, ierrx) ! workspace query by work(1)
		lwork = ceiling(work(1)) ! optimal workspace size
		deallocate(work)
		allocate(work(lwork)) ! re-allocate
		call DGESDD('A', px+1, px+1, Gx, px+1, Lx, Ux, px+1, VxT, px+1, work, lwork, iwork, ierrx) ! will destroy Gx
		!<---end for Gx--->
		!<---begin for Gy--->
		call DGESDD('A', py+1, py+1, Gy, py+1, Ly, Uy, py+1, VyT, py+1, work, -1, iwork, ierry) ! workspace query by work(1)
		lwork = ceiling(work(1)) ! optimal workspace size
		deallocate(work)
		allocate(work(lwork)) ! re-allocate
		call DGESDD('A', py+1, py+1, Gy, py+1, Ly, Uy, py+1, VyT, py+1, work, lwork, iwork, ierry) ! will destroy Gy
		!<---end for Gy--->
		deallocate(work)
		if (ierrx/=0 .or. ierry/=0) error stop 'SVD failed!'
        
        ! save the results
        this%Ux(:,:) = Ux(:,:)
        this%Uy(:,:) = Uy(:,:)
        this%Lx(:) = Lx(:)
        this%Ly(:) = Ly(:)
        
		end block
		
		do j = 1, py+1
			do i = 1, px+1
				if (i<=this%px .or. j<=this%py) then
					this%Ldagger(i,j) = 1.0 / (this%Lx(i)**2 + this%Ly(j)**2)
				else
					this%Ldagger(i,j) = 0.0
				end if
			end do
		end do
		
		this%ready = .true.
	end subroutine init
	
	pure subroutine apply(this, uu)
		!! maps the Legendre modal coefficients (2D vector on reference element [-1,1]^2) to those of the locally divergence-free projected one
		implicit none
		class(LDF_projector), intent(in) :: this
		real, dimension(this%px+1, this%py+1, 2), intent(inout) :: uu ! Legendre coefficients on [-1,1]^2
		
		! real, dimension(this%px+1, this%py+1, 2) :: ucorr
		real, dimension(this%px+1, this%py+1) :: udiv
		
		if (.not. this%ready) error stop 'LDF_projector: not ready!'
		
		! ! apply R
		! udiv = this%sqrtm2d * (matmul(this%Dx, uu(:,:,1)) + matmul(uu(:,:,2), transpose(this%Dy))) ! apply sqrt(M)*div
		! udiv = matmul(matmul(transpose(this%Ux), udiv), this%Uy) ! inner product with U
		
		! ! apply pseudo inverse of L
		! udiv = this%Ldagger * udiv
		
		! ! apply bigMinv * R^T
		! udiv = this%sqrtm2d * (matmul(matmul(this%Ux, udiv), transpose(this%Uy))) ! apply sqrt(M)*U
		
		! ucorr(:,:,1) = this%invm2d * matmul(transpose(this%Dx), udiv) ! apply bigMinv*div^T
		! ucorr(:,:,2) = this%invm2d * matmul(udiv, this%Dy) ! apply bigMinv*div^T
		
		! ! perform correction
		! uu = uu - ucorr
		
		! alternative: combined version, faster
		udiv = this%sqrtm2d * (matmul(matmul(this%Ux, this%Ldagger * matmul(matmul(transpose(this%Ux), this%sqrtm2d * (matmul(this%Dx, uu(:,:,1)) + matmul(uu(:,:,2), transpose(this%Dy)))), this%Uy)), transpose(this%Uy)))
		uu(:,:,1) = uu(:,:,1) - this%invm2d * matmul(transpose(this%Dx), udiv)
		uu(:,:,2) = uu(:,:,2) - this%invm2d * matmul(udiv, this%Dy)
	end subroutine apply
	
	pure subroutine div(this, uu, udiv)
		!! maps the Legendre modal coefficients (2D vector on reference element [-1,1]^2) to those of the divergence
		implicit none
		class(LDF_projector), intent(in) :: this
		real, dimension(this%px+1, this%py+1, 2), intent(in) :: uu ! Legendre coefficients on [-1,1]^2
		real, dimension(this%px+1, this%py+1), intent(out) :: udiv ! Legendre coefficients on [-1,1]^2
		if (.not. this%ready) error stop 'LDF_projector: not ready!'
		udiv = matmul(this%Dx, uu(:,:,1)) + matmul(uu(:,:,2), transpose(this%Dy))
	end subroutine div
	
	subroutine destroy(this)
		!! default destructor (automatically runned; cannot be explicitly called)
		use, intrinsic :: iso_fortran_env, only: stdin  => input_unit, stdout => output_unit, stderr => error_unit ! need Fortran 2003
		implicit none
		type(LDF_projector), intent(inout) :: this
		this%ready = .false.
		deallocate(this%mx, this%my, this%invmx, this%invmy, this%Lx, this%Ly, this%sqrtm2d, this%invm2d, this%Dx, this%Dy, this%Ux, this%Uy, this%Ldagger)
	end subroutine destroy
end module LDF_m