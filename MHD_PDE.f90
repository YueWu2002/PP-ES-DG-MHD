module MHD_PDE_m
	implicit none
	real, protected, save :: gamma = 5.0/3.0 ! ratio of specific heat for ideal polytopic gas
contains
	subroutine set_gamma(gamma_)
		implicit none
		real, intent(in) :: gamma_
		if (gamma_ /= gamma) print '(A)', 'MHD_PDE: Warning: gamma is changed!'
		gamma = gamma_
	end subroutine set_gamma

	pure function Euler_advective_flux(Uc, p, normal, dim) result(F)
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: Uc ! assume contiguous storage
		real, intent(in) :: p
		real, dimension(:), contiguous, intent(in) :: normal
		real, dimension(dim+2) :: F
		real :: vn
		vn = dot_product(Uc(2:dim+1), normal(1:dim)) / Uc(1)
		F(1) = vn * Uc(1)
		F(2:dim+1) = vn * Uc(2:dim+1) + p * normal(1:dim)
		F(dim+2) = vn * (Uc(dim+2) + p)
	end function Euler_advective_flux
	
	pure function Euler_pressure(Uc, dim) result(p)
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: Uc ! assume contiguous storage
		real :: p
		p = (gamma-1.0) * (Uc(dim+2) - 0.5 * sum(Uc(2:dim+1)**2) / Uc(1))
	end function Euler_pressure

	pure function Euler_con2pri(Uc, dim) result(Up)
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: Uc ! assume contiguous storage
		real, dimension(dim+2) :: Up
		Up(1) = Uc(1)
		Up(2:dim+1) = Uc(2:dim+1) / Uc(1)
		Up(dim+2) = Euler_pressure(Uc, dim)
	end function Euler_con2pri
	
	pure function MHD_con2pri(uc) result(up)
		implicit none
		real, dimension(8), intent(in) :: uc
		real, dimension(8) :: up
		up(1) = uc(1)
		up(2) = uc(2)/uc(1)
		up(3) = uc(3)/uc(1)
		up(4) = uc(4)/uc(1)
		up(5) = uc(5)
		up(6) = uc(6)
		up(7) = uc(7)
		up(8) = (gamma-1.0)*(uc(8) - 0.5*((uc(2)**2 + uc(3)**2 + uc(4)**2)/uc(1) + (uc(5)**2 + uc(6)**2 + uc(7)**2)))
	end function MHD_con2pri
	
	pure function MHD_pri2con(up) result(uc)
		implicit none
		real, dimension(8), intent(in) :: up
		real, dimension(8) :: uc
		uc(1) = up(1)
		uc(2) = up(1)*up(2)
		uc(3) = up(1)*up(3)
		uc(4) = up(1)*up(4)
		uc(5) = up(5)
		uc(6) = up(6)
		uc(7) = up(7)
		uc(8) = up(8)/(gamma-1.0) + 0.5*(up(1)*(up(2)**2 + up(3)**2 + up(4)**2) + (up(5)**2 + up(6)**2 + up(7)**2))
	end function MHD_pri2con
	
	pure function Euler_pri2con(Up, dim) result(Uc)
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: Up ! assume contiguous storage
		real, dimension(dim+2) :: Uc
		Uc(1) = Up(1)
		Uc(2:dim+1) = Up(1) * Up(2:dim+1)
		Uc(dim+2) = Up(dim+2) / (gamma-1.0) + 0.5 * Up(1) * sum(Up(2:dim+1)**2)
	end function Euler_pri2con

	pure subroutine Euler_Roe_average(U_l, U_r, normal, dim, U_roe, c_roe)
		! compute the (modified) Roe's average of primitive variables and the sonic speed on the normal direction
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: U_l, U_r ! assume contiguous storage
		real, dimension(:), contiguous, intent(in) :: normal ! assume contiguous storage
		real, dimension(dim+2), intent(out) :: U_roe ! assume contiguous storage
		real, intent(out), optional :: c_roe
		real :: rho_l, rho_r, v_l(dim), v_r(dim), p_l, p_r, c_l_sq, c_r_sq, R_roe, deno, rho_roe, v_roe(dim), c_roe_sq, p_roe

		rho_l = U_l(1)
		v_l(:) = U_l(2:dim+1)
		p_l = U_l(dim+2)
		c_l_sq = gamma * abs(p_l / rho_l) ! stable implementation

		rho_r = U_r(1)
		v_r(:) = U_r(2:dim+1)
		p_r = U_r(dim+2)
		c_r_sq = gamma * abs(p_r / rho_r) ! stable implementation

		R_roe = sqrt( abs(rho_r / rho_l) )
		deno = 1.0 + R_roe
		rho_roe = sqrt( abs(rho_l * rho_r) )
		v_roe(:) = (v_l(:) + R_roe * v_r(:)) / deno
		c_roe_sq = (0.5*(gamma-1.0) * R_roe / deno**2) * dot_product(v_l(:) - v_r(:), normal(:))**2 + (c_l_sq + R_roe * c_r_sq) / deno
		p_roe = rho_roe * c_roe_sq / gamma
		U_roe = [rho_roe, v_roe(:), p_roe]
		if(present(c_roe)) c_roe = sqrt(c_roe_sq)
	end subroutine Euler_Roe_average
	
	pure subroutine Euler_Roe_average_1D(U_l, U_r, U_roe, c_roe)
		! compute the (modified) Roe's average of primitive variables and the sonic speed on the normal direction
		implicit none
		real, dimension(3), intent(in) :: U_l, U_r ! assume contiguous storage
		real, dimension(3), intent(out) :: U_roe ! assume contiguous storage
		real, intent(out) :: c_roe
		real :: R_roe, deno, c_roe_sq, wl, wr

		R_roe = sqrt( abs( U_r(1) / U_l(1) ) )
		deno = 1.0 + R_roe
		wl = 1.0 / deno
		wr = R_roe / deno
		
		U_roe(1) = abs(R_roe * U_l(1))
		U_roe(2) = wl * U_l(2) + wr * U_r(2)
		c_roe_sq = (0.5*(gamma-1.0) * R_roe / deno**2) * (U_l(2) - U_r(2))**2 + (wl * abs(gamma * U_l(3) / U_l(1)) + wr * abs(gamma * U_r(3) / U_r(1)))
		U_roe(3) = U_roe(1) * c_roe_sq / gamma
		c_roe = sqrt(c_roe_sq)
	end subroutine Euler_Roe_average_1D

	pure subroutine Euler_eigenmat_con(U_p, normal, dim, Lmat, Rmat, eigs)
		! compute the left and right eigen matrices in the normal direction for compressible Euler equations (for ideal gas only)
		! We only support dim==1 and dim==2. Otherwise unexpected errors might occur!
		implicit none
		integer, intent(in) :: dim
		real, dimension(dim+2), intent(in) :: U_p ! assume contiguous storage
		real, dimension(:), contiguous, intent(in) :: normal ! assume contiguous storage
		real, dimension(dim+2, dim+2), intent(out) :: Lmat, Rmat ! assume contiguous storage
		real, dimension(dim+2), intent(out), optional :: eigs ! assume contiguous storage
		real :: rho, v(dim), vn, ek, p, T, c_sq, c, H
		integer :: i

		rho = U_p(1)
		v(:) = U_p(2:dim+1)
		vn = dot_product(v(:), normal(:))
		ek = 0.5 * sum(v(:)**2)
		p = U_p(dim+2)
		T = abs(p/rho)
		c_sq = gamma * T
		c = sqrt(c_sq)
		H = (gamma/(gamma-1.0)) * T + ek

		if(present(eigs)) then
			eigs(:) = vn
			eigs(1) = vn - c
			eigs(3) = vn + c
		end if

		Rmat(:,1) = (1.0/c) * [1.0, v(:), H] - [0.0, normal(:), vn]
		Rmat(:,2) = rho * [1.0, v(:), ek]
		Rmat(:,3) = (1.0/c) * [1.0, v(:), H] + [0.0, normal(:), vn]

		Lmat(1,:) = ((gamma-1.0)/(2*c)) * [ek, -v(:), 1.0] + 0.5 * [vn, -normal(:), 0.0] ! left sonic wave
		Lmat(2,:) = (1.0/rho) * [1.0, [(0.0, i=1,dim)], 0.0] - ((gamma-1.0)/(rho*c_sq)) * [ek, -v(:), 1.0] ! entropy contact wave
		Lmat(3,:) = ((gamma-1.0)/(2*c)) * [ek, -v(:), 1.0] - 0.5 * [vn, -normal(:), 0.0] ! right sonic wave

		select case (dim)
		case(1)
			! do nothing
		case(2)
			! the orthogonal unit vector to normal is [-normal(2), normal(1)]
			Rmat(:,4) = [0.0, [-normal(2), normal(1)], dot_product([-normal(2), normal(1)], v(:))]
			Lmat(4,:) = [-Rmat(4,4), [-normal(2), normal(1)], 0.0]
		case default
			error stop "Euler_eigenmat_con not implemented for dim>2"
		end select
	end subroutine Euler_eigenmat_con
end module MHD_PDE_m