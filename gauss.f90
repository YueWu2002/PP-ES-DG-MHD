module gauss_m
	! modules for orthogonal polynomials, quadratures and nodal spectral element methods
	! Refs: 
	! [1] Shen, J., Tang, T., Wang, LL. (2011). Orthogonal Polynomials and Related Approximation Results. In: Spectral Methods. Springer Series in Computational Mathematics, vol 41. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-71041-7_3.
	! [2] https://mathworld.wolfram.com/LobattoQuadrature.html
	! [3] Teukolsky, S. A. (2015). Short note on the mass matrix for Gauss–Lobatto grid points. Journal of Computational Physics, 283, 408-413.
	implicit none
	real, private :: &
		sqrt2 = 1.4142135623730950488016887242096980785696718753769, &
		sqrt5 = 2.2360679774997896964091736687312762354406183596115, &
		sqrt3 = 1.7320508075688772935274463415058723669428052538104, &
		sqrt7 = 2.6457513110645905905016157536392604257102591830825, &
		sqrt11 = 3.3166247903553998491149327366706866839270885455894, &
		sqrt15 = 3.8729833462074168851792653997823996108329217052916, &
		sqrt21 = 4.5825756949558400065880471937280084889844565767680
contains
	pure subroutine gauss_lobatto(n, x, w, interval)
	!!   compute the n-point Gauss-Lobatto quadrature points and weights on [-1,1]
	!!   algebraic exactness: 2n-3, convergence rate when using composite integration: 2n-2
	!! 
	!! input: 
	!!   n:  number of points, n>=1
	!! 
	!! output: (exact for 1<=n<=7)
	!!   x:  quadrature nodes, size = [n, 1]
	!!   w:  quadrature weights, size = [n, 1]
	!!   interval: optional, if provided, map the points from [-1,1] to [a,b]
		implicit none
		integer, intent(in) :: n
		real, intent(out) :: x(n), w(n)
		real, intent(in), optional :: interval(2)

		! extract the reference Legendre-Gauss-Lobatto points and weights on [-1,1]
		select case (n)
		case (1)
			x = [0.0]
			w = [2.0]
		case (2)
			x = [-1.0, 1.0]
			w = [1.0, 1.0]
		case (3)
			x = [-1.0, 0.0, 1.0]
			w = [1.0/3.0, 4.0/3.0, 1.0/3.0]
		case (4)
			x = [-1.0, -1.0/sqrt5, 1.0/sqrt5, 1.0]
			w = [1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0]
		case (5)
			x = [-1.0, -sqrt21/7.0, 0.0, sqrt21/7.0, 1.0]
			w = [0.1, 49.0/90.0, 32.0/45.0, 49.0/90.0, 0.1]
		case (6)
			x = [-1.0, -sqrt((7.0 + 2*sqrt7)/21.0), -sqrt((7.0 - 2*sqrt7)/21.0), sqrt((7.0 - 2*sqrt7)/21.0), sqrt((7.0 + 2*sqrt7)/21.0), 1.0]
			w = [1.0/15.0, (14.0 - sqrt7)/30.0, (14.0 + sqrt7)/30.0, (14.0 + sqrt7)/30.0, (14.0 - sqrt7)/30.0, 1.0/15.0]
		case (7)
			x = [-1.0, -sqrt((15.0 + 2*sqrt15)/33.0), -sqrt((15.0 - 2*sqrt15)/33.0), 0.0, sqrt((15.0 - 2*sqrt15)/33.0), sqrt((15.0 + 2*sqrt15)/33.0), 1.0]
			w = [1.0/21.0, (124.0 - 7*sqrt15)/350.0, (124.0 + 7*sqrt15)/350.0, 256.0/525.0, (124.0 + 7*sqrt15)/350.0, (124.0 - 7*sqrt15)/350.0, 1.0/21.0]
		case default
			error stop "gauss_lobatto: n must be in [1,7], otherwise not implemented."

			! See lobpts.m and jacpts.m in Chebfun for higher n implementation.
			!
			! The points are the roots of (1-x^2)*P'_(n-1)(x)=0, where P'_(n-1) is the first derivative of the (n-1)th Legendre polynomial.
			! The interior points are also the same as the roots of P^(1,1)_{n-2}(x), so roots of Jacobi polynomials can be used to find them.
			! We focus on N<100 here for simplicity. The Newton iteration method is used.
		end select
		
		! map to the interval [a,b] if specified
		if (present(interval)) then
			x = (0.5*(interval(2)-interval(1))) * x + 0.5*(interval(1)+interval(2))
			w = (0.5*(interval(2)-interval(1))) * w
		end if
	end subroutine gauss_lobatto

	pure function gauss_lobatto_test(interval) result(out)
	!! test Gauss Lobatto quadrature with monomials
		implicit none
		real, intent(in), optional :: interval(2)
		real, dimension(2) :: interval_
		real :: out
		integer :: n, p
		real :: exact, approx, error
		real, allocatable :: x(:), w(:)

		if (present(interval)) then
			interval_ = interval
		else
			interval_ = [-1.0, 1.0]
		end if

		error = 0.0
		do n = 1, 7
			allocate(x(n), w(n))
			if (present(interval)) then
				call gauss_lobatto(n, x, w, interval)
			else
				call gauss_lobatto(n, x, w)
			end if
			do p = 0, 2*n-3
				exact = (interval_(2)**(p+1) - interval_(1)**(p+1))/(p+1)
				approx = sum(x**p * w)
				error = error + abs(exact - approx) / (abs(exact)+1.0)
			end do
			deallocate(x, w)
		end do
		out = error
	end function gauss_lobatto_test

	pure function gauss_lobatto_test_legendre() result(out)
		!! test algebraic degree of exactness of Gauss Lobatto quadrature with Legendre polynomials
		implicit none
		real :: out
		integer :: n, p
		real :: exact, approx, error
		real, allocatable :: x(:), w(:), y(:)

		error = 0.0
		do n = 1, 7
			allocate(x(n), w(n), y(n))
			call gauss_lobatto(n, x, w)
			do p = 0, 2*n-3
				if (p==0) then
					exact = 2.0
				else
					exact = 0.0
				end if
				call legendre(p, x, y)
				approx = sum(y * w)
				error = error + abs(exact - approx) / (abs(exact)+1.0)
			end do
			deallocate(x, w, y)
		end do
		out = error
	end function gauss_lobatto_test_legendre

	pure subroutine gauss_lobatto_matrix(p, M, Minv, D, Pmn, Pnm)
		!! extract the vectors and matrices for p-th order Gauss--Lobatto nodal space on the 1D interval [-1,1]
		!! using Lobatto nodal basis and Legendre modal basis
		!! 
		!! input: 
		!!   p:      polynommial order
		!! 
		!! output: 
		!!   M:      mass matrix with respect to nodal DOFs, size = [p+1, p+1]
		!!   Minv:   inverse of the mass matrix, size = [p+1, p+1]
		!!   D:      differentiation matrix (each column is the coefficient vector of a basis function), size = [p+1, p+1]
		!!   Tl:     left trace matrix, size = [1, p+1] (not outputed)
		!!   Tr:     right trace matrix, size = [1, p+1] (not outputed)
		!!   Pmn:    maps Legendre modal coefs to nodal values (through evaluation), size = [p+1, p+1]
		!!   Pnm:    maps nodal values to Legendre modal coefs (through interpolation), size = [p+1, p+1]
		!! 
		!! Remark: 
		!!   SBP property 1: (D'*diag(qw) + diag(qw)*D) - (Tr'*Tr - Tl'*Tl) = 0
		!!   SBP property 2: (D'*M + M*D) - (Tr'*Tr - Tl'*Tl) = 0
		implicit none
		integer, intent(in) :: p
		real, dimension(p+1,p+1), intent(out) :: M, Minv, D, Pmn, Pnm
		integer :: i,j

		! We can reduce nearly half of the element by using symmetric.
		! We only specify:
		!   M and invM by the left triangular part (M(i,j)==M(j,i)==M(p+2-i,p+2-j)==M(p+2-j,p+2-i))
		!   D by the lower triangular part (D(i,j)=-D(p+2-i,p+2-j))
		!   Pmn by its upper half
		!   Pnm by its left half
		select case (p)
		case (0)
			M(1,1) = 2.0
			Minv(1,1) = 0.5
			D(1,1) = 0.0
			Pmn(1,1) = 1.0
			Pnm(1,1) = 1.0

		case (1)
			M(1:2,1) = [2.0/3.0, 1.0/3.0]

			Minv(1:2,1) = [2.0, -1.0]

			D(1:2,1) = [-0.5, -0.5]
			D(2,2) = 0.5

			Pmn(1,:) = [1.0, -1.0]

			Pnm(:,1) = [0.5, -0.5]

		case (2)
			M(1:3,1) = [4.0/15.0, 2.0/15.0, -1.0/15.0]
			M(2,2) = 16.0/15.0

			Minv(1:3,1) = [4.5, -0.75, 1.5]
			Minv(2,2) = 1.125

			D(1:3,1) = [-1.5, -0.5, 0.5]
			D(2:3,2) = [0.0, -2.0]
			D(3,3) = 1.5

			Pmn(1,:) = [1.0, -1.0, 1.0]
			Pmn(2,:) = [1.0, 0.0, -0.5]

			Pnm(:,1) = [1.0/6.0, -0.5, 1.0/3.0]
			Pnm(:,2) = [2.0/3.0, 0.0, -2.0/3.0]

		case (3)
			M(1:4,1) = [1.0/7.0, sqrt5/42.0, -sqrt5/42.0, 1.0/42.0]
			M(2:3,2) = [5.0/7.0, 5.0/42.0]

			Minv(1:4,1) = [8.0, -2.0/sqrt5, 2.0/sqrt5, -2.0]
			Minv(2:3,2) = [1.6, -0.4]

			D(1,1) = -3.0
			D(2,1:2) = [0.25*(-1.0 - sqrt5), 0.0]
			D(3,1:3) = [0.25*(-1.0 + sqrt5), -0.5*sqrt5, 0.0]
			D(4,1:4) = [-0.5, 1.25*(-1.0 + sqrt5), -1.25*(1.0 + sqrt5), 3.0]

			Pmn(1,:) = [1.0, -1.0, 1.0, -1.0]
			Pmn(2,:) = [1.0, -1.0/sqrt5, -0.2, 1.0/sqrt5]

			Pnm(:,1) = [1.0/12.0, -0.25, 5.0/12.0, -0.25]
			Pnm(:,2) = [5.0/12.0, -0.25*sqrt5, -5.0/12.0, 0.25*sqrt5]

		case (4)
			M(1:5,1) = [4.0/45.0, 7.0/270.0, -4.0/135.0, 7.0/270.0, -1.0/90.0]
			M(2:4,2) = [196.0/405.0, 28.0/405.0, -49.0/810.0]
			M(3,3) = 256.0/405.0

			Minv(1:5,1) = [12.5, -15.0/14.0, 0.9375, -15.0/14.0, 2.5]
			Minv(2:4,2) = [225.0/98.0, -45.0/112.0, 45.0/98.0]
			Minv(3,3) = 225.0/128.0

			D(1,1) = -5.0
			D(2,1:2) = [-3.0/28.0*(7.0 + sqrt21), 0.0]
			D(3,1:3) = [0.375, -0.875*sqrt21/3.0, 0.0]
			D(4,1:4) = [3.0/28.0*(-7.0 + sqrt21), 0.5*sqrt21/3.0, -8.0/sqrt21, 0.0]
			D(5,1:5) = [0.5, 7.0/12.0*(-7.0 + sqrt21), 8.0/3.0, -7.0/12.0*(7.0 + sqrt21), 5.0]

			Pmn(1,:) = [1.0, -1.0, 1.0, -1.0, 1.0]
			Pmn(2,:) = [1.0, -sqrt21/7.0, 1.0/7.0, 3.0/49.0*sqrt21, -3.0/7.0]
			Pmn(3,:) = [1.0, 0.0, -0.5, 0.0, 0.375]
			
			Pnm(:,1) = [0.05, -0.15, 0.25, -0.35, 0.2]
			Pnm(:,2) = [49.0/180.0, -0.35*sqrt21/3.0, 7.0/36.0, 0.35*sqrt21/3.0, -7.0/15.0]
			Pnm(:,3) = [16.0/45.0, 0.0, -8.0/9.0, 0.0, 8.0/15.0]

		case (5)
			M(1:6,1) = [2.0/33.0, sqrt(28.0 - 2*sqrt7)/330.0, -sqrt(7.0 + 0.5*sqrt7)/165.0, sqrt(7.0 + 0.5*sqrt7)/165.0, -sqrt(28.0 - 2*sqrt7)/330.0, 1.0/165.0]
			M(2:5,2) = [(14.0 - sqrt7)/33.0, sqrt21/110.0, -sqrt21/110.0, (14.0 - sqrt7)/330.0]
			M(3:4,3) = [(14.0 + sqrt7)/33.0, (14.0 + sqrt7)/330.0]

			Minv(1:6,1) = [18.0, -sqrt((28.0 + 2*sqrt7)/21.0), sqrt((28.0 - 2*sqrt7)/21.0), -3.0/sqrt(7.0 + 0.5*sqrt7), sqrt((28.0 + 2*sqrt7)/21.0), -3.0]
			Minv(2:5,2) = [(56.0 + 4*sqrt7)/21.0, -2.0/sqrt21, 2.0/sqrt21, -(28.0 + 2*sqrt7)/63.0]
			Minv(3:4,3) = [(56.0 - 4*sqrt7)/21.0, (-28.0 + 2*sqrt7)/63.0]

			D(1,1) = -7.5
			D(2,1:2) = [-(2.0 + sqrt7 + sqrt(21.0 + 6*sqrt7))/6.0, 0.0]
			D(3,1:3) = [-(2.0 - sqrt7 - sqrt(21.0 - 6*sqrt7))/6.0, -9.0/sqrt(58.0*(1.0 - sqrt3/sqrt7) - 8*(sqrt3 - sqrt7)), 0.0]
			D(4,1:4) = [-(2.0 - sqrt7 + sqrt(21.0 - 6*sqrt7))/6.0, sqrt(58*(7.0 - sqrt21) + 56*(sqrt3 - sqrt7))/12.0, -0.5*sqrt(7.0 + 2*sqrt7), 0.0]
			D(5,1:5) = [-(2.0 + sqrt7 - sqrt(21.0 + 6*sqrt7))/6.0, -0.5*sqrt(7.0 - 2*sqrt7), sqrt(58*(7.0 - sqrt21) - 56*(sqrt3 - sqrt7))/12.0, -9.0/sqrt(58*(1.0 - sqrt3/sqrt7) + 8*(sqrt3 - sqrt7)), 0.0]
			D(6,1:6) = [-0.5, 0.25*(7.0 + 4*sqrt7)*(-1.0 + sqrt(7.0 - 2*sqrt7)), -0.25*(-7.0 + 4*sqrt7)*(-1.0 + sqrt(7.0 + 2*sqrt7)), 0.25*(-7.0 + 4*sqrt7)*(1.0 + sqrt(7.0 + 2*sqrt7)), -0.25*(7.0 + 4*sqrt7)*(1.0 + sqrt(7.0 - 2*sqrt7)), 7.5]

			Pmn(1,:) = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
			Pmn(2,:) = [1.0, -1.0/sqrt(7.0 - 2*sqrt7), 1.0/sqrt7, 1.0/sqrt(637.0 + 238*sqrt7), 2.0/(7.0 - 5*sqrt7), 1.0/sqrt(7.0 - 0.5*sqrt7)]
			Pmn(3,:) = [1.0, -1.0/sqrt(7.0 + 2*sqrt7), -1.0/sqrt7, 1.0/sqrt(637.0 - 238*sqrt7), 2.0/(7.0 + 5*sqrt7), -1.0/sqrt(7.0 + 0.5*sqrt7)]

			Pnm(:,1) = [1.0/30.0, -0.1, 1.0/6.0, -7.0/30.0, 0.3, -1.0/6.0]
			Pnm(:,2) = [(14.0 - sqrt7)/60.0, -0.05*sqrt(49.0 + 10*sqrt7), (-1.0 + 2*sqrt7)/12.0, sqrt(931.0 - 350*sqrt7)/60.0, -0.15*(1.0 + sqrt7), sqrt(7.0 - 0.5*sqrt7)/6.0]
			Pnm(:,3) = [(14.0 + sqrt7)/60.0, -0.05*sqrt(49.0 - 10*sqrt7), (-1.0 - 2*sqrt7)/12.0, sqrt(931.0 + 350*sqrt7)/60.0, 0.15*(-1.0 + sqrt7), -sqrt(7.0 + 0.5*sqrt7)/6.0]

		case (6)
			M(1:7,1) = [4.0/91.0, (-3.0 + 7*sqrt15)/2730.0, (-3.0 - 7*sqrt15)/2730.0, 16.0/1365.0, (-3.0 - 7*sqrt15)/2730.0, (-3.0 + 7*sqrt15)/2730.0, -1.0/273.0]
			M(2:6,2) = [(744.0 - 42*sqrt15)/2275.0, 121.0/4550.0, (24.0 - 56*sqrt15)/6825.0, 121.0/4550.0, (-124.0 + 7*sqrt15)/4550.0]
			M(3:5,3) = [(744.0 + 42*sqrt15)/2275.0, (24.0 + 56*sqrt15)/6825.0, -(124.0 + 7*sqrt15)/4550.0]
			M(4,4) = 1024.0/2275.0

			Minv(1:7,1) = [24.5, (-105.0 - 245*sqrt15)/726.0, (-105.0 + 245*sqrt15)/726.0, -35.0/32.0, (-105.0 + 245*sqrt15)/726.0, (-105.0 - 245*sqrt15)/726.0, 3.5]
			Minv(2:6,2) = [1225.0*(124.0 + 7*sqrt15)/43923.0, -175.0/363.0, 175.0*(3.0 + 7*sqrt15)/11616.0, -175.0/363.0, 175.0*(124.0 + 7*sqrt15)/43923.0]
			Minv(3:5,3) = [1225.0*(124.0 - 7*sqrt15)/43923.0, 175.0*(3.0 - 7*sqrt15)/11616.0, 175.0*(124.0 - 7*sqrt15)/43923.0]
			Minv(4,4) = 1225.0/512.0

			D(1,1) = -10.5
			D(2,1:2) = [10.0/(3.0 - 7*sqrt15 + sqrt(300.0 + 26*sqrt15)), 0.0]
			D(3,1:3) = [10.0/(3.0 + 7*sqrt15 - sqrt(300.0 - 26*sqrt15)), -3.0/44.0*(124.0 + 13*sqrt11 + (10*sqrt11 - 7.0)*sqrt15)/sqrt(45.0 + 6*sqrt15), 0.0]
			D(4,1:4) = [-0.3125, 0.0625*sqrt(621.0 - 529.5*sqrt3/sqrt5), -0.0625*sqrt(621.0 + 529.5*sqrt3/sqrt5), 0.0]
			D(5,1:5) = [10.0/(3.0 + 7*sqrt15 + sqrt(300.0 - 26*sqrt15)), -3.0/44.0*(-124.0 + 13*sqrt11 + (10*sqrt11 + 7.0)*sqrt15)/sqrt(45.0 + 6*sqrt15), 0.5*sqrt(3.0 + 2*sqrt3/sqrt5), -16.0/sqrt(75.0 - 6.5*sqrt15), 0.0]
			D(6,1:6) = [10.0/(3.0 - 7*sqrt15 - sqrt(300.0 + 26*sqrt15)), 0.5*sqrt(3.0 - 2*sqrt3/sqrt5), 3.0/44.0*(-124.0 - 13*sqrt11 + (10*sqrt11 - 7.0)*sqrt15)/sqrt(45.0 - 6*sqrt15), 16.0/sqrt(75.0 + 6.5*sqrt15), -3.0/44.0*(124.0 - 13*sqrt11 + (10*sqrt11 + 7.0)*sqrt15)/sqrt(45.0 - 6*sqrt15), 0.0]
			D(7,1:7) = [0.5, 0.05*(-39.0 - 30*sqrt15 + sqrt(8955.0 + 1974*sqrt15)), 0.05*(-39.0 + 30*sqrt15 - sqrt(8955.0 - 1974*sqrt15)), -3.2, 0.05*(-39.0 + 30*sqrt15 + sqrt(8955.0 - 1974*sqrt15)), 0.05*(-39.0 - 30*sqrt15 - sqrt(8955.0 + 1974*sqrt15)), 10.5]

			Pmn(1,:) = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
			Pmn(2,:) = [1.0, -1.0/sqrt(3.0 - 2*sqrt3/sqrt5), (2.0 + sqrt15)/11.0, -7.0/sqrt(837.0 + 762*sqrt3/sqrt5), (-67.0 + 5*sqrt15)/363.0, 5.0/sqrt(621.0 - 529.5*sqrt3/sqrt5), -5.0/363.0*(3.0 + 7*sqrt15)]
			Pmn(3,:) = [1.0, -1.0/sqrt(3.0 + 2*sqrt3/sqrt5), (2.0 - sqrt15)/11.0, 7.0/sqrt(837.0 - 762*sqrt3/sqrt5), (-67.0 - 5*sqrt15)/363.0, -5.0/sqrt(621.0 + 529.5*sqrt3/sqrt5), 5.0/363.0*(-3.0 + 7*sqrt15)]
			Pmn(4,:) = [1.0, 0.0, -0.5, 0.0, 0.375, 0.0, -0.3125]

			Pnm(:,1) = [1.0/42.0, -1.0/14.0, 5.0/42.0, -1.0/6.0, 3.0/14.0, -11.0/42.0, 1.0/7.0]
			Pnm(:,2) = [(124.0 - 7*sqrt15)/700.0, -sqrt(51705.0 + 1686*sqrt15)/700.0, (13.0 + 10*sqrt15)/140.0, -sqrt(21855.0 - 4894*sqrt15)/300.0, 3.0*(-73.0 + 9*sqrt15)/700.0, 11.0/210.0*sqrt(75.0 + 6.5*sqrt15), (3.0 - 7*sqrt15)/70.0]
			Pnm(:,3) = [(124.0 + 7*sqrt15)/700.0, -sqrt(51705.0 - 1686*sqrt15)/700.0, (13.0 - 10*sqrt15)/140.0, sqrt(21855.0 + 4894*sqrt15)/300.0, -3.0*(73.0 + 9*sqrt15)/700.0, -11.0/210.0*sqrt(75.0 - 6.5*sqrt15), (3.0 + 7*sqrt15)/70.0]
			Pnm(:,4) = [128.0/525.0, 0.0, -64.0/105.0, 0.0, 144.0/175.0, 0.0, -16.0/35.0]

		case default
			error stop "gauss_lobatto_prep_table: p must be in [0,6], otherwise not implemented."

		end select

		! recover the lower triangular part by the left-triangular part (with diagonals)
		! M(i,j)==M(j,i)==M(p+2-i,p+2-j)==M(p+2-j,p+2-i)
		do j = 1, (p+1)/2
			! reflect column j
			do i = j, p+1-j
				M(p+2-j, p+2-i) = M(i,j)
				Minv(p+2-j, p+2-i) = Minv(i,j)
			end do
		end do

		do i = 1, p+1
			do j = i+1, p+1
				M(i,j) = M(j,i)
				Minv(i,j) = Minv(j,i)
				D(i,j) = -D(p+2-i, p+2-j)
			end do

			if (iand(i,1)/=0) then
				! symmetric
				do j = 1, (p+1)/2
					Pnm(i,p+2-j) = Pnm(i,j)
					Pmn(p+2-j,i) = Pmn(j,i)
				end do
			else
				! skew-symmetric
				do j = 1, (p+1)/2
					Pnm(i,p+2-j) = -Pnm(i,j)
					Pmn(p+2-j,i) = -Pmn(j,i)
				end do
			end if
		end do
	end subroutine gauss_lobatto_matrix

	pure function eye(n) result(ans)
		implicit none
		integer, intent(in) :: n
		integer :: i
		real, dimension(n,n) :: ans
		ans(:,:) = 0.0
		do i = 1, n
			ans(i,i) = 1.0
		end do
	end function eye
	
	pure subroutine mat_fast_pow(mat, n, ans)
		implicit none
		real, dimension(:,:), intent(in) :: mat
		integer, value :: n
		real, dimension(size(mat,1), size(mat,2)) :: mat_
		real, dimension(size(mat,1), size(mat,2)), intent(out) :: ans
		integer :: i, m
		m = size(mat,1)
		ans = 0.0
		do i = 1, m
			ans(i,i) = 1.0
		end do
		if (n <= 0) return
		mat_ = mat
		do while(n > 0)
			if (mod(n,2) == 1) ans = matmul(mat_, ans)
			mat_ = matmul(mat_, mat_)
			n = n/2
		end do
	end subroutine mat_fast_pow

	subroutine test_lobatto_mat(p, x, w, M, Minv, D, Pmn, Pnm)
		use, intrinsic :: iso_fortran_env, only: stdin  => input_unit, stdout => output_unit, stderr => error_unit ! need Fortran 2003
		implicit none
		integer, intent(in) :: p
		real, dimension(0:p), intent(in) :: x, w
		real, dimension(0:p, 0:p), intent(in) :: M, Minv, D, Pmn, Pnm
		real, dimension(0:p,0:p) :: D_pow, mat
		integer :: i
		real :: err
		err = 0.0
		err = max(err, maxval(abs(matmul(Pmn, Pnm) - eye(p+1)))) ! inverse relation error
		err = max(err, maxval(abs(matmul(M, Minv) - eye(p+1)))) ! inverse relation error
		D_pow = eye(p+1)
		call mat_fast_pow(D, p+1, D_pow)
		err = max(err, maxval(abs(D_pow))) ! nilpotent error

		! SBP error
		mat = matmul(M, D) + matmul(transpose(D), M)
		mat(p,p) = mat(p,p) - 1.0
		mat(0,0) = mat(0,0) + 1.0
		err = max(err, maxval(abs(mat)))

		write(stdout,*) err, maxval(abs(mat))
	end subroutine test_lobatto_mat
	
	pure elemental subroutine legendre(n, x, y, dy)
		!! compute values and derivatives of the standard Legendre polynomials
		!! recurrence: (n+1)*P_{n+1}(x) = (2n+1)*x*P_n(x) - n*P_{n-1}(x) for n>=0
		!! input:
		!!   n: degree, >=0
		!!   x: evaluation point
		!! output:
		!!   y: P_n(x)
		!!   dy: optional, P_n'(x)
		implicit none
		integer, intent(in) :: n
		real, intent(in) :: x
		real, intent(out) :: y
		real, intent(out), optional :: dy
		real :: t1, t2, alpha
		integer :: k

		if (present(dy)) then
			dy = 0.0
		end if

		if (n < 0) then
			y = 0.0
		else
			! n>=0 and k=0
			t1 = 0.0    ! p_{-2}(x)
			t2 = 0.0    ! p_{-1}(x)
			y = 1.0     ! p_0(x)
			if (present(dy) .and. iand(n,1)/=0) then
				dy = 1.0
			end if

			! n>0
			do k = 1, n
				! compute p_{k}(x), and contribute to p_{k}'(x)
				t1 = t2 ! t1 is p_{k-2}(x)
				t2 = y  ! t2 is p_{k-1}(x)
				alpha = 1.0/k
				y = (2.0-alpha)*x*t2 - (1.0-alpha)*t1 ! y is p_{k}(x)
				if (present(dy) .and. iand(n-k,1)/=0) then
					dy = dy + (2*k+1)*y
				end if
			end do
		end if
	end subroutine legendre
	
	pure subroutine gauss_legendre_matrix(p, m, minv, D, tl, tr)
		!! extract the vectors and matrices for p-th order Gauss--Legendre modal space on the 1D interval [-1,1]
		!! 
		!! input: 
		!!   p:      polynommial order
		!! 
		!! output:
		!!   m:      diagonal mass matrix entries w.r.t. modal DOFs, size = [p+1]
		!!   minv:   inverse of the mass matrix, size = [p+1]
		!!   tl:     left trace vector, size = [p+1]
		!!   tr:     right trace vector, size = [p+1]
		!!   D:      differentiation matrix (each column is the coefficient vector of a basis function), size = [p+1, p+1]
		implicit none
		integer, intent(in) :: p
		real, dimension(0:p), intent(out) :: m, minv, tl, tr
		real, dimension(0:p,0:p), intent(out) :: D
		integer :: i,j
		
		do j = 0, p
			m(j) = 1.0 / (j + 0.5)
			
			minv(j) = j + 0.5
			
			tr(j) = 1.0
			if (iand(j,1)==0) then
				tl(j) = 1.0
			else
				tl(j) = -1.0
			end if
			
			do i = 0, j-1
				if (iand(j-i, 1)==1) then
					D(i, j) = 2.0*i + 1.0
				else
					D(i, j) = 0.0
				end if
			end do
			D(j:p, j) = 0.0
		end do
	end subroutine gauss_legendre_matrix
	
	subroutine test_legendre_mat(p, m, D, tl, tr)
		!! test the SBP error of Legendre matrices
		use, intrinsic :: iso_fortran_env, only: stdin  => input_unit, stdout => output_unit, stderr => error_unit ! need Fortran 2003
		implicit none
		integer, intent(in) :: p
		real, dimension(0:p), intent(in) :: m, tl, tr
		real, dimension(0:p, 0:p), intent(in) :: D
		real, dimension(0:p, 0:p) :: D_pow, mat
		integer :: i, j
		real :: err
		
		err = 0.0
		D_pow = eye(p+1)
		call mat_fast_pow(D, p+1, D_pow)
		err = max(err, maxval(abs(D_pow))) ! nilpotent error

		! SBP error
		do j = 0, p
			do i = 0, p
				mat(i,j) = m(i) * D(i,j) + D(j,i) * m(j) - (tr(i)*tr(j) - tl(i)*tl(j))
			end do
		end do
		err = max(err, maxval(abs(mat)))

		write(stdout,*) err, maxval(abs(mat))
	end subroutine test_legendre_mat
	
	pure function interpolation_matrix(x0, x1) result(mat)
		!! generate an polynomial interpolation matrix from point set x0 to x1
		implicit none
		real, dimension(:), intent(in) :: x0, x1
		real, dimension(size(x1), size(x0)) :: mat
		integer :: m, n, i, j, k
		m = size(x1)
		n = size(x0)
		mat(:,:) = 1.0
		do j = 1, n
			do k = 1, n
				if (k /= j) then
					do i = 1, m
						mat(i,j) = mat(i,j) * ( (x1(i) - x0(k)) / (x0(j) - x0(k)) )
					end do
				end if
			end do
		end do
	end function interpolation_matrix
	
	pure function Horner_Taylor(list, x) result(ans)
		!! Horner's method / Clenshaw formula for evaluating a Taylor / monomial expansion of a polynomial
		!! list(:) should contain the coefficients of the polynomial in ascending order, i.e., list(1) is the constant term, list(2) is the coefficient of x, etc.
		implicit none
		real, dimension(:), intent(in) :: list
		real, intent(in) :: x
		real :: ans
		
		integer :: i
		ans = 0.0
		do i = size(list), 1, -1
			ans = ans * x + list(i)
		end do
	end function Horner_Taylor
	
	pure function Clenshaw_Legendre(list, x) result(b0)
		!! Horner's method / Clenshaw formula for evaluating a Legendre expansion of a polynomial
		!! list(:) should contain the coefficients of the polynomial in ascending order, i.e., list(1) is the coefficient of P_0(x), list(2) is the coefficient of P_1(x), etc.
		implicit none
		real, dimension(:), intent(in) :: list
		real, intent(in) :: x
		
		integer :: k
		real :: b1, b2, b0
		
		b0 = 0.0
		b1 = 0.0
		b2 = 0.0
		do k = size(list), 1, -1
			b2 = b1
			b1 = b0
			b0 = list(k) + ((2.0*k-1.0)/k) * x * b1 - (k/(k+1.0)) * b2
		end do
	end function Clenshaw_Legendre
end module gauss_m