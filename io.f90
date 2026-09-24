module io
	implicit none
contains
	subroutine save_tecbin_finite_volume(u_avg, nx, ny, grid_x, grid_y, id, time)
		!! output DG-MHD computation result in finite volume style (tecplot binary format .plt)
		use, intrinsic :: iso_fortran_env, only: int32, real32, real64
		use MHD_PDE_m, only: MHD_con2pri
		integer,intent(in)	:: nx, ny, id
		real,	intent(in)	:: grid_x(nx+1), grid_y(ny+1), time
		real,	intent(in)	:: u_avg(8, 0:nx+1, 0:ny+1)
		integer, parameter :: ncoor = 2, nvar = 8
		
		integer(kind=int32), parameter :: zero=0, one=1, neg_one=-1, two=2, neg_two=-2 ! use automatic type casting
		character(len=100) :: filename, zonename
		real(kind=real32), dimension(:,:,:), allocatable :: buffer_xy, buffer_var
		real(kind=real64), dimension(nvar) :: var_min, var_max
		integer :: file_unit, ios, i, j, k, l
		
		write(filename, '(A, "_", g0, ".plt")') 'solution', id
		open(newunit=file_unit, file=trim(filename), form='unformatted', status='replace', access='stream', iostat=ios)
		if (ios/=0) error stop 'save_result_tecbin: failed to open file!'
		
		! HEADER SECTION
		write(file_unit) "#!TDV191"
		write(file_unit) one ! Integer value of 1. This is used to determine the byte order of the reader, relative to the writer.
		write(file_unit) zero ! FileType: 0 = FULL, 1 = GRID, 2 = SOLUTION
		call write_string("MHD result") ! The TITLE.
		write(file_unit) int(ncoor + nvar, kind=int32)
		call write_string("x")
		call write_string("y")
		call write_string("Density")
		call write_string("Vx")
		call write_string("Vy")
		call write_string("Vz")
		call write_string("Bx")
		call write_string("By")
		call write_string("Bz")
		call write_string("Pressure")
		write(file_unit) real(299.0, kind=real32) ! Zone marker for v112 or v191
		write(zonename, '(A, "_", g0)') 'zone', id
		call write_string(zonename)
		write(file_unit) neg_one ! no longer used parent zone marker, set to -1
		write(file_unit) neg_two ! strand ID
		write(file_unit) real(time, kind=real64) ! solution time (double)
		write(file_unit) neg_one ! defautl zone color
		write(file_unit) zero ! zone type, 0 = ordered
		write(file_unit) one  ! 1 = specify data location
		do k = 1, ncoor
			write(file_unit) zero ! x and y data are at nodes
		end do
		do k = 1, nvar
			write(file_unit) one ! vars are cell centered
		end do
		write(file_unit) zero ! ordered data must specify 0 for raw local 1-to-1 face connection
		write(file_unit) zero ! number of misc user-defined face neighbor connections
		write(file_unit) int(nx+1, kind=int32) ! maxI
		write(file_unit) int(ny+1, kind=int32) ! maxJ
		write(file_unit) one ! maxK
		write(file_unit) zero ! 0 = no more auxiliary name/value pairs
		! no geometry
		! no text
		! no custom label
		! no user rec
		! no dataset auxiliary data
		! no variable auxiliary data
		
		! EOHMARKER
		write(file_unit) real(357.0, kind=real32)
		
		
		allocate(buffer_xy(nx+1, ny+1, ncoor), buffer_var(nx+1, ny, nvar))
		var_min(:) = +huge(1.0_real64)
		var_max(:) = -huge(1.0_real64)
		do j = 1, ny+1
			do i = 1, nx+1
				buffer_xy(i,j,1) = grid_x(i)
				buffer_xy(i,j,2) = grid_y(j)
			end do
		end do
		do j = 1, ny
			do i = 1, nx
				buffer_var(i,j,:) = MHD_con2pri(u_avg(:,i,j))
				do l = 1, nvar
					var_min(l) = min(var_min(l), real(buffer_var(i,j,l), kind=real64))
					var_max(l) = max(var_max(l), real(buffer_var(i,j,l), kind=real64))
				end do
			end do
		end do
		
		! DATA SECTION
		write(file_unit) real(299.0, kind=real32) ! Zone marker for v112
		do i = 1, ncoor + nvar
			write(file_unit) one ! 1 = single precision float
		end do
		write(file_unit) zero ! 0 = no passive variables
		write(file_unit) zero ! 0 = no variable sharing
		write(file_unit) neg_one ! -1 = no sharing connectivity list
		write(file_unit) real(minval(grid_x), kind=real64)
		write(file_unit) real(maxval(grid_x), kind=real64)
		write(file_unit) real(minval(grid_y), kind=real64)
		write(file_unit) real(maxval(grid_y), kind=real64)
		do i = 1, nvar
			write(file_unit) var_min(i), var_max(i)
		end do
		write(file_unit) buffer_xy
		write(file_unit) buffer_var ! 少了一个都不行，说明size(nx+1,ny,nvar)是正好的！
		close(file_unit)
		deallocate(buffer_xy, buffer_var)
	contains
		subroutine write_string(string)
			character(len=*), intent(in) :: string
			integer(kind=int32), dimension(len_trim(string)+1) :: ASCII_buffer
			integer :: i
			do i = 1, len_trim(string)
				ASCII_buffer(i) = iachar(string(i:i), kind=int32)
			end do
			ASCII_buffer(len_trim(string)+1) = int(0, kind=int32)
			write(file_unit) ASCII_buffer
		end subroutine write_string
	end subroutine save_tecbin_finite_volume
    
    subroutine save_tecbin_entlog(my_log, n, filename, zonename)
		!! output DG-MHD computation log of (time, entropy) (tecplot binary format .plt)
		use, intrinsic :: iso_fortran_env, only: int32, real32, real64
		integer,intent(in)	:: n
		real,	intent(in)	:: my_log(:,:) ! the first dimension must be nvar
		integer,  parameter :: nvar = 2
        character(len=*), intent(in) :: filename, zonename
		
		integer(kind=int32), parameter :: zero=0, one=1, neg_one=-1, two=2, neg_two=-2 ! use automatic type casting
		real(kind=real64), dimension(nvar) :: var_min, var_max
		integer :: file_unit, ios, i, j
		
		! write(filename, '(A, ".plt")') 'entropy'
		open(newunit=file_unit, file=trim(filename), form='unformatted', status='replace', access='stream', iostat=ios)
		if (ios/=0) error stop 'save_result_tecbin: failed to open file!'
		
		! HEADER SECTION
		write(file_unit) "#!TDV191"
		write(file_unit) one ! Integer value of 1. This is used to determine the byte order of the reader, relative to the writer.
		write(file_unit) zero ! FileType: 0 = FULL, 1 = GRID, 2 = SOLUTION
		call write_string("MHD entropy log") ! The TITLE.
		write(file_unit) int(nvar, kind=int32)
		call write_string("t")
		call write_string("total entropy")
		write(file_unit) real(299.0, kind=real32) ! Zone marker for v112 or v191
		! write(zonename, '(A, "_", g0)') 'zone', 0 ! zone title
		call write_string(zonename)
		write(file_unit) neg_one ! no longer used parent zone marker, set to -1
		write(file_unit) neg_two ! strand ID
		write(file_unit) real(0.0, kind=real64) ! (pseudo) solution time (double)
		write(file_unit) neg_one ! defautl zone color
		write(file_unit) zero ! zone type, 0 = ordered
		write(file_unit) one  ! 1 = specify data location
		do i = 1, nvar
			write(file_unit) zero ! vars are at nodes
		end do
		write(file_unit) zero ! ordered data must specify 0 for raw local 1-to-1 face connection
		write(file_unit) zero ! number of misc user-defined face neighbor connections
		write(file_unit) int(n, kind=int32) ! maxI
		write(file_unit) one ! maxJ
		write(file_unit) one ! maxK
		write(file_unit) zero ! 0 = no more auxiliary name/value pairs
		! no geometry
		! no text
		! no custom label
		! no user rec
		! no dataset auxiliary data
		! no variable auxiliary data
		
		! EOHMARKER
		write(file_unit) real(357.0, kind=real32)
		
		! DATA SECTION
		write(file_unit) real(299.0, kind=real32) ! Zone marker for v112
		do i = 1, nvar
			write(file_unit) one ! 1 = single precision float
		end do
		write(file_unit) zero ! 0 = no passive variables
		write(file_unit) zero ! 0 = no variable sharing
		write(file_unit) neg_one ! -1 = no sharing connectivity list
        var_min(:) = minval(my_log(:, 1:n), 2)
        var_max(:) = maxval(my_log(:, 1:n), 2)
        do i = 1, nvar
            write(file_unit) var_min(i), var_max(i)
        end do
        do i = 1, nvar
            do j = 1, n
                write(file_unit) real(my_log(i,j), kind=real32)
            end do
        end do
		close(file_unit)
	contains
		subroutine write_string(string)
			character(len=*), intent(in) :: string
			integer(kind=int32), dimension(len_trim(string)+1) :: ASCII_buffer
			integer :: i
			do i = 1, len_trim(string)
				ASCII_buffer(i) = iachar(string(i:i), kind=int32)
			end do
			ASCII_buffer(len_trim(string)+1) = int(0, kind=int32)
			write(file_unit) ASCII_buffer
		end subroutine write_string
	end subroutine save_tecbin_entlog
end module io