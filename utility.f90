module utility
	use, intrinsic :: iso_fortran_env, only: stdin  => input_unit, stdout => output_unit, stderr => error_unit ! need Fortran 2003
	implicit none
	real, parameter :: PI_math = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
	
	character(len=80), parameter, private :: paragraph_sep = repeat('-',80)
contains
	subroutine compiler_para_output(file_unit)
		!! output the compiler information to the file represented by 'file_unit'
		use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options ! need Fortran 2008 standard support
		implicit none
		integer, intent(in) :: file_unit
		write(file_unit, '(A)') paragraph_sep
		write(file_unit, '(A)') 'compiler_version:'
		write(file_unit, '(A)') compiler_version()
		write(file_unit, '(A)') 'compiler_options:'
		write(file_unit, '(A)') compiler_options()
		write(file_unit, '(A)') paragraph_sep
	end subroutine compiler_para_output
	
	subroutine omp_para_output(file_unit)
		!! output the OpenMP parameters to the file represented by 'file_unit'
#ifdef _OPENMP
		use omp_lib
#endif
		integer, intent(in) :: file_unit
		integer :: i
#ifdef _OPENMP
		write(file_unit, '(A)') paragraph_sep
		write(file_unit, '(A, g0)') "OMP_version     = ", _OPENMP ! in YYYYMM format
		write(file_unit, '(A, g0)') "OMP_max_threads = ", omp_get_max_threads() ! max number of threads able to use
		write(file_unit, '(A, g0)') "OMP_num_procs   = ", omp_get_num_procs() ! number of visible processors
		write(file_unit, '(A, g0, A)') "OMP_proc_bind   = ", omp_get_proc_bind(), ' (0: false, 1: true, 2: master, 3: close, 4: spread)'
		write(file_unit, '(A, g0)') "OMP_num_places  = ", omp_get_num_places() ! number of places available to the exe in the place list
		write(file_unit, '(A, g0)') "OMP_dynamic     = ", omp_get_dynamic() ! create a new thread pool each time
		write(file_unit, '(A, g0)') "OMP_max_active_levels = ", omp_get_max_active_levels() ! max parallel level
		!$omp parallel do ordered schedule(static,1)
		do i = 1, omp_get_max_threads()
			!$omp ordered
			block
			integer :: affinity_info_len
			character (len=0) :: null_string
			character (len=:), allocatable :: affinity_buffer
			affinity_info_len = omp_capture_affinity(null_string, '')
			allocate(character(len=affinity_info_len) :: affinity_buffer)
			affinity_info_len = omp_capture_affinity(affinity_buffer, '')
			write(file_unit, '(A)') affinity_buffer
			deallocate(affinity_buffer)
			end block
			!$omp end ordered
		end do
		!$omp end parallel do
		write(file_unit, '(A)') paragraph_sep
#else
		write(file_unit, '(A)') paragraph_sep
		write(file_unit, '(A)') 'OpenMP is not enabled.'
		write(file_unit, '(A)') paragraph_sep
#endif
	end subroutine omp_para_output
end module utility