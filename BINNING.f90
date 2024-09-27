program BINNING
     use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
     character(len=256) :: inputfile  !specific simulation results
     character(len=256) :: outputfile  
     integer(kind=dpi) :: n_cols
     real(kind=dpr), dimension(:,:), allocatable :: data
     real(kind=dpr) :: mean, err
     integer(kind=dpi) :: n, m, i, j, k

     !get params and read data
     call PARAM_READER(inputfile, outputfile, n_cols)
     call FILE_READER(inputfile, data, n_cols)

     open(unit=26, file=outputfile, action='write', status='replace')
     ! Itera sobre cada columna
     do j = 1, n_cols
       write(26,*) 'Columna', j
       write(26,*) 'Tamaño bloque Mediana error'
       ! Itera sobre los tamaños de bloque
       k = dint(log(0.5_dpr*n)/log(2.0_dpr))
       do i = 0, k
         m = 2**i
         call BLOCK_AVERAGE(data(:, j), m, mean, err)
         write(26,*) m, mean, err
       end do
       write(26,*) ''
     end do
     close(26)

 contains

  !!      READERS       !!
  !read file with data organized in columns (magnitudes) at each time(row)
  subroutine FILE_READER(file2read, data, n)
     use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
     character(len=256), intent(in) :: file2read  !gnral file containing data from simulation
     real(kind=dpr), dimension(:,:), allocatable :: data
     integer(kind=dpi), intent(out) :: n
     !local
     real(kind=dpr) :: x, mi
     integer(kind=dpi) :: io, i

     !get n lines
     n=0
     open(unit=25, file=file2read, action='read', status='old')
     iterate: do
       read(25, *, iostat=io) x, mi
         if (io .LT. 0) then
                 exit iterate
         end if
         n= n+1
     end do iterate

     !allocate + re-read to get data
     allocate(data(n, n_cols))
     rewind(25)
     do i=1, n
       read(25,*) data(i,:)
     end do
     close(25)
 end subroutine FILE_READER

 !Parameters reader from user input
 subroutine PARAM_READER(inputfile, outputfile, n_cols)
     use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
     character(len=256), intent(out) :: inputfile  !specific simulation results
     character(len=256), intent(out) :: outputfile !specific simulation results ONCE BLOCK AVERAGED
     integer(kind=dpi), intent(out) :: n_cols    !n columns used to BLOCK AVERAGE
     !local
     character(len=256)::input 


     write(*,*) 'Write the main name (e.q: simulat1) of the file (w.o. extension) to perform BLOCK AVERAGE'
     write(*,*) 'and the numbers of bins /columns for statistics separated by a spacebar' 
     write(*,*) 'DATA will be read from the simulat1.dat file and '
     write(*,*)  'AVEAGED RESULTS written in the simulat1_AVERAGED.dat file'
     read(*,*)  input ,  n_cols

     inputfile= trim(input)//".dat"
     outputfile= trim(input)//"_AVERAGED.dat"
 end subroutine PARAM_READER

 !!   BLOCK AVERAGE METHOD for a single column   !!
 subroutine BLOCK_AVERAGE(column_data, m, mean, err)
     use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
     real(kind=dpr), dimension(:), intent(in) :: column_data
     integer(kind=dpi), intent(in) :: m
     real(kind=dpr), intent(out) :: mean, err
     !local
     real(kind=dpr), dimension(:), allocatable :: block_vals
     integer(kind=dpi) :: N, block_num, i

     N=size(column_data)
     block_num= ceiling(dble(N/m))
     allocate(block_vals(block_num))
     err= 0.0_dpr
     do i=1, block_num
       block_vals(i) = sum(column_data((i-1)*m+1: min(i*m,N)))/dble(min(m,N))
     end do

     mean=sum(block_vals)/block_num
     do i=1, block_num
       err= err+ (block_vals(i)-mean)**2
     end do
     err= dsqrt(err/(block_num*(block_num-1)))
end subroutine BLOCK_AVERAGE

end program BINNING
