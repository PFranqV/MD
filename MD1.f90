 program MD1
  use, intrinsic :: iso_fortran_env, only : dpr=> real64, dpi=> int64
  use MDmod

  integer(kind=dpi) :: N
  real(kind=dpr) :: rho, L, a, Rc, T, dt, V_LJ, KE 
  integer(kind=dpi) :: n_step, n_distr, writestep, i
  real(kind=dpr), dimension(:,:), allocatable :: pos, vels
  real(kind=dpr), dimension(:,:), allocatable :: V_distr

  !!! DATA of the system !!!
  N=125               !N particles
  rho= 0.7            !reduced density
  T= 100.0            !kbT
  dt= 1e-4                !timestep
  Rc= 2.5             !cutoff radii
  L= (N/rho)**(1.0_dpr/3) !length of simulation box

  allocate(pos(N,3), vels(N,3))

  !!! INITIALIZATION SC lattice + initial Vels !!!
  call SC_INIT(N, L, rho, pos, a)
  if (Rc .LT. a) Rc= 1.25_dpr*a   !ensure min. conv. img.
  call BIM_DISTR(T, vels, KE)
 !call Vxyz_DISTR(vels,20_dpi,'MD1_VxyzDISTR_init.dat')
  
  open(unit=20, file='MD1termo_VVerlet.dat', action='write', status='replace')
  write(20,*) '#time     Pressure       T        EK        V_LJ        Et'
  write(20,*) 0.0_dpr, f_Pmom(vels), 2*KE/(3*N), KE, V_LJ, V_LJ+KE
  n_step=2e5
  writestep= 100
  n_distr= 1000
  allocate(V_distr(n_distr*N,3))
  do i=1, n_step
    call VVERLET(L,Rc,dt,pos,vels,V_LJ)
    if (modulo(i,writestep) .EQ. 0) then
            call KINETIC(vels,KE)
            write(20,*) i*dt, f_Pmom(vels), 2*KE/(3*N), KE, V_LJ, KE+V_LJ
    end if
    !close to the final, save v values for statistics
    if (i .GE. n_step-n_distr) then
            V_distr((i-n_step+n_distr+1)*N: (i-n_step+n_distr+2)*N, :)= vels
    end if
  end do
  call VNORM_DISTR(V_distr,20_dpi,'MD1_VnormDISTR_fin.dat')
  call Vxyz_DISTR(V_distr,20_dpi,'MD1_VxyzDISTR_fin.dat')
  close(20)
  write(*,*) '*** FINISHED ***'

 end program MD1
