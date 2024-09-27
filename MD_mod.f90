 module MDmod

      use, intrinsic :: iso_fortran_env, only: dpr => real64, dpi=> int64
 contains
        !!!! INITIALIZATION OF SYSTEM !!!! 
 !! SC lattice initialization !!
 subroutine SC_INIT(N, L, rho, pos, a)
      use, intrinsic :: iso_fortran_env, only: dpr => real64, dpi=> int64
      integer(kind=dpi), intent(in) :: N
      real(kind=dpr), intent(in) :: L, rho
      real(kind=dpr), dimension(:,:), intent(out) :: pos
      real(kind=dpr), intent(out) :: a
      !local
      integer(kind=dpi) :: num, M, i, j, k
      real(kind=dpr) :: x, y

      M= int((N/rho)**(1.0_dpr/3))
      a = L/M
      !1st particle in 0,0,0
      do i=1, M
        x = a*(i-1)
        do j=1, M
          y= a*(j-1)
          do k=1, M
            num= M**2 *(i-1)+ M*(j-1) +k
            pos(num,:) = (/x,y,a*(k-1)/)
          end do
        end do
      end do
 end subroutine SC_INIT

 !! FCC lattice inicialization !!
 subroutine FCC_INIT(N, L, rho, pos, a)
      use, intrinsic :: iso_fortran_env, only: dpr => real64, dpi=> int64
      real(kind=dpr), intent(in) :: L, rho
      integer(kind=dpi), intent(in) :: N
      real(kind=dpr), dimension(:,:), intent(inout) :: pos
      real(kind=dpr), intent(out) :: a
      ! local
      real(kind=dpr) :: x, y, z
      integer(kind=dpi) :: i, j, k, num, M

      M= nint((N/4)**(1.0_dpr/3))
      a= L/M
      do i=1, M
        x= a*(i-1)
        do j=1,M
          y= a*(j-1)
          do k=1,M
            z= a*(k-1)
            num= M**2*(i-1)+ M*(j-1) +k
            pos(num,:)= (/x, y, z/)
            pos(M**3+num,:) = (/x+a/2, y+a/2, z/)
            pos(2*M**3+num,:) = (/x, y+a/2, z+a/2/)
            pos(3*M**3+num,:) = (/x+a/2, y, z+a/2/)
          end do
        end do
      end do
end subroutine FCC_INIT

 !! PBC implementation !!
 subroutine PBC(pos,L)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), intent(in) :: L
      real(kind=dpr), dimension(:,:), intent(inout) :: pos
      !local
      integer(dpi) ::  N, D, i, j

      N= size(pos(:,1))
      D= size(pos(1,:))
      do i=1, N
        do j=1, D
          if (pos(i,j) .GT. L/2) then
                  pos(i,j) = pos(i,j) - L
          else if (pos(i,j) .LT. -L/2) then
                  pos(i,j) = pos(i,j) + L
          end if
        end do
      end do
 end subroutine PBC

 !! BIMODAL DISTRIBUTION INITIAL VELOCITY for v_tot=0 !!
 subroutine BIM_DISTR(T, vels, KE)
      use, intrinsic :: iso_fortran_env, only: dpr => real64, dpi=> int64
      real(kind=dpr), intent(in) :: T
      real(kind=dpr), dimension(:,:), intent(inout) :: vels
      real(kind=dpr), intent(out) :: KE
      !local
      real(kind=dpr), dimension(size(vels(1,:))) :: aleat
      integer(kind=dpi) :: N, i

      N= size(vels(:,1))
      do i=1,N
        call random_number(aleat)
        aleat = aleat - 0.5_dpr
        vels(i,:) = dsign(sqrt(T),aleat)
      end do
      KE = 0.5_dpr*sum(vels*vels)
 end subroutine BIM_DISTR

        !!!! ----------------------------------------------
        !!!!           FORCES, EoM INTEGRATORS         !!!!

 !! LJ_FORCES !!
 subroutine LJ_FORCES(pos, L, Rc, V_LJ, F)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:,:), intent(in) :: pos
      real(kind=dpr), intent(in) :: L, Rc
      real(kind=dpr), intent(out) :: V_LJ
      real(kind=dpr), dimension(size(pos(:,1)), size(pos(1,:))), intent(out) :: F
      !local
      real(kind=dpr) :: d_ij2, d_ij4, d_ij8, d_ij12, Rc2, V_Rc, F_ij
      real(kind=dpr), dimension(1,size(pos(1,:)))  :: R_ij

      N=size(pos(:,1))
      Rc2= Rc*Rc
      V_LJ= 0.0_dpr
      F= 0.0_dpr
      V_Rc = 4.0_dpr*(1/Rc**12 - 1.0_dpr/Rc**6)
      do i=1, N
        do j=i+1, N
          R_ij(1,:)= pos(i,:)- pos(j,:)
          call PBC(R_ij,L)
          d_ij2= sum(R_ij*R_ij)
          if (d_ij2 .LT. Rc2) then
                  d_ij4=d_ij2*d_ij2
                  d_ij8=d_ij4*d_ij4
                  d_ij12=d_ij8*d_ij4
                  F_ij= 48.0_dpr/(d_ij8*d_ij4*d_ij2)- 24.0_dpr/d_ij8
                  !pairwise sum -> substract in 2nd part so sum over N counts once
                  F(i,:) = F(i,:) + F_ij*R_ij(1,:)
                  F(j,:) = F(j,:) - F_ij*R_ij(1,:)
                  V_LJ= V_LJ+ 4.0_dpr/d_ij12- 4.0_dpr/(d_ij4*d_ij2)- V_Rc
          end if
        end do
      end do
 end subroutine LJ_FORCES

 !! Equations of Motion INTEGRATORS !!
 ! Euler
 subroutine EULER(L, Rc, dt, pos, vels, V_LJ)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), intent(in) :: L, Rc, dt
      real(kind=dpr), dimension(:,:), intent(inout) :: pos, vels
      real(kind=dpr), intent(out) :: V_LJ
      !local
      real(kind=dpr), dimension(size(pos(:,1)), size(pos(1,:))) :: F
      integer(kind=dpi) :: i, j

      call LJ_FORCES(pos, L, Rc, V_LJ, F)
      pos= pos+ vels*dt+ 0.5_dpr*F*dt*dt
      call PBC(pos,L)
      vels= vels+ F*dt
 end subroutine EULER
  
 ! Velocity Verlet
 subroutine VVERLET(L, Rc, dt, pos, vels, V_LJ)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), intent(in) :: L, Rc, dt
      real(kind=dpr), dimension(:,:), intent(inout) :: pos, vels
      real(kind=dpr), intent(out) :: V_LJ
      !local
      real(kind=dpr), dimension(size(pos(:,1)), size(pos(1,:))) :: F
      integer(kind=dpi) :: i, j

      call LJ_FORCES(pos, L, Rc, V_LJ, F)
      pos= pos+ vels*dt+ 0.5_dpr*dt*dt*F
     !write(*,*) pos
      call PBC(pos, L)    !Evaluated on L for continuity of the system
      vels= vels+ 0.5_dpr*dt*F
     !write(*,*) vels
      call LJ_FORCES(pos, L, Rc, V_LJ, F)
      vels= vels+ 0.5_dpr*dt*F
      !write(*,*) vels
 end subroutine VVERLET

        !!!! ----------------------------------------------
        !!!!                   STATISTICS                !!

 subroutine VNORM_DISTR(vels, nbox, fileout)
 !module distribution
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:,:), intent(in) :: vels
      integer(kind=dpi), intent(in) :: nbox
      character(len=*), intent(in) :: fileout
      !local
      real(kind=dpr), dimension(size(vels(:,1))) :: V_n
      real(kind=dpr) :: widthbox
      integer(kind=dpi), dimension(nbox) :: V_distr
      integer(kind=dpi) :: N, num, i

      V_distr= 0
      N= size(vels(:,1))
      do i=1, N
        V_n(i) = dsqrt(sum(vels(i,:)**2))
      end do
      !classify between nbox intervals
      widthbox= (maxval(vels)*1.01)/dble(nbox)
      do i=1, N
        num= 1+ int(V_n(i)/widthbox)
        V_distr(num) = 1+ V_distr(num)
      end do
      !WRITE
      open(21, file=fileout, action='write', status='replace')
      do i=1, nbox
        write(21,*) widthbox*(i-1+0.5), V_distr(i)
      end do
      close(21)
 end subroutine VNORM_DISTR

 !components distribution
 subroutine Vxyz_DISTR(vels, nbox, fileout1)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:,:), intent(in) :: vels
      integer(kind=dpi), intent(in) :: nbox
      character(len=*), intent(in) :: fileout1
      !local
      real(kind=dpr) :: widthbox
      integer(kind=dpi), dimension(nbox,3) :: V_distr
      integer(kind=dpi) :: N, num, comp, i, j

      N= size(vels(:,1))
      comp= size(vels(1,:))
      widthbox= 2.01_dpr*maxval(abs(vels))/nbox
      V_distr= 0
      do i=1, N
        do j=1, comp
          if (vels(i,j) .LT. 0.0) then  !bin shift so that midbin approx v=0
                  num= 1+ int(vels(i,j)/widthbox)+ 0.5*nbox
          else
                  num= int(vels(i,j)/widthbox)+ 0.5*nbox
          end if
          V_distr(num,j)= V_distr(num,j)+1
        end do
      end do
      !WRITE
      open(unit=22, file=fileout1, action='write', status='replace')
      do i=1, nbox
        write(22,*) 1.05_dpr*widthbox*(i-0.5*nbox-0.5), V_distr(i,1), V_distr(i,2), V_distr(i,3)
      end do
      close(22)
      
 end subroutine Vxyz_DISTR

 !! Box-Muller Transformation: Uniform to Normal distr
 subroutine BOX_MULLER(vel, sigma)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:), intent(inout) :: vel
      real(kind=dpr), intent(in) :: sigma
      !local
      real(kind=dpr), dimension(2) :: U_aleat   !2 numbers from an uniform distribution
      real(kind=dpr), parameter :: pi=acos(-1.0_dpr)
      integer(kind=dpi) :: i, comp

      comp=size(vel)
      call random_number(U_aleat)
      do i=1, comp
        !We use only 1 generated random number on the N distr
        vel(i)= sigma*dsqrt(-2.0_dpr*dlog(U_aleat(1)))*dcos(2*pi*U_aleat(2))
        !vel(i)= sigma*dsqrt(-2.0_dpr*dlog(U_aleat(1)))*dsin(2*pi*U_aleat(2))
      end do
 end subroutine BOX_MULLER


        !!!!----------------------------------------------
        !!!!              THERMALIZATION              !!!!
 subroutine ANDERSEN(prob, sigma, vels)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), intent(in) :: prob, sigma
      real(kind=dpr), dimension(:,:), intent(inout) :: vels
      !local
      real(kind=dpr)  :: Rand
      integer(kind=dpi) :: N, i

      N=size(vels(:,1))
      do i=1, N
        call random_number(Rand)
        if (Rand .LT. prob) then
                call BOX_MULLER(vels(i,:), sigma)
        end if
      end do
 end subroutine ANDERSEN

 subroutine KINETIC(vels,KE)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:,:), intent(in) :: vels
      real(kind=dpr), intent(out) :: KE

      KE= 0.5_dpr*sum(vels*vels)
end subroutine KINETIC

        !!!!----------------------------------------------
        !!!!                  FUNCTIONS                !!!!

 function f_Pressure(L, Rc, pos) result (Pressure)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), intent(in) :: L, Rc
      real(kind=dpr) :: Pressure
      real(kind=dpr), dimension(:,:), intent(in) :: pos
      !local
      real(kind=dpr) :: Rc2, d_ij2, d_ij4, d_ij8, d_ij14
      real(kind=dpr), dimension(1,size(pos(1,:))) :: R_ij
      integer(kind=dpi) :: N, i, j

      Rc2= Rc*Rc
      N= size(pos(:,1))
      Pressure= 0.0_dpr
      !Pairwise sum
      do i=1, N
        do j=i+1, N
          R_ij(1,:)= pos(i,:)- pos(:,j)
          call PBC(R_ij,L)  !min image convention
          d_ij2=sum(R_ij*R_ij)
          if (d_ij2 .LT. Rc) then
                  d_ij4= d_ij2*d_ij2
                  d_ij8= d_ij4*d_ij4
                  d_ij14= d_ij8*d_ij4*d_ij2
                  Pressure= Pressure+ (48.0_dpr/d_ij14- 24.0_dpr/d_ij8)*sum(R_ij*R_ij)
          end if
        end do
      end do
 end function f_Pressure
 
 function f_Pmom(vels) result (Pmom)
      use, intrinsic :: iso_fortran_env, only: dpr=> real64, dpi=> int64
      real(kind=dpr), dimension(:,:), intent(in) :: vels
      real(kind=dpr) :: Pmom
      !local
      real(kind=dpr), dimension(size(vels(1,:))) :: P_comp
      integer(kind=dpi) :: comp, i

      comp= size(vels(1,:))
      P_comp= 0.0
      do i=1, comp
        P_comp(i)= sum(vels(:,i)*vels(:,i))
      end do
      Pmom = dsqrt(sum(P_comp))
 end function f_Pmom

 end module MDmod
