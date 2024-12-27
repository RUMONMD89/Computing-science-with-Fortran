! _/_/_/ Finite Difference Method for 3D Thermal Conduction in Eulerian _/_/_/
! 2022.08.23-29 Written by Y.Hirokawa
!
! Governing Eauation: cv*dT/dt = rambda*d^2T/dx^2
! IMAT: 0=SUS304, 1=Cu

program main
  implicit none
  integer, parameter :: IDIM=10, JDIM=10, KDIM=10, MAXSTEP=6000, NFREQ=600, IMAT=1
  integer            :: istep
  double precision   :: t0, th, tl
  double precision   :: cv, dt, dx, dy, dz, rambda
  double precision   :: t(IDIM,JDIM,KDIM)

  ! Initialization
  write(*,*) "[INFO] Initialze Field..."
  call sub_init(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM, IMAT, istep)

  ! Output Initial States
  write(*,*) "[INFO] File Output of Initialze Field..."
  call sub_output(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM, IMAT, istep)

  ! Time Marching
  do istep = 1, MAXSTEP
    ! Set Boundary Condition
    call sub_boundary(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM)

    ! Calculate the Equation
    call sub_calc(t, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM)

    if(mod(istep, NFREQ) == 0  .or.  istep == MAXSTEP) then
      ! Output Initial States
      write(*,*) "[INFO] Step =", istep, ", Time =", dble(istep)*dt, " [s]"
      call sub_output(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM, IMAT, istep)
    endif
  enddo

  write(*,*) "[INFO] Successfully Completed. "

contains

  subroutine sub_init(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM, IMAT, istep)
    integer, intent(in)             :: IDIM, JDIM, KDIM, IMAT
    integer, intent(inout)          :: istep
    double precision, intent(inout) :: t(IDIM,JDIM,KDIM), t0, th, tl, cv, rambda, dt, dx, dy, dz
    double precision                :: rho, sh, cdiff, ds

    if(IMAT == 0) then
      ! Quantities of SUS304
      ! https://www.jssa.gr.jp/contents/about_stainless/key_properties/comparison/
      rho    = 7.93d3  ! [kg/m^3]
      sh     = 0.50d3  ! [J/kg・K]
      cv     = sh*rho  ! [J/m^3・K]
      rambda = 1.60d1  ! [W/m・K]
      write(*,*) "[INFO] Material = SUS304"
    else
      ! Quantities of Cu
      ! https://www.hakko.co.jp/qa/qakit/html/h01020.htm
      rho    = 8.96d3  ! [kg/m^3]
      sh     = 3.85d2  ! [J/kg・K]
      cv     = sh*rho  ! [J/m^3・K]
      rambda = 3.86d2  ! [W/m・K]
      write(*,*) "[INFO] Material = Cu"
    endif

    ! Discrete Time
    dt = 1.0d-1     ! [s]

    ! Discrete Space
    dx = 5.0d-2      ! [m]
    dy = 5.0d-2      ! [m]
    dz = 5.0d-2      ! [m]

    ! Diffusion Number Check for Numerical Stability
    ds = min(dx, dy, dz)
    cdiff = (rambda/cv)*dt/(ds**2)
    write(*,*) "[INFO] Diffusion Number =", cdiff

    ! Initial Condition
    t0 = 3.0d2
    t(:,:,:) = t0

    ! Timestep
    istep = 0

    ! Boundary Condition for Hot, Cool
    th = 4.0d2
    tl = 3.0d2

    return
  end subroutine sub_init


  subroutine sub_boundary(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM)
    integer, intent(in)             :: IDIM, JDIM, KDIM
    double precision, intent(inout) :: t(IDIM,JDIM,KDIM), t0, th, tl, cv, rambda, dt, dx, dy, dz
    double precision                :: rho, sh
    integer                         :: i, j ,k

    ! Finite Difference Method (1st order, Explicit)

    ! Left: Adiabatic boundary (Neumann boudary)
    t(1,:,:) = t(2,:,:)

    ! Right: Adiabatic boundary (Neumann boudary)
    t(IDIM,:,:) = t(IDIM-1,:,:)

    ! Front: Adiabatic boundary (Neumann boudary)
    t(:,1,:) = t(:,2,:)

    ! Rear: Adiabatic boundary (Neumann boudary)
    t(:,JDIM,:) = t(:,JDIM-1,:)

    ! Bottom: Dirichlet boundary
    k = 1
    t(:,:,k) = th

    ! Top: Dirichlet boundary
    k = KDIM
    t(:,:,k) = tl

    return
  end subroutine sub_boundary


  subroutine sub_calc(t, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM)
    integer, intent(in)             :: IDIM, JDIM, KDIM
    double precision, intent(inout) :: t(IDIM,JDIM,KDIM), cv, rambda, dt, dx, dy, dz
    double precision                :: rho, sh, to(IDIM,JDIM,KDIM)
    integer                         :: i, j ,k

    ! Backup Current Data
    to(:,:,:) = t(:,:,:)

    ! Finite Difference Method (1st order, Explicit)
    do k = 2, KDIM-1
      do j = 2, JDIM-1
        do i = 2, IDIM-1
          t(i,j,k) = to(i,j,k) + (                                               &
                         (to(i+1,j  ,k  ) - 2*to(i,j,k) + to(i-1,j  ,k  ))/dx**2 &
                       + (to(i  ,j+1,k  ) - 2*to(i,j,k) + to(i  ,j-1,k  ))/dy**2 &
                       + (to(i  ,j  ,k+1) - 2*to(i,j,k) + to(i  ,j  ,k-1))/dz**2 &
                     )*(rambda*dt/cv)
        enddo
      enddo
    enddo

    return
  end subroutine sub_calc


  subroutine sub_output(t, t0, th, tl, cv, rambda, dt, dx, dy, dz, IDIM, JDIM, KDIM, IMAT, istep)
    integer, intent(in)             :: IDIM, JDIM, KDIM, IMAT, istep
    double precision, intent(inout) :: t(IDIM,JDIM,KDIM), t0, th, tl, cv, rambda, dt, dx, dy, dz
    double precision                :: rho, sh, x, y, z, time
    integer                         :: i, j ,k, iu
    character                       :: cfname*19, cstep*8, cprefix*7, cextension*4

    ! Casting integer to character
    write(cstep,'(I8.8)') istep
    ! Prefix
    if(IMAT == 0) then
      cprefix = 'SUS304_'
    else
      cprefix = 'Cu_'
    endif

    ! Extension
    cextension = '.csv'

    ! File Name
    cfname = trim(cprefix) // adjustl(cstep) // cextension

    ! File Open (Reserved Unit Number = 5:Keyboard, 6:Screen)
    iu = 10
    open(unit=iu, file=cfname, form='formatted', status='replace')
    write(iu,*) "time, x, y, z, t"

    ! Calculate Time (Step starts from 0)
    time = istep*dt

    ! Calculate Position and Write Data (i,j,k starts from 1)
    do k = 1, KDIM
      do j = 1, JDIM
        do i = 1, IDIM
          x = (i-1)*dx
          y = (j-1)*dy
          z = (k-1)*dz
          write(iu,*) time,', ',x,', ',y,', ',z,', ',t(i,j,k)
        enddo
      enddo
    enddo
    close(iu)

    return
  end subroutine sub_output

end program main
