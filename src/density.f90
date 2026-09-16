
subroutine density


! ******************************************************
! ***   FIND DENSITY AND CURRENT IN PHYSICAL SPACE   ***
! ******************************************************
!
! This subroutine integrates over momentum space
! to find rho and curr.  These quantities are
! defined as:
!
!                /
! rho  =  1/r**2 | f dp 
!                /
!
!                /
! curr =  1/r**2 | p f dp 
!                /
!

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor,average_rho,mass
  real(8) :: cutoff_rho,cutoff_avg,sval
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)

  character(20) filename ! Name of outupt file.


  smallpi = acos(-1.0d0)


! Zero Angular Momentum
  if (Lfix == 0.0d0) then

     factor = 1.0d0

! Include Angular Momentum
  else

     factor = 2.0*smallpi*Lfix*drc*dpc!0.25D0/smallpi*drc*dpc

  endif  

  rho = 0.D0
  avg_rho = 0.D0
  curr = 0.D0

! The B-spline shape functions have compact support, a few cells wide,
! so only particles lying near grid point "i" can contribute to it.
! build_cell_list (utils.f90) groups the particles by radial cell in
! O(Nr+Npart), which turns the deposit below from a scan over all
! Npart particles per grid point into a scan over a few neighbouring
! cells.

  call build_cell_list(cell_start,particle_order)

  cutoff_rho = (dble(bsplineorder)+1.0d0)*drc
  cutoff_avg = (dble(bsplineorder)+1.0d0)*dr
  Wcell = ceiling(max(cutoff_rho,cutoff_avg)/dr) + 1

! The parallel loop runs over "i" only, never collapsed with "j":
! rho(i), curr(i) and avg_rho(i) accumulate over all particles, so one
! grid point must belong to a single thread for the whole inner loop
! to keep the accumulation race-free without atomics.
!
! The shape function depends on r_part alone, so it is evaluated once
! into "sval" and reused by the density and the current.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,sval)

  do i=1,Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        if (abs(r(i)-r_part(j))<=cutoff_rho) then

          sval = Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
          rho(i) = rho(i) + f(j)*sval
          curr(i) = curr(i)+f(j)*p_part(j)*sval
        end if

        if (abs(r(i)-r_part(j))<=cutoff_avg) then

          avg_rho(i) = avg_rho(i) + f(j)/(dr+dr**3/(12.d0*r(i)**2))*Wn(bsplineorder,(r(i)-r_part(j))/dr)

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)

! Ghost points using the reflection symmetry f(r,p) = f(-r,-p),
! which for scalars integrated over p (rho, avg_rho) is even:
! rho(-r) = rho(r).  The grid is staggered by dr/2 to avoid the
! r=0 singularity (r(i) = (i-0.5)*dr), so the ghost point with
! index (1-k) sits at -r(k) and must mirror the physical point k,
! i.e. rho(1-k) = rho(k) (NOT rho(k-1) = rho(k), which instead
! overwrites the physical point at index (k-1) with the value at
! k, corrupting rho near the origin).
  do i=1,ghost
      rho(1-i) = rho(i)
      avg_rho(1-i) = avg_rho(i)
  end do


  rho = factor*m0*rho/r**2
  avg_rho = factor*m0*avg_rho/r**2

  average_rho = 0.D0

! Integrate with the trapezoidal rule. Second order accurate.

!  do i=1,Nr
    do j=1,Npart
!      if (r1<=r(i) .and. r(i)<=r2) then
      if (r_part(j)>=r1 .and. r_part(j)<= r2) then
      
 !       average_rho = average_rho + 1.0D0/(r2**2-r1**2)*0.5D0*(avg_rho(i)*r(i)**2+avg_rho(i+1)*r(i+1)**2)*dr
        average_rho = average_rho + 1.0D0/(r2**2-r1**2)*f(j)*0.25D0/smallpi

      end if
    end do
!  end do

  filename = 'vlasov_rhomix'
  call save0Ddata(directory,filename,t,average_rho)


end subroutine density



subroutine avg_density

! In order to make the Poisson subroutine more efficient, 
! I will separate the subroutine that calculates 
! the average density from the rest of integrals.

! ******************************************************
! ***   FIND DENSITY AND CURRENT IN PHYSICAL SPACE   ***
! ******************************************************
!
! This subroutine integrates over momentum space
! to find only rho which is defined as:
!                /
! rho  =  1/r**2 | f dp 
!                /

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor,diff,contribution
  real(8) :: cutoff_avg
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)

  smallpi = acos(-1.0d0)


! Zero Angular Momentum
  if (Lfix == 0.0d0) then

     factor = 1.0d0

! Include Angular Momentum
  else

     factor = 2.0*smallpi*Lfix*drc*dpc!0.25D0/smallpi*drc*dpc

  endif  

  avg_rho = 0.D0

! Same cell-list deposit as in density() above, and likewise parallel
! over the grid index alone so each accumulator belongs to a single
! thread. poisson_rk calls this routine on every time step of a
! self-gravitating run, so it is the hottest loop in that case.

  call build_cell_list(cell_start,particle_order)

  cutoff_avg = (dble(bsplineorder) + 1.0d0)*dr
  Wcell = ceiling(cutoff_avg/dr) + 1

! Parallelizing only over "i" (no collapse) means each i is owned by
! exactly one thread for its whole inner loop, so the accumulation
! into avg_rho(i) is race-free without needing !$OMP ATOMIC.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,diff,contribution)

  do i = 1, Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        diff = abs(r(i) - r_part(j))
        if (diff <= cutoff_avg) then
            contribution = f(j) / (r(i)**2*dr+dr**3/12.d0) * Wn(bsplineorder, diff / dr)
            avg_rho(i) = avg_rho(i) + contribution
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)


! Ghost points using the reflection symmetry avg_rho(-r) = avg_rho(r)
! (see the matching note in subroutine density above).
  do i=1,ghost
      avg_rho(1-i) = avg_rho(i)
  end do

  avg_rho = factor*m0*avg_rho

end subroutine avg_density
