
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
  real(8) :: smallpi,factor,average_rho
  real(8) :: cutoff_rho,cutoff_avg,sval,simg,wd,wi,vol
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)
  logical :: images

  character(20) filename ! Name of outupt file.


  smallpi = acos(-1.0d0)


! Zero Angular Momentum
! Unreachable for now: validate (paramfile.f90) stops a run with Lfix = 0 (E1).
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

! Images. With rmin = 0 the distribution obeys f(r,p) = f(-r,-p), so every
! particle has an image at (-r_j,-p_j), and near the origin its weight on
! grid point i, W((r_i+r_j)/dr), is added (to the current with the opposite
! sign, since p changes sign). Without it the part of a particle's weight
! that falls on the ghost points was lost: up to half its mass at r_j = 0,
! 32 % at r_j = dr/4 with bsplineorder = 3 (AUDITORIA_L0_2026-09-21.md, E9).
! It also makes the deposit of a particle that is at r_j < 0 inside a time
! step (main.f90 reflects it only at the end) the same as at -r_j. An image
! only reaches the first points, whose cell range already includes the cells
! near the origin, where build_cell_list files such particles.

  images = (ghost > 0)

! The parallel loop runs over "i" only, never collapsed with "j":
! rho(i), curr(i) and avg_rho(i) accumulate over all particles, so one
! grid point must belong to a single thread for the whole inner loop
! to keep the accumulation race-free without atomics.
!
! The shape function depends on r_part alone, so it is evaluated once
! into "sval" and reused by the density and the current.

! avg_rho is the mass deposited on grid point i over the volume the weight
! covers,
!
!   V_i = Int W_n((r_i-r)/dr) 4 pi r**2 dr = 4 pi dr (r_i**2 + (n+1) dr**2/12),
!
! since W_n (bsplineorder = n) has unit area and second moment (n+1)/12 in
! units of dr**2. With it a uniform density is reproduced exactly at every
! point. The cell volume 4 pi (r_i**2 dr + dr**3/12) used before is right
! for a cell but not for the region the weight spreads the mass over: the
! density came out (r_i**2 + (n+1) dr**2/12)/(r_i**2 + dr**2/12) times the
! true one, +25/50/75 % at the first point for n = 1, 2, 3 whatever dr
! (AUDITORIA_L0_2026-09-21.md, E21). "vol" below is V_i/(4 pi r_i**2).
! This density is only written out; poisson_rk takes its own from
! avg_density.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,sval,simg,wd,wi,vol)

  do i=1,Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)
    vol = dr + dble(bsplineorder+1)*dr**3/(12.d0*r(i)**2)

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        sval = 0.0d0
        simg = 0.0d0
        if (abs(r(i)-r_part(j))<=cutoff_rho) sval = Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
        if (images .and. abs(r(i)+r_part(j))<=cutoff_rho) simg = Sn(bsplineorder,(r(i)+r_part(j))/drc,drc)

        if (sval /= 0.0d0 .or. simg /= 0.0d0) then
          rho(i) = rho(i) + f(j)*(sval + simg)
          curr(i) = curr(i)+f(j)*p_part(j)*(sval - simg)
        end if

        wd = 0.0d0
        wi = 0.0d0
        if (abs(r(i)-r_part(j))<=cutoff_avg) wd = Wn(bsplineorder,(r(i)-r_part(j))/dr)
        if (images .and. abs(r(i)+r_part(j))<=cutoff_avg) wi = Wn(bsplineorder,(r(i)+r_part(j))/dr)

        if (wd /= 0.0d0 .or. wi /= 0.0d0) then
          avg_rho(i) = avg_rho(i) + f(j)/vol*(wd + wi)
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



subroutine avg_density(dens)

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
!
! It is the input of poisson_rk, and is returned in "dens" rather than in
! avg_rho so that the density written out is always the one of density()
! (at t = 0 poisson_rk runs between density() and the output). Here the
! deposited mass is divided by the cell volume 4 pi (r**2 dr + dr**3/12);
! poisson_rk only uses the mass, which it recovers from it, so the
! convention of density() (E21) does not enter the field.

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor,diff,contribution
  real(8) :: cutoff_avg,wd,wi
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)
  logical :: images
  real(8), intent(out) :: dens(1-ghost:Nr)

  smallpi = acos(-1.0d0)


! Zero Angular Momentum
! Unreachable for now: validate (paramfile.f90) stops a run with Lfix = 0 (E1).
  if (Lfix == 0.0d0) then

     factor = 1.0d0

! Include Angular Momentum
  else

     factor = 2.0*smallpi*Lfix*drc*dpc!0.25D0/smallpi*drc*dpc

  endif  

  dens = 0.D0

! Same cell-list deposit as in density() above, and likewise parallel
! over the grid index alone so each accumulator belongs to a single
! thread. poisson_rk calls this routine on every time step of a
! self-gravitating run, so it is the hottest loop in that case.

  call build_cell_list(cell_start,particle_order)

  cutoff_avg = (dble(bsplineorder) + 1.0d0)*dr
  Wcell = ceiling(cutoff_avg/dr) + 1

! Images at -r_j, as in density() above.
  images = (ghost > 0)

! Parallelizing only over "i" (no collapse) means each i is owned by
! exactly one thread for its whole inner loop, so the accumulation
! into dens(i) is race-free without needing !$OMP ATOMIC.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,diff,contribution,wd,wi)

  do i = 1, Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        wd = 0.0d0
        wi = 0.0d0
        diff = abs(r(i) - r_part(j))
        if (diff <= cutoff_avg) wd = Wn(bsplineorder, diff / dr)
        if (images) then
            diff = abs(r(i) + r_part(j))
            if (diff <= cutoff_avg) wi = Wn(bsplineorder, diff / dr)
        end if
        if (wd /= 0.0d0 .or. wi /= 0.0d0) then
            contribution = f(j) / (r(i)**2*dr+dr**3/12.d0) * (wd + wi)
            dens(i) = dens(i) + contribution
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)


! Ghost points using the reflection symmetry dens(-r) = dens(r)
! (see the matching note in subroutine density above).
  do i=1,ghost
      dens(1-i) = dens(i)
  end do

  dens = factor*m0*dens

end subroutine avg_density
