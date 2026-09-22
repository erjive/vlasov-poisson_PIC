  subroutine poisson_rk

! Erik: This subroutine was originally written by Miguel Alcubierre
! I adapted it for the vlasov_PIC code. 
! WARNING: This subroutine does not execute in parallel
! *******************
! ***   POISSON   ***
! *******************

! This subroutine solves the Poisson equation
! for the self-gravitating case:
!
! Laplacian (pot)  =  4 pi rho
!
! with pot the gravitational field and rho the mass
! density.  Once we have pot, the gravitational
! force is calculated as:
!
! force = - grad (pot)
!
! In spherical symmetry the above equations reduce to:
!
!  2
! d pot  +  (2/r) d pot  =  4 pi rho
!  r               r
!
! force  =  - d pot /dr

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j

  real(8) ra,rb,slope,c0,c3,c4
  real(8) rho0,pi
  real(8), allocatable :: rt(:)
  real(8) cutoff_interp,wgt
  integer :: Wgrid,jc,m
  real(8) :: rj,pj,fj,rm,sgn

! *******************
! ***   NUMBERS   ***
! *******************

  pi = acos(-1.d0)


! ****************************************
! ***   FIND GRAVITATIONAL POTENTIAL   ***
! ****************************************

! Initialize arrays to zero.

  pot = 0.0d0
  dev_pot = 0.0d0

! The equation is integrated in the enclosed mass rather than in dPhi/dr,
!
!   dM/dr = 4 pi r**2 rho,      dPhi/dr = M/r**2,
!
! which removes the 2/r term and with it the only source of stiffness at the
! origin. The second order Runge-Kutta this file is named after integrated
! dPhi/dr directly and lost the mass of the particles near the origin: a
! particle within dr/2 of it produced no field at all, one at 1.5 dr 1.21
! times its mass (AUDITORIA_L0_2026-09-21.md, E8). Same scheme as
! VlasovPoisson_PIC_sp (5c0e07f).

! First calculate the density
  call avg_density

! Between two grid points the density is the straight line through its two
! values. avg_density divides the mass deposited on point k, m_k, by the cell
! volume 4 pi (r_k**2 dr + dr**3/12) (the convention of the manuscript). The
! straight line through values rt_k holds 4 pi Sum_k rt_k dr (r_k**2 + dr**2/6),
! because the linear hat of width dr has second moment dr**2/6, so with
!
!   rt_k = m_k / (4 pi dr (r_k**2 + dr**2/6))
!
! the field sees exactly the deposited mass, Sum_k m_k, for any bsplineorder
! and at any radius, the first cell included (see below). rt only lives here:
! the density written out keeps the cell volume.

  allocate(rt(1:Nr))
  do i=1,Nr
     rt(i) = avg_rho(i)*(r(i)**2 + dr**2/12.d0)/(r(i)**2 + dr**2/6.d0)
  end do

! Between the origin and r(1) the density is even, so the ghost point carries
! the same value as r(1) and rho is constant there:
!
!   M(r1) = (4/3) pi rho(1) r1**3,   Phi(r1) - Phi(0) = (2/3) pi rho(1) r1**2.
!
! Phi(0) = 0 is arbitrary; the additive constant is fixed at the end. Note 2
! of the original code still holds: with rmin = 0 there are two ghost points
! to the left of the origin and the first point with positive r is i=1.

  rho0 = rt(1)

! Between r(i-1) and r(i), with rho(x) = rt(i-1) + slope (x - ra), M is a
! quartic with no linear or quadratic term,
!
!   M(x) = c0 + c3 x**3 + c4 x**4,   c3 = 4 pi (rt(i-1) - slope ra)/3,
!                                    c4 = pi slope,
!
! and both integrals are done in closed form, with no quadrature error:
!
!   Int M/x**2 dx = -c0/x + c3 x**2/2 + c4 x**3/3.
!
! dev_pot carries M while the integration runs, and is divided by r**2 at the
! end. It is a local role: outside this loop dev_pot is always dPhi/dr.

  dev_pot(1) = 4.d0/3.d0*pi*rho0*r(1)**3
  pot(1)     = 2.d0/3.d0*pi*rho0*r(1)**2

  do i=2,Nr

     ra = r(i-1)
     rb = r(i)

     slope = (rt(i) - rt(i-1))/dr

     c3 = 4.d0*pi*(rt(i-1) - slope*ra)/3.d0
     c4 = pi*slope
     c0 = dev_pot(i-1) - c3*ra**3 - c4*ra**4

     dev_pot(i) = c0 + c3*rb**3 + c4*rb**4

     pot(i) = pot(i-1) + (- c0/rb + c3*rb**2/2.d0 + c4*rb**3/3.d0) &
                       - (- c0/ra + c3*ra**2/2.d0 + c4*ra**3/3.d0)

  end do

  deallocate(rt)

! From the enclosed mass to dPhi/dr.

  dev_pot(1:Nr) = dev_pot(1:Nr)/r(1:Nr)**2

! Ghost points using symmetries.

  pot(-1) = pot(2)
  pot( 0) = pot(1)

  dev_pot(-1) = - dev_pot(2)
  dev_pot( 0) = - dev_pot(1)

! Calculate the force.

  force = - dev_pot

! We now need to substract a constant to the solution
! to make sure that far away the potential behaves as
! pot ~ 1/r.  This condition implies that we should
! have pot + r*dev_pot = 0 far away.  But we won't since
! we arbitrarily fixed pot(r=0)=0 above. So now we
! find the value of pot + r*dev_pot at the outer boundary
! and just substract it from the whole solution.

  pot = pot - (pot(Nr) - force(Nr)*r(Nr))

! *******************************
! ***   ADD ANGULAR MOMENTUM  ***
! *******************************
!  if(Lfix /= 0.0d0) then
!     pot   = pot + 0.5d0*Lfix*Lfix/(r*r + eps*eps)
!     force = force + Lfix*Lfix*r/(r*r + eps*eps)**2
!  end if

! So far we have solved for the potential and the force 
! felt on the mesh. In order to calculate the force that 
! particles felt, we need to interpolate the potential 
! and force on each of them.

  pot_part   = 0.0D0
  force_part = 0.0D0

! Interpolate the grid potential and force back to the particles, with the
! same weight W_n used in the deposit, over every grid point j within its
! support. The grid is uniform and staggered, r_j = (j - 1/2) dr, so the few
! points within the support are found from the particle's radius instead of
! scanning all Nr points. The field is known beyond the stored points:
!
!   j <= 0   the mirror of point 1-j (r_j = -r_(1-j)): the potential is even
!            and the force odd, the symmetry f(r,p) = f(-r,-p) at the origin;
!   beyond   a point past r(Nr), stored or mirrored, takes the exterior
!            solution Phi = Phi(r_Nr) r_Nr/r and F = F(r_Nr) (r_Nr/r)**2 of
!            the mass on the grid.
!
! So the weights of every particle add up to one at any radius. Only points
! 1..Nr were used before: a particle closer to the origin than 1.5 dr felt a
! force 30 to 50 times the exact one, and one beyond r(Nr) + dr no
! self-gravity at all (AUDITORIA_L0_2026-09-21.md, E10). For a particle
! whose support lies within 1..Nr nothing changes. Same scheme as
! VlasovPoisson_PIC_sp (a81f4cf).
!
! The parallel loop runs over particles only, never collapsed with the
! grid index: pot_part(i) and force_part(i) accumulate over j, so each
! particle must stay with one thread to avoid a race. The interpolation
! weight depends only on the particle-grid distance, so it is formed once
! into "wgt" and shared by the potential and the force.

  Wgrid = (bsplineorder+2)/2
  cutoff_interp = 0.5d0*dble(bsplineorder+1)*dr

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,jc,wgt,rj,pj,fj,m,rm,sgn)

  do i=1,Npart

    jc  = nint((r_part(i)-r(1))/dr) + 1

    do j=jc-Wgrid,jc+Wgrid

!     The same expression as construct_grid, so r_j = r(j) for 1 <= j <= Nr.
      rj = (dble(j)-0.5d0)*dr
      if (abs(r_part(i)-rj) >= cutoff_interp) cycle

      if (j <= 0) then
        m = 1-j
        sgn = -1.0d0
      else
        m = j
        sgn = 1.0d0
      end if
      if (m > Nr) then
        rm = (dble(m)-0.5d0)*dr
        pj = pot(Nr)*r(Nr)/rm
        fj = sgn*force(Nr)*(r(Nr)/rm)**2
      else
        pj = pot(m)
        fj = sgn*force(m)
      end if

      wgt = Wn(bsplineorder,(r_part(i)-rj)/dr)

      pot_part(i)   = pot_part(i) + pj*wgt

      force_part(i) = force_part(i) + fj*wgt

    end do
  end do
  !$OMP END PARALLEL DO

!  filename = 'vlasov_potpart'
!  call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,pot_part)

!     if (mod(l,spatial_output).eq.0) then

!  filename = 'vlasov_potpoisson'
!  call save1Ddata(directory,filename,Nr,t,r,pot)

!  filename = 'vlasov_potpoisson_r0'
!  call save0Ddata(directory,filename,t,pot(1))

!     end if



! ***************
! ***   END   ***
! ***************

  end subroutine poisson_rk


