! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays
    use utils
    use distribution

    implicit none

    integer :: i,j,indx
    integer :: unit_ic,ios              ! state="checkpoint"
    real(8) :: smallpi,f_max
    real(8) :: raux,paux
    real(8) :: rand3(3)
    real(8) :: gaussian_fixedL
    real(8) :: halton
    real(8) :: w,x,y,z
    real(8) :: energy

!   Auxiliary variables for a distribution function
!   that depends on action-angle varialbes
    real(8) :: Jr, Qr, s, s1, s2, er1, er2, eta,argaux
!   Auxiliary variables for the (Q3,J3) quadrature grid + Newton-Raphson
!   inversion back to (r,p_r) -- see state "aa_quad" below.
    real(8) :: Jgrid, Qgrid, Egrid, ecc, etaNR, sgrid, rgrid, paux2
    real(8) :: Jminc, Jmaxc, dJc, dQc

    smallpi = acos(-1.0d0)

    call df0_report
    call init_rng

! At fixed L, particles are drawn from an arbitrary distribution function
! f(r,p_r,L) by acceptance-rejection: with fmax the maximum of f, draw
! (x,y,z) uniformly in (rmin,rmax) x (pmin,pmax) x (0,fmax), evaluate
! W = f(x,y,L) and accept the point if z <= W, otherwise draw again.


! Initial data for the density function. Notice that
! we add two copies of the function in order to guarantee 
! that the ! boundary condition f(-r,-p) = f(r,p) is satisfied.

! Find the size of the cell

    drc = (rmaxc-rminc)/dble(Nrc)
    dpc = (pmaxc-pminc)/dble(Npc)

    print *, "(drc,dpc)=",drc,dpc
    if(state.eq."gaussian") then

!     Nodes at the midpoints of the Nrc x Npc cells of the box
!     [rminc,rmaxc] x [pminc,pmaxc], as in the other grid states. This
!     branch used to be named "gaussian1" while read_parameters only accepts
!     "gaussian" (also the default), so the state was never reached: every
!     particle stayed at r = 0 with f = 0 and the run gave NaN. Its nodes were
!     also shifted, r by a whole cell (the last one outside the box) and p to
!     the right edge of each cell (AUDITORIA_L0_2026-09-21.md, E2 and E15).

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+(dble(j)-0.5D0)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = gaussian_fixedL(1.0D0,r0,p0,Lfix,raux,paux,sr,sp)
        end do
      end do
      !$OMP END PARALLEL DO

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc
    
    else if(state .eq."aa") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) SHARED(r_part,p_part,f) &
      !$OMP PRIVATE(j,raux,paux,energy,er1,er2,s,s1,s2,argaux,eta,Qr,Jr)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
!         Midpoint of the cell in p, as in r. The node sat on the right edge,
!         p_j = pminc + j dpc, so the grid ran from pminc+dpc to pmaxc and was
!         not symmetric about p = 0 in a symmetric box
!         (AUDITORIA_L0_2026-09-21.md, E15).
          paux = pminc+(dble(j)-0.5D0)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux     

          energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
          er1 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)-dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          er2 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)+dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
!         On a circular orbit s1 = s2 and the phase is undefined (0/0 gave
!         NaN); guards as in analysish (E13).
          if (s2 > s1) then
            argaux = (s1+s2-2.0*s)/(s2-s1)
          else
            argaux = 0.0d0
          end if

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - sqrt((-2.d0*energy)**3)*sqrt(max(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy,0.0d0))/(-2.d0*energy)*sin(eta)
          Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(Lfix+sqrt(Lfix**2+4.d0))

          f((i-1)*Npc+j) = df0(Qr,Jr)

          if ((f((i-1)*Npc+j) /= f((i-1)*Npc+j) )) then

            f((i-1)*Npc+j) = 0.D0
            r_part((i-1)*Npc+j) = 10000.D0
          end if

        end do
      end do
      !$OMP END PARALLEL DO
      
      f_max = maxval(f)
      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
           if ((r0-0.00)*f_max<=f((i-1)*Npc+j) .and. f((i-1)*Npc+j)<= r0*f_max ) then
              f((i-1)*Npc+j)=0.0D0
           else if (f((i-1)*Npc+j)<= (r0-0.00)*f_max ) then
              r_part((i-1)*Npc+j) = 100000.D0
           end if
!          !print *, Qr,Jr
        end do
      end do 
      !$OMP END PARALLEL DO
      call reduce_arrays(always=.true.)

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      !f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_halton") then
!     Same as "aa", but each grid point is displaced within its own cell
!     by a 2D low-discrepancy (Halton, bases 2 and 3) sequence. This
!     breaks the near-regular lattice in J3 that makes all particles
!     rephase at the same time, the PIC recurrence artifact in h_k(t),
!     while keeping the sample far more uniform than random jitter.

      !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(GUIDED) SHARED(r_part,p_part,f) &
      !$OMP PRIVATE(j,indx,raux,paux,energy,er1,er2,s,s1,s2,argaux,eta,Qr,Jr)
      do i=1,Nrc
        do j=1,Npc
          indx = (i-1)*Npc+j
          raux = rminc+(dble(i)-0.5D0)*drc + (halton(indx,2)-0.5d0)*drc
!         Each point moves within its own cell, centred on the midpoint as in
!         r; it was centred on the right edge, so the cells ran half a cell
!         past pmaxc (E15).
          paux = pminc+(dble(j)-0.5D0)*dpc + (halton(indx,3)-0.5d0)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux

          energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
          er1 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)-dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          er2 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)+dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
!         On a circular orbit s1 = s2 and the phase is undefined (0/0 gave
!         NaN); guards as in analysish (E13).
          if (s2 > s1) then
            argaux = (s1+s2-2.0*s)/(s2-s1)
          else
            argaux = 0.0d0
          end if

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))
          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - sqrt((-2.d0*energy)**3)*sqrt(max(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy,0.0d0))/(-2.d0*energy)*sin(eta)
          Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(Lfix+sqrt(Lfix**2+4.d0))

          f((i-1)*Npc+j) = df0(Qr,Jr)

          if ((f((i-1)*Npc+j) /= f((i-1)*Npc+j) )) then
            f((i-1)*Npc+j) = 0.D0
            r_part((i-1)*Npc+j) = 10000.D0
          end if

        end do
      end do
      !$OMP END PARALLEL DO

      f_max = maxval(f)
      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
           if ((r0-0.00)*f_max<=f((i-1)*Npc+j) .and. f((i-1)*Npc+j)<= r0*f_max ) then
              f((i-1)*Npc+j)=0.0D0
           else if (f((i-1)*Npc+j)<= (r0-0.00)*f_max ) then
              r_part((i-1)*Npc+j) = 100000.D0
           end if
        end do
      end do
      !$OMP END PARALLEL DO
      call reduce_arrays(always=.true.)

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_quad") then
!     Clean tensor-product QUADRATURE rule directly in (Q3,J3): a regular
!     midpoint grid in both variables, NO jitter and NO randomization,
!     inverted back to (r,p_r) by Newton-Raphson on the Kepler-like
!     equation Qr = eta - ecc*sin(eta).
!
!     Without self-interaction J3 is exactly conserved and the phase is
!     exactly Q(t)=Q(0)+omega(J)t, so h_k(t) is not a statistical
!     sampling problem at all -- it is a QUADRATURE of a smooth but
!     increasingly oscillatory integral, with
!        n_osc(t) = k*Delta_omega*t/(2*pi)
!     oscillations across the support in J. For a smooth integrand a
!     designed quadrature beats Monte Carlo by orders of magnitude: the
!     error collapses as soon as the J resolution crosses Nyquist,
!        Nrc  >~ 4*k*Delta_omega*t_max/(2*pi),
!     and then degrades catastrophically (aliasing = the Birdsall &
!     Langdon recurrence, seen from the other side) beyond that. So
!     pick Nrc for the t you intend to run and do not trust results
!     past t ~ T_rec/2. In Q the rule is the periodic trapezoid, which
!     converges spectrally -- Npc ~ 40 is already enough, so nearly all
!     particles should go into resolving J, not Q.
!
!     Jitter of any kind must be avoided here: it destroys the ordered
!     node placement the quadrature error estimate relies on, leaving
!     only the 1/sqrt(N) Monte Carlo rate.
!
!     Normalization: every quadrature weight dQc*dJc is the same, so the
!     constant is absorbed by the mass normalization below and f holds
!     the raw value of the distribution function, as in "aa". The
!     normalization uses drc*dpc because density, energy and analysish
!     all multiply f by drc*dpc, so the factor cancels.

      call df0_Jrange(Jminc,Jmaxc)
      dJc = (Jmaxc-Jminc)/dble(Nrc)
      dQc = 2.0d0*smallpi/dble(Npc)

      !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(GUIDED) &
      !$OMP PRIVATE(j,Jgrid,Qgrid,Egrid,ecc,etaNR,sgrid,rgrid,paux2,raux,paux,er1,er2,s1,s2,argaux,Jr,Qr) &
      !$OMP SHARED(r_part,p_part,f)
      do i=1,Nrc          ! index over J3 (resolves the oscillation)
        do j=1,Npc        ! index over Q3 (periodic trapezoid)

          Jgrid = Jminc + (dble(i)-0.5D0)*dJc
          Qgrid = (dble(j)-0.5D0)*dQc

!         Invert Jr(E) for E at fixed L=Lfix:
          Egrid = -1.d0/(2.d0*(Jgrid+0.5d0*(Lfix+sqrt(Lfix**2+4.d0)))**2)

!         Radicands guarded as in invert_QJ_to_rp (E13).
          er1 = dsqrt(max((1.d0+Egrid*(2.d0+Lfix**2)-dsqrt(max(1.d0+2.d0*Egrid*(2.d0+2.d0*Egrid+Lfix**2),0.0d0))) &
                /(2.d0*Egrid**2),0.0d0))
          er2 = dsqrt(max((1.d0+Egrid*(2.d0+Lfix**2)+dsqrt(max(1.d0+2.d0*Egrid*(2.d0+2.d0*Egrid+Lfix**2),0.0d0))) &
                /(2.d0*Egrid**2),0.0d0))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)

          ecc = sqrt((-2.d0*Egrid)**3)*sqrt(max(-Lfix**2-2.d0*Egrid-2.d0-0.5D0/Egrid,0.0d0))/(-2.d0*Egrid)

!         Solve Qgrid = etaNR - ecc*sin(etaNR): safeguarded Newton-Raphson,
!         kepler_eta in utils.f90 (E12).
          etaNR = kepler_eta(Qgrid,ecc,1.0d-13)

          argaux = cos(etaNR)
          sgrid = (s1+s2-argaux*(s2-s1))/2.0d0
          rgrid = sqrt(max((sgrid-1.d0)**2-1.d0,0.0d0))

          paux2 = 2.d0*(Egrid + 1.d0/(1.d0+dsqrt(1.d0+rgrid**2)) - 0.5d0*Lfix**2/max(rgrid**2,1.0d-12))
          paux2 = sqrt(max(paux2,0.0d0))
          if (mod(etaNR,2.0d0*smallpi) > smallpi) then
            paux = -paux2
          else
            paux = paux2
          end if
          raux = rgrid

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux

          Jr = Jgrid
          Qr = Qgrid
          f((i-1)*Npc+j) = df0(Qr,Jr)

          if ((f((i-1)*Npc+j) /= f((i-1)*Npc+j)) .or. (raux /= raux)) then
            f((i-1)*Npc+j) = 0.D0
            r_part((i-1)*Npc+j) = 10000.D0
          end if

        end do
      end do
      !$OMP END PARALLEL DO

!     NOTE: deliberately NO r0 cutoff and NO reduce_arrays here. Dropping
!     the low-f nodes would truncate the quadrature rule, and at late t
!     the integral survives only through near-total cancellation, so even
!     a 1%-level truncation of the tails can dominate the answer.

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_random") then

      Npart = Nrc*Npc
      f_max = df0_max()

      i = 1
      do while(i<= Npart)
!      do i=1,Npart
!        accepted = .false.
!        do while (.not. accepted)
! The subroutine random_number creates a pseudo-random number in the interval (0,1],
! we would like to include the zero, so here I just create a triad of three uniform 
! random numbers in the interval [0,1]

        call random_number(rand3)

        x = 1.d0-rand3(1)
        y = 1.d0-rand3(2)
        z = 1.d0-rand3(3)
  
! Random position and momentum of the particle.

        raux = (rmaxc-rminc)*x + rminc
        paux = (pmaxc-pminc)*y + pminc
        z = f_max*z 

        energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
        if (energy < 0.0D0) then

          er1 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)-dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          er2 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)+dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))) &
                /(2.d0*energy**2),0.0d0))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          !eta = dacos(sign(min(abs(2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)),1.0),2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)))
!         On a circular orbit s1 = s2 and the phase is undefined (0/0 gave
!         NaN); guards as in analysish (E13).
          if (s2 > s1) then
            argaux = (s1+s2-2.0*s)/(s2-s1)
          else
            argaux = 0.0d0
          end if

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - dsqrt((-2.d0*energy)**3)*dsqrt(max(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy,0.0d0))/(-2.d0*energy)*dsin(eta)
          Jr = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))
           w = df0(Qr,Jr)

          if (z <= w ) then

            r_part(i) = raux
            p_part(i) = paux
!           Equal weights: rejection sampling already places particles
!           with number density proportional to F, so the shape of the
!           distribution comes from where they land. Weighting them by F
!           again would represent F^2.
            f(i)      = 1.0d0
            !accepted = .true.
            i = i+1
          end if
        end if
        !end do
      end do




!      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
!      !$OMP END PARALLEL DO
!     Normalize to the requested total mass a0. The drc*dpc factor stays
!     even though these particles are not on a grid: density, energy and
!     analysish always multiply f by drc*dpc, so dividing by it here lets
!     it cancel and the sums recover the Monte Carlo estimate.
      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc



    else if(state.eq."checkpoint") then
!     Particles read from a file, one line "r p_r f" per particle, exactly
!     Npart = Nrc*Npc lines. This is how an initial state that the code
!     cannot build itself enters the run -- for instance a self-consistent
!     equilibrium plus a perturbation, placed on a quadrature grid in the
!     action-angle variables of the total potential
!     (paper_runs/scripts/equilibrio.py).
!
!     f holds raw values of the distribution function on a grid of equal
!     cells, as in "aa_quad", so the same normalization to the mass a0
!     applies: density, energy and analysish multiply f by drc*dpc, which
!     cancels here.

      open(newunit=unit_ic,file=trim(CheckPointfile),status='old',action='read',iostat=ios)
      if (ios /= 0) then
         print *
         print *, 'state="checkpoint": cannot open ',trim(CheckPointfile)
         print *, 'Aborting ...'
         print *
         stop 1
      end if
      do i=1,Npart
         read(unit_ic,*,iostat=ios) r_part(i),p_part(i),f(i)
         if (ios /= 0) then
            print *
            print *, 'state="checkpoint": ',trim(CheckPointfile),' has fewer than Npart=Nrc*Npc =',Npart,' lines'
            print *, 'Aborting ...'
            print *
            stop 1
         end if
      end do
      read(unit_ic,*,iostat=ios) raux
      if (ios == 0) then
         print *
         print *, 'state="checkpoint": ',trim(CheckPointfile),' has more than Npart=Nrc*Npc =',Npart,' lines'
         print *, 'Aborting ...'
         print *
         stop 1
      end if
      close(unit_ic)

      f = a0/(drc*dpc*8.0d0*smallpi**2*Lfix*sum(f))*f
      print *, "Read",Npart," particles from ",trim(CheckPointfile)
      print *, "Initial total mass=",sum(f)*8.0d0*smallpi**2*Lfix*drc*dpc

    else

!      read_parameters only lets through the states it lists; this catches a
!      state added there without its branch here, which is how "gaussian"
!      went silently to NaN.
       print *
       print *, 'Unknown initial state: ',trim(state)
       print *, 'Aborting ...'
       print *
       stop 1

    endif

  end subroutine initial_data

  function gaussian_fixedL(a0,r0,p0,L,x,y,sr,sp)

  implicit none

  real(8) a0
  real(8) gaussian_fixedL
  real(8) L,r0,p0,x,y,sr,sp 
  real(8) smallpi

  smallpi = acos(-1.0d0)

  if (L == 0) then
    gaussian_fixedL = dble(a0)/(smallpi*sr*sp)*&
               (dexp(-(x-r0)**2/sr**2)*dexp(-(y-p0)**2/sp**2)+ &
                dexp(-(x+r0)**2/sr**2)*dexp(-(y+p0)**2/sp**2))
  else
    gaussian_fixedL = dble(a0)/(2.0d0*smallpi*L*smallpi*sr*sp)*&
               (dexp(-(x-r0)**2/sr**2)*dexp(-(y-p0)**2/sp**2)+ &
                dexp(-(x+r0)**2/sr**2)*dexp(-(y+p0)**2/sp**2))
  end if

  end function gaussian_fixedL

  function halton(idx,base) result(h)
!   Radical-inverse Halton sequence value for index idx>=1 in the given
!   base (use coprime bases, e.g. 2 and 3, for a genuine 2D low-discrepancy
!   sequence). Result in [0,1).

  implicit none

  integer :: idx,base
  real(8) :: h
  real(8) :: f
  integer :: n

  h = 0.0d0
  f = 1.0d0/dble(base)
  n = idx
  do while (n > 0)
    h = h + f*dble(mod(n,base))
    n = n/base
    f = f/dble(base)
  end do

  end function halton


