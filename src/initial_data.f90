! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays
    use utils

    implicit none

    logical :: accepted
    integer :: i,j,indx
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
    real(8) :: Jgrid, Qgrid, Egrid, ecc, etaNR, gNR, gpNR, sgrid, rgrid, paux2
    real(8) :: Jminc, Jmaxc, dJc, dQc
    integer :: iterNR

    smallpi = acos(-1.0d0)

! 
! For a fixed value of L, We generate particles for an 
! arbitrary distribution function f(r,pr,L) via an acceptance-rejection method.
! Let fmax the maximum value of f. We generate arbitrary (x,y,z) numbers 
! in the range of (rmin,rmax), (pmin,pmax), (0,fmax) respectively. 
! Then evaluate W=f(x,y,L), if z<=W, accept the point, 
! otherwise, repeat until the condition in fulfilled.


! Initial data for the density function. Notice that
! we add two copies of the function in order to guarantee 
! that the ! boundary condition f(-r,-p) = f(r,p) is satisfied.

! Find the size of the cell

    drc = (rmaxc-rminc)/dble(Nrc)
    dpc = (pmaxc-pminc)/dble(Npc)

    print *, "(drc,dpc)=",drc,dpc
    if(state.eq."gaussian1") then

!      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j)
!      do i=1,Nrc
!        if (rminc == 0.d0) then
!          raux = (dble(i)+0.5d0)*drc
!          raux = (dble(i)-0.5d0)*dr

!        else
!          raux = rminc+(dble(i)+0.5)*drc
!          raux = rminc+(dble(i)-0.5D0)*dr
!        end if

!        do j=1,Npc
!          paux = pminc+dble(j)*dpc

!          r_part((i-1)*Npc+j) = raux
!          p_part((i-1)*Npc+j) = paux          
!          f((i-1)*Npc+j)      = gaussian_fixedL(a0,r0,p0,Lfix,raux,paux,sr,sp)
!        end do
!      end do
!      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)+0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = gaussian_fixedL(1.0D0,r0,p0,Lfix,raux,paux,sr,sp)
        end do
      end do
      !$OMP END PARALLEL DO

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc
    
    else if(state.eq."gaussian2") then

! For a gaussian distribution, the above choice 
! makes the normalization analytic, so that the normalization factor 
! "N0" now corresponds directly to the total initial number of particles.

       if (Lfix == 0.0d0) then ! Zero Angular Momentum
          f_max = Npart/(smallpi*sr*sp)
       else ! Include Angular Momentum
          f_max = Npart/(2.0d0*smallpi*Lfix*smallpi*sr*sp)
       endif

! Populate the phase space with particles.

      i = 1

      do while (i<=Npart)

! The subroutine random_number creates a pseudo-random number in the interval (0,1],
! we would like to include the zero, so here I just create a triad of three uniform 
! random numbers in the interval [0,1]

        call random_number(rand3)

        x = 1.d0-rand3(1)
        y = 1.d0-rand3(2)
        z = 1.d0-rand3(3)
  
! Random position and momentum of the particle.

        x = (rmaxc-rminc)*x + rminc
        y = (pmaxc-pminc)*y + pminc
        z = f_max*z 
        w = gaussian_fixedL(a0,r0,p0,Lfix,x,y,sr,sp)
        if (z <= w ) then
        
          r_part(i) = x
          p_part(i) = y
          f(i)      = w
          i = i+1          
        end if



      end do

    else if (state.eq."Plummer") then 

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux,energy)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i))*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux

          energy = -1.0/(1.0D0+sqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2 + eps*eps) + 0.5D0*paux**2



          if (energy<0.d0) then
            f((i-1)*Npc+j) = (-energy)**3.5
          else
            f((i-1)*Npc+j) = 0.0D0
          end if
        end do
      end do
      !$OMP END PARALLEL DO
    
    else if(state == "compact") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = a0/(32.D0*smallpi**2*sr*sp)*(1.D0+dcos(smallpi/sr*(raux-r0))) * &
                                                              (1.D0+dcos(smallpi/sp*(paux-p0))) 
        end do
      end do
      !$OMP END PARALLEL DO

      f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)

    else if(state == "compact2") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = 2.0D0*a0/(9.D0*smallpi**2*sr*sp)*(dcos(0.5D0*smallpi/sr*(raux-r0)))**4* &
                                                            (dcos(0.5D0*smallpi/sp*(paux-p0)))**4 
        end do
      end do
      !$OMP END PARALLEL DO

      f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)

    else if(state .eq."aa") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux) SHARED(r_part,p_part,f)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux     

          energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
          er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          argaux = (s1+s2-2.0*s)/(s2-s1)

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - sqrt((-2.d0*energy)**3)*sqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*sin(eta)
          Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(Lfix+sqrt(Lfix**2+4.d0))

          f((i-1)*Npc+j) = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

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
      call reduce_arrays

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      !f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_halton") then
!     Same as "aa" but the regular (r,p) grid is jittered with a genuine
!     2D low-discrepancy (Halton, bases 2 and 3) sequence instead of the
!     original regular grid, to break the near-regular J3 lattice that
!     causes the PIC recurrence artifact in h_k(t). See BUGS_TODO.md.

      !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(GUIDED) PRIVATE(j,indx,raux,paux) SHARED(r_part,p_part,f)
      do i=1,Nrc
        do j=1,Npc
          indx = (i-1)*Npc+j
          raux = rminc+(dble(i)-0.5D0)*drc + (halton(indx,2)-0.5d0)*drc
          paux = pminc+dble(j)*dpc         + (halton(indx,3)-0.5d0)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux

          energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
          er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          argaux = (s1+s2-2.0*s)/(s2-s1)

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))
          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - sqrt((-2.d0*energy)**3)*sqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*sin(eta)
          Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(Lfix+sqrt(Lfix**2+4.d0))

          f((i-1)*Npc+j) = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

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
      call reduce_arrays

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_quad") then
!     Clean tensor-product QUADRATURE rule directly in (Q3,J3): a regular
!     midpoint grid in both variables, NO jitter and NO randomization,
!     inverted back to (r,p_r) by Newton-Raphson on the Kepler-like
!     equation Qr = eta - ecc*sin(eta).
!
!     Why this (and why it is NOT the same as the reverted "aa_qj"):
!     without self-interaction J3 is exactly conserved and the phase is
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
!     "aa_qj" (reverted) failed for unrelated reasons: a NaN-propagation
!     bug that froze the J range at its +-1e30 sentinels, and then, once
!     that was fixed, Weyl jitter on J plus random Q -- which is exactly
!     what destroys the quadrature property this state relies on. The
!     earlier claim that a uniform Q grid "must" cancel the Fourier sum
!     was wrong: that argument applies to an UNWEIGHTED sum of roots of
!     unity, whereas here the sum is weighted by F(Q), which converges
!     spectrally to the true Fourier coefficient.
!
!     Normalization: all quadrature weights dQc*dJc are EQUAL here, so
!     the constant is absorbed by the mass normalization below and f can
!     simply hold the raw DF value, exactly as in "aa". Note this must
!     use drc*dpc (not dJc*dQc): every consumer of f() -- density.f90,
!     energy.f90, analysish.f90 -- multiplies by drc*dpc unconditionally,
!     so using the same factor here makes it cancel. ("aa_qj" normalized
!     with dJc*dQc instead, which left a spurious drc*dpc/(dJc*dQc)
!     factor in h_k -- that is exactly the unexplained 2.87x bias in its
!     h_0, see BUGS_TODO.md.)

      Jminc = 1.0d-4*sr
      Jmaxc = 6.0d0*sr
      dJc = (Jmaxc-Jminc)/dble(Nrc)
      dQc = 2.0d0*smallpi/dble(Npc)

      !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(GUIDED) &
      !$OMP PRIVATE(j,Jgrid,Qgrid,Egrid,ecc,etaNR,gNR,gpNR,iterNR,sgrid,rgrid,paux2,raux,paux,er1,er2,s1,s2,argaux,Jr,Qr) &
      !$OMP SHARED(r_part,p_part,f)
      do i=1,Nrc          ! index over J3 (resolves the oscillation)
        do j=1,Npc        ! index over Q3 (periodic trapezoid)

          Jgrid = Jminc + (dble(i)-0.5D0)*dJc
          Qgrid = (dble(j)-0.5D0)*dQc

!         Invert Jr(E) for E at fixed L=Lfix:
          Egrid = -1.d0/(2.d0*(Jgrid+0.5d0*(Lfix+sqrt(Lfix**2+4.d0)))**2)

          er1 = dsqrt((1.d0+Egrid*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*Egrid*(2.d0+2.d0*Egrid+Lfix**2)))/(2.d0*Egrid**2))
          er2 = dsqrt((1.d0+Egrid*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*Egrid*(2.d0+2.d0*Egrid+Lfix**2)))/(2.d0*Egrid**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)

          ecc = sqrt((-2.d0*Egrid)**3)*sqrt(-Lfix**2-2.d0*Egrid-2.d0-0.5D0/Egrid)/(-2.d0*Egrid)

!         Newton-Raphson solve of Qgrid = etaNR - ecc*sin(etaNR).
          etaNR = Qgrid
          do iterNR=1,50
            gNR  = etaNR - ecc*sin(etaNR) - Qgrid
            gpNR = 1.d0 - ecc*cos(etaNR)
            etaNR = etaNR - gNR/gpNR
            if (abs(gNR) < 1.0d-13) exit
          end do

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
          f((i-1)*Npc+j) = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

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
      f_max = dexp(-dsin(0.5d0*0.0)**2/sp**2)*dexp(-sr**2/sr**2)*sr**2

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

          er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          !eta = dacos(sign(min(abs(2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)),1.0),2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)))
          argaux = (s1+s2-2.0*s)/(s2-s1)

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - dsqrt((-2.d0*energy)**3)*dsqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*dsin(eta)
          Jr = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))
           w = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

          if (z <= w ) then

            r_part(i) = raux
            p_part(i) = paux
!           BUG FIX: this used to set f(i)=w, double-counting F -- the
!           accepted particles from rejection sampling are ALREADY
!           distributed with density proportional to F(Qr,Jr) (that is
!           what rejection sampling means), so weighting them by F again
!           on top of that biases the reconstructed distribution towards
!           already-dense regions (effectively ~F^2, renormalized) instead
!           of representing F itself. Equal-weight macroparticles is the
!           correct MC representation; the density *shape* comes entirely
!           from where particles land, not from their individual weight.
            f(i)      = 1.0d0
            !accepted = .true.
            i = i+1
          end if
        end if
        !end do
      end do




!      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
!      !$OMP END PARALLEL DO
!     drc*dpc IS needed here even though these particles are not on a
!     regular grid: every consumer of f() (density.f90, energy.f90,
!     analysish.f90) treats it as "F evaluated at that point" and
!     multiplies by the drc*dpc quadrature weight itself, unconditionally,
!     regardless of how the particle was placed -- so f() must be
!     pre-divided by drc*dpc here for that multiplication to reconstruct
!     the correct MC mass estimate downstream. (An earlier version of
!     this fix dropped drc*dpc, reasoning these are already-random
!     samples that need no cell-size weight on their own -- true in
!     isolation, but wrong given analysish.f90 always multiplies by it;
!     verified via h_0, which per phase-mixing theory must be constant
!     in time and equal to the exact value ~1.4998e-8 from
!     paper_runs/notebooks/hk_exact.ipynb -- dropping drc*dpc gave
!     ~1.56e-10 (off by ~drc*dpc), restoring it gives the right order.)
      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc



    else if(state.eq."checkpoint") then !NO IMPLEMENTED

!       open(101,file=CheckPointFile)

!       do i=0,Nr
!          do j=0,Np
!             read(101,*) aux1, aux2, f(i,j)
!          end do
!       end do

!       close(101)

    else if(state.eq."other3") then !NO IMPLEMENTED
       f = 0.0d0       
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


