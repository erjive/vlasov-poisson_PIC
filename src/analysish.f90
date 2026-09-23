! Projects the particle distribution onto two Gaussian test functions
! Phi_1, Phi_2 and records the first five angular Fourier modes of each.
!
! For a distribution sharply peaked at a single angular momentum,
! F(r,p_r,L) = F_0(r,p_r) delta(L-L0), the projection reduces to
!
!   h_k = 8 pi^2 L0  Int F_0(Q,J) Phi(Q,J) e^{-ikQ} dQ dJ,
!
! with (Q,J) the radial angle-action pair. The 8 pi^2 L0 factor is the
! measure of the angular degrees of freedom that the delta collapses.
!
! Particles carry (r,p_r), so each one is first mapped to (Q,J); the
! integral then becomes a weighted sum over particles, dQ dJ = dr dp_r
! because the transformation is canonical.
!
! Output: hk1.tl / hk2.tl hold |h_k|, and hk1_complex.tl / hk2_complex.tl
! hold Re and Im separately. The phase matters when averaging several
! runs: discretisation noise has a random phase between realisations and
! cancels in a complex average, whereas averaging magnitudes never lets
! it cancel.

  subroutine analysish

    use parameters
    use arrays
!$  use omp_lib

    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8) :: energy,s, s1, s2, er1, er2,argaux,disc
    logical,dimension(1:Npart) :: bound     !E < 0: the particle has action-angle variables
    complex(8) :: ii                          !imaginary unit
    complex(8) :: expv                        !exp(-ii*Qr(j)) for the current particle
    complex(8),dimension(0:4) :: hk1,hk2           !h_k mode
    complex(8),dimension(0:4) :: lh1,lh2           !Partial sums of one thread
    complex(8),allocatable :: parth(:,:,:)
    integer  :: nth,tid
    real(8),dimension(0:4) :: abs_hk1,abs_hk2       !Magnitude h_k mode
    integer  :: i,j,k                         !Counters
    integer  :: mode = 4                          !Number of modes

    real(8) :: eta
    real(8) :: smallpi

    integer, parameter :: nquad = 512         !Simpson intervals for the phik integral (must be even)
                                               !kept identical to the original per-call phik() precision
    real(8) :: quadQ(0:nquad),quadW(0:nquad)  !Quadrature nodes/weights on [0,pi], built once
    real(8) :: hstep,quadnorm
    real(8) :: aq1(0:nquad),aq2(0:nquad)      !Q-only factor of the integrand, per quadrature node
    real(8) :: cq1(0:4),cq2(0:4)              !particle-INDEPENDENT quadrature constant, per mode
    real(8) :: bj1,bj2                        !J-only factor, for the current particle

    character(20) filestatus


    ! Constants
    smallpi =  acos(-1.0d0)


! Unbound particles (E >= 0) have no action-angle variables: the formulas
! give NaN, and a single one turned every h_k into NaN. They are left out,
! which is the continuous extension of the test functions (J -> infinity as
! E -> 0-, and B(J) -> 0). Every radicand is clamped at zero: all of them
! vanish on a circular orbit, where rounding alone can make them negative.
! Where a radicand is positive the clamp returns it unchanged, bit for bit
! (AUDITORIA_L0_2026-09-21.md, E7).

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,energy,s,s1,s2,er1,er2,argaux,eta,disc) SHARED(Qr,Jr,bound)
    do i = 1,Npart
      energy = -1.0/(1.0D0+dsqrt(1.0D0+r_part(i)**2)) + 0.5d0*Lfix**2/(r_part(i)**2) + 0.5D0*p_part(i)**2
      bound(i) = (energy < 0.0d0)
      if (.not. bound(i)) then
        Qr(i) = 0.0d0
        Jr(i) = 0.0d0
        cycle
      end if
      disc = dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2),0.0d0))
      er1 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)-disc)/(2.d0*energy**2),0.0d0))
      er2 = dsqrt(max((1.d0+energy*(2.d0+Lfix**2)+disc)/(2.d0*energy**2),0.0d0))
      s1 = 1.d0 + dsqrt(1.d0+er1**2)
      s2 = 1.d0 + dsqrt(1.d0+er2**2)
      s  = 1.d0 + dsqrt(1.d0+r_part(i)**2)
      if (s2 > s1) then
        argaux = (s1+s2-2.0*s)/(s2-s1)
      else
        argaux = 0.0d0
      end if
      Jr(i) = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))

  !  do i=1,Npart

          if (p_part(i)>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0D0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0D0),argaux))+smallpi
          end if

          Qr(i) = eta - dsqrt((-2.d0*energy)**3)*sqrt(max(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy,0.0d0))/(-2.d0*energy)*dsin(eta)

    end do
    !$OMP END PARALLEL DO

    ii = (0.d0,1.d0)

! Simpson nodes and weights on [0,pi] with nquad (even) sub-intervals.
! The test functions are even and 2pi-periodic in Q, so the integral over
! the full period is twice the one over [0,pi] and only cosines survive.

    hstep = smallpi/dble(nquad)
    do k=0,nquad
      quadQ(k) = dble(k)*hstep
    end do
    quadW(0)     = 1.d0
    quadW(nquad) = 1.d0
    do k=1,nquad-1,2
      quadW(k) = 4.d0
    end do
    do k=2,nquad-2,2
      quadW(k) = 2.d0
    end do
    quadnorm = hstep/(3.d0*smallpi)

! The test function is a product of a Q-part and a J-part,
!
!   Phi(Q,J) = exp(-sin(Q/2)^2/sq^2) * exp(-(J-J0)^2/sj^2) * J^2
!            =       A(Q)            *          B(J),
!
! so its Fourier coefficient separates,
!
!   Phi_k(J) = B(J) * [ quadnorm * sum_n w_n A(Q_n) cos(k Q_n) ]
!            = B(J) * C(k).
!
! C(k) involves no particle data, so the whole Q-quadrature is evaluated
! once here (cq1, cq2) instead of once per particle. Each particle then
! costs only its own B(J).

    do k=0,nquad
      aq1(k) = exp(-sin(0.5d0*quadQ(k))**2/sq1**2)
      aq2(k) = exp(-sin(0.5d0*quadQ(k))**2/sq2**2)
    end do

    do i=0,mode
      cq1(i) = 0.d0
      cq2(i) = 0.d0
      do k=0,nquad
        cq1(i) = cq1(i) + quadW(k)*aq1(k)*cos(dble(i)*quadQ(k))
        cq2(i) = cq2(i) + quadW(k)*aq2(k)*cos(dble(i)*quadQ(k))
      end do
      cq1(i) = quadnorm*cq1(i)
      cq2(i) = quadnorm*cq2(i)
    end do

    hk1 = (0.d0,0.d0)
    hk2 = (0.d0,0.d0)

! Accumulate the sum over particles. Each thread adds its share into local
! accumulators, combined afterwards in thread order, with a static split so
! every run gives each thread the same particles. An OpenMP reduction
! combines the parts in the order the threads finish, which varies from run
! to run and changed the last bits: two runs of the same binary differed by
! ~1e-15 in h_k (AUDITORIA_L0_2026-09-21.md, E17). Same scheme as
! VlasovPoisson_PIC_sp (e6d2011). The result still depends on the number of
! threads. The loop is split over particles rather than over modes because
! there are only five modes but 10^3-10^6 particles.

    nth = 1
!$  nth = omp_get_max_threads()
    allocate(parth(0:mode,2,0:nth-1))
    parth = (0.d0,0.d0)

    !$OMP PARALLEL PRIVATE(j,i,bj1,bj2,expv,lh1,lh2,tid)
    lh1 = (0.d0,0.d0)
    lh2 = (0.d0,0.d0)
    !$OMP DO SCHEDULE(STATIC)
    do j=1,Npart

      if (.not. bound(j)) cycle

!     J-dependent factor of each test function; the Q-quadrature is the
!     precomputed cq1/cq2.
      bj1 = exp(-(Jr(j)-j1)**2/sj1**2)*Jr(j)**2
      bj2 = exp(-(Jr(j)-j2)**2/sj2**2)*Jr(j)**2

      lh1(0) = lh1(0) + f(j)*bj1*cq1(0)
      lh2(0) = lh2(0) + f(j)*bj2*cq2(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        lh1(i) = lh1(i) + f(j)*bj1*cq1(i)*expv**i
        lh2(i) = lh2(i) + f(j)*bj2*cq2(i)*expv**i
      end do

    end do
    !$OMP END DO
    tid = 0
!$  tid = omp_get_thread_num()
    parth(:,1,tid) = lh1(0:mode)
    parth(:,2,tid) = lh2(0:mode)
    !$OMP END PARALLEL

    do tid=0,nth-1
      hk1(0:mode) = hk1(0:mode) + parth(:,1,tid)
      hk2(0:mode) = hk2(0:mode) + parth(:,2,tid)
    end do
    deallocate(parth)

    hk1 = drc*dpc*hk1
    hk2 = drc*dpc*hk2

! drc*dpc is the phase-space weight each particle represents; f holds the
! distribution function sampled there. The loop above only evaluates the
! double integral over (r,p_r) at fixed L, so the measure of the angular
! degrees of freedom removed by delta(L-L0), namely 8 pi^2 L0, is applied
! here to recover the physical h_k.
    hk1 = 8.0d0*smallpi**2*Lfix*hk1
    hk2 = 8.0d0*smallpi**2*Lfix*hk2
    abs_hk1 = abs(hk1)
    abs_hk2 = abs(hk2)


! *****************
! *** SAVE DATA ***
! *****************

! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim("hk1")//'.tl',form='formatted',status=filestatus)
     open(102,file=trim(directory)//'/'//trim("hk2")//'.tl',form='formatted',status=filestatus)
     open(103,file=trim(directory)//'/'//trim("hk1_complex")//'.tl',form='formatted',status=filestatus)
     open(104,file=trim(directory)//'/'//trim("hk2_complex")//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim("hk1")//'.tl',form='formatted',status=filestatus,position='append')
     open(102,file=trim(directory)//'/'//trim("hk2")//'.tl',form='formatted',status=filestatus,position='append')
     open(103,file=trim(directory)//'/'//trim("hk1_complex")//'.tl',form='formatted',status=filestatus,position='append')
     open(104,file=trim(directory)//'/'//trim("hk2_complex")//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(7ES24.16)") t,abs_hk1(:)
  write(102,"(7ES24.16)") t,abs_hk2(:)

! Columns: t, Re(h_0), Im(h_0), Re(h_1), Im(h_1), ... , Re(h_4), Im(h_4).
! The phase is kept because it is what allows discretisation noise to
! cancel when several independent runs are averaged.

  write(103,"(11ES24.16)") t,(real(hk1(k)),aimag(hk1(k)),k=0,mode)
  write(104,"(11ES24.16)") t,(real(hk2(k)),aimag(hk2(k)),k=0,mode)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)
  close(102)
  close(103)
  close(104)


  end subroutine analysish
