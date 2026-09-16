  subroutine analysish

    use parameters
    use arrays
    !use utils


    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8) :: energy,s, s1, s2, er1, er2,argaux
    complex(8) :: ii                          !imaginary unit
    complex(8) :: expv                        !exp(-ii*Qr(j)) for the current particle
    complex(8),dimension(0:4) :: hk1,hk2           !h_k mode
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


    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,energy,s, s1, s2, er1, er2,argaux) SHARED(Qr,Jr)
    do i = 1,Npart
      energy = -1.0/(1.0D0+dsqrt(1.0D0+r_part(i)**2)) + 0.5d0*Lfix**2/(r_part(i)**2) + 0.5D0*p_part(i)**2
      er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
      er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
      s1 = 1.d0 + dsqrt(1.d0+er1**2)
      s2 = 1.d0 + dsqrt(1.d0+er2**2)
      s  = 1.d0 + dsqrt(1.d0+r_part(i)**2)
      argaux = (s1+s2-2.0*s)/(s2-s1)
      Jr(i) = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))

  !  do i=1,Npart

          if (p_part(i)>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0D0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0D0),argaux))+smallpi
          end if

          Qr(i) = eta - dsqrt((-2.d0*energy)**3)*sqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*dsin(eta)

    end do
    !$OMP END PARALLEL DO

    ii = (0.d0,1.d0)

! Quadrature nodes/weights for Simpson's rule on [0,pi] with nquad
! (even) sub-intervals -- the same grid for every particle and every
! mode, so build it once instead of inside phik() on every one of the
! (mode+1)*Npart*2 calls it used to get (this loop's own cost is
! negligible, O(nquad), done once per analysish() call).

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

! The integrand of phik factorises exactly:
!
!   g(Q,J) = exp(-sin(Q/2)^2/sq^2) * exp(-(J-J0)^2/sj^2) * J^2
!          =        A(Q)           *            B(J)
!
! so the whole Q-quadrature separates from the particle:
!
!   phik(J,i) = quadnorm * sum_k w_k A(Q_k) B(J) cos(i Q_k)
!             = B(J) * [ quadnorm * sum_k w_k A(Q_k) cos(i Q_k) ]
!             = B(J) * C(i)
!
! and C(i) does not depend on the particle at all -- it is a constant of
! the whole run. It used to be recomputed inside the particle loop, so
! every particle paid a full (nquad+1)*(mode+1) quadrature: 2*513 exp,
! 2*513 sin and 2*5*513 cos/multiply-add EACH, ~500x more arithmetic
! than needed. Hoisted out here, each particle now costs 2 exp and 2*5
! multiplies. (Same class of redundancy as the Wn/Sn duplication noted
! in BUGS_TODO.md, but that one measured ~1x and this one ~500x.)
!
! Note this changes the summation order, so results are no longer
! bit-identical to the previous version -- they agree to roundoff
! (~1e-16 relative, verified), far below the integrator and quadrature
! errors that actually limit h_k.

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

! Each particle's contribution to hk1(0:4)/hk2(0:4) is independent and
! only summed, so this parallelizes over particles with a plain
! reduction on the two (tiny, 5-element) accumulators -- no atomics
! needed. Parallelizing over particles (rather than over the 5 modes, as
! an earlier disabled directive did) scales with Npart, which is what
! actually matters at N_c ~ 10^3-10^5.
!
! Historical note: phik() used to be a function called once per
! (particle, mode, test function) that rebuilt the whole Q-quadrature on
! every call. Two rounds of cleanup removed that: first hoisting the
! mode-independent part out of the mode loop, then (see the
! factorisation above) hoisting the entire Q-quadrature out of the
! particle loop, which is where the 10x speedup came from.

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,i,bj1,bj2,expv) REDUCTION(+:hk1,hk2)
    do j=1,Npart

!     Only the J-dependent factor is per-particle now; the Q-quadrature
!     lives in cq1/cq2, computed once above.
      bj1 = exp(-(Jr(j)-j1)**2/sj1**2)*Jr(j)**2
      bj2 = exp(-(Jr(j)-j2)**2/sj2**2)*Jr(j)**2

      hk1(0) = hk1(0) + f(j)*bj1*cq1(0)
      hk2(0) = hk2(0) + f(j)*bj2*cq2(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        hk1(i) = hk1(i) + f(j)*bj1*cq1(i)*expv**i
        hk2(i) = hk2(i) + f(j)*bj2*cq2(i)*expv**i
      end do

    end do
    !$OMP END PARALLEL DO

    hk1 = drc*dpc*hk1
    hk2 = drc*dpc*hk2

! Ec. 44 del paper, reducida via F(r,pr,L) = F0(r,pr)*delta(L-L0), trae un
! factor global 8*pi^2*L0 que el bucle de arriba no incluye (solo calcula
! la integral doble en (r,pr) a L0 fijo). Sin este factor, hk1/hk2 son
! proporcionales al h_k fisico real, no iguales -- confirmado comparando
! contra los datos originales del articulo (ver
! VlasovPoisson_PIC_sp/Vlasov_Poisson_evolutions/h0_normalization_check.md
! SS9-10): con el factor aplicado, h_0 coincide al 0.01% con el valor
! "Analytical" publicado.
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

! hk1/hk2 above only ever save the magnitude |h_k|, which is fine for
! looking at a single run but useless for telling real decay-to-zero
! apart from a discreteness-noise floor that never shrinks: the noise
! has a random phase from run to run (different seed), so averaging
! |h_k| over several independent runs does NOT cancel it (same effect
! as Rayleigh-distributed magnitude noise never averaging below its
! own scale) -- only averaging the COMPLEX h_k first, then taking the
! magnitude of THAT average, lets random-phase noise cancel while a
! real physical (phase-coherent) signal survives. Saving Re/Im here
! (columns: t, Re(h_0),Im(h_0), Re(h_1),Im(h_1), ..., Re(h_4),Im(h_4))
! makes that kind of ensemble averaging possible in post-processing.
! See BUGS_TODO.md.

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
