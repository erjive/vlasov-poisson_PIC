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
    real(8) :: gval1(0:nquad),gval2(0:nquad)  !mode-independent integrand for Phi_1/Phi_2, per node
    real(8) :: contrib1(0:4),contrib2(0:4)    !phik(.,mode) for mode=0..4, current particle, Phi_1/Phi_2

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

    hk1 = (0.d0,0.d0)
    hk2 = (0.d0,0.d0)

! phik(J,l,J0,sq,sj) = (1/pi) * Simpson[ g(J,J0,sq,sj,Q)*cos(l*Q) dQ, Q=0..pi ],
! with g(...,Q) = exp(-sin(Q/2)^2/sq^2)*exp(-(J-J0)^2/sj^2)*J^2 the part
! of the original phi(J,Q,l,J0,sq,sj) integrand that does NOT depend on
! l (the mode). The previous implementation called phik() once per
! (particle, mode, test function) -- (mode+1)*2 = 10 times per particle
! -- and each call recomputed g(Q) at all 513 quadrature points from
! scratch, i.e. 10x more exp() evaluations than necessary, since g only
! depends on the particle (through Jr(j)) and on which test function
! (Phi_1 or Phi_2) is being evaluated, not on which of the 5 modes is
! being accumulated. Computing gval1/gval2 once per particle and
! reusing them for all 5 modes removes that redundancy. It also
! replaces the original complex phi(Q) = g(Q)*exp(-i*l*Q) with
! g(Q)*cos(l*Q) directly: phik's own result was always real anyway (the
! previous code built it from "real(auxsum)*2", silently discarding the
! imaginary part of the Simpson sum every time -- see the symmetry
! argument in BUGS_TODO.md), so the sin(l*Q) part it implicitly threw
! away is simply never computed now.
!
! Each particle's contribution to hk1(0:4)/hk2(0:4) is independent and
! only summed, so this parallelizes over particles with a plain
! reduction on the two (tiny, 5-element) accumulators -- no atomics
! needed. This also replaces the disabled "!!$OMP" that used to
! parallelize the outer do-i-over-modes loop (only 5 iterations, a poor
! match for 8 threads): parallelizing over particles instead scales
! with Npart, which is what actually matters at N_c ~ 10^3-10^5.

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,k,i,gval1,gval2,contrib1,contrib2,expv) REDUCTION(+:hk1,hk2)
    do j=1,Npart

      do k=0,nquad
        gval1(k) = exp(-sin(0.5d0*quadQ(k))**2/sq1**2)*exp(-(Jr(j)-j1)**2/sj1**2)*Jr(j)**2
        gval2(k) = exp(-sin(0.5d0*quadQ(k))**2/sq2**2)*exp(-(Jr(j)-j2)**2/sj2**2)*Jr(j)**2
      end do

      do i=0,mode
        contrib1(i) = 0.d0
        contrib2(i) = 0.d0
        do k=0,nquad
          contrib1(i) = contrib1(i) + quadW(k)*gval1(k)*cos(dble(i)*quadQ(k))
          contrib2(i) = contrib2(i) + quadW(k)*gval2(k)*cos(dble(i)*quadQ(k))
        end do
        contrib1(i) = quadnorm*contrib1(i)
        contrib2(i) = quadnorm*contrib2(i)
      end do

      hk1(0) = hk1(0) + f(j)*contrib1(0)
      hk2(0) = hk2(0) + f(j)*contrib2(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        hk1(i) = hk1(i) + f(j)*contrib1(i)*expv**i
        hk2(i) = hk2(i) + f(j)*contrib2(i)*expv**i
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
  else
     open(101,file=trim(directory)//'/'//trim("hk1")//'.tl',form='formatted',status=filestatus,position='append')
     open(102,file=trim(directory)//'/'//trim("hk2")//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(7ES16.8)") t,abs_hk1(:)
  write(102,"(7ES16.8)") t,abs_hk2(:)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)
  close(102)


  end subroutine analysish
