! ===========================================================================
! distribution.f90
! ===========================================================================
!> Initial distribution function in action-angle variables, F0(Q,J).
!!
!! The parameter "dftype" chooses which one is used, so the rest of the code
!! never mentions a particular functional form: the samplers in
!! initial_data.f90 only ask for F0 at a point and for an upper bound on it.
!!
!! All of them are written with the same convention for a width, exp(-x^2/s^2),
!! and are evaluated on Q in [0,2pi) and J >= 0.

module distribution

  use parameters

  implicit none

! --- "bimodal": two groups of actions, three angular harmonics ---
!
!   F0 = [1 + b1 cos Q + b2 cos 2Q] * [exp(-(J-Ja)^2/sa^2) + wb exp(-(J-Jb)^2/sb^2)]
!
! The two groups sit at different J, so they mix at different rates and h_1
! beats at their frequency difference. The angular part is a finite Fourier
! series: it has no harmonic above k=2, so h_3 and h_4 must vanish exactly.

  real(8), parameter :: bim_b1 = 0.8d0    ! amplitude of cos(Q)
  real(8), parameter :: bim_b2 = 0.3d0    ! amplitude of cos(2Q)
  real(8), parameter :: bim_Ja = 0.08d0   ! first group: centre in J
  real(8), parameter :: bim_sa = 0.03d0   !              width in J
  real(8), parameter :: bim_Jb = 0.22d0   ! second group: centre in J
  real(8), parameter :: bim_sb = 0.04d0   !               width in J
  real(8), parameter :: bim_wb = 0.7d0    ! weight of the second group

! --- "spiral": not separable, wound up already at t=0 ---
!
!   F0 = exp(-(J-J0)^2/s0^2) * exp(-sin((Q - beta J)/2)^2/sq^2)
!
! The angular profile is centred on Q = beta*J, so the distribution starts as
! a spiral in the (Q,J) plane. Phase mixing carries the centre to
! Q = beta*J + omega(J) t, and the spread of that across J closes at
! t = -beta/omega'(J0) > 0, where the distribution unwinds and |h_k| peaks
! before decaying for good.

  real(8), parameter :: spi_J0   = 0.15d0  ! centre in J
  real(8), parameter :: spi_s0   = 0.04d0  ! width in J
  real(8), parameter :: spi_beta = 50.0d0  ! winding of the initial spiral
  real(8), parameter :: spi_sq   = 0.5d0   ! width in Q

 contains

!> F0 at one point of the action-angle plane.

  function df0(Q,J) result(F)

    implicit none

    real(8), intent(in) :: Q,J
    real(8) :: F

    select case (dftype)

    case ('gauss')
       F = exp(-sin(0.5d0*Q)**2/sp**2)*exp(-J**2/sr**2)*J**2

    case ('bimodal')
       F = (1.0d0 + bim_b1*cos(Q) + bim_b2*cos(2.0d0*Q)) &
           *(exp(-(J-bim_Ja)**2/bim_sa**2) + bim_wb*exp(-(J-bim_Jb)**2/bim_sb**2))

    case ('spiral')
       F = exp(-(J-spi_J0)**2/spi_s0**2) &
           *exp(-sin(0.5d0*(Q-spi_beta*J))**2/spi_sq**2)

    case default
       F = 0.0d0
       call df0_unknown

    end select

  end function df0


!> An upper bound on F0 over the whole plane, for the rejection sampling.
!! It only has to be an upper bound: a loose one costs rejected trials, one
!! that is too small biases the sample.

  function df0_max() result(Fmax0)

    implicit none

    real(8) :: Fmax0

    select case (dftype)

    case ('gauss')
!      Exact maximum: at Q=0 and J=sr.
       Fmax0 = exp(-1.0d0)*sr**2

    case ('bimodal')
!      Each factor bounded separately: the angular part at Q=0, and the two
!      gaussians in J by their peaks.
       Fmax0 = (1.0d0 + bim_b1 + bim_b2)*(1.0d0 + bim_wb)

    case ('spiral')
!      Both factors are exponentials of a non-positive number.
       Fmax0 = 1.0d0

    case default
       Fmax0 = 0.0d0
       call df0_unknown

    end select

  end function df0_max


  subroutine df0_unknown

    implicit none

    print *
    print *, 'Unknown dftype: ',trim(dftype)
    print *, 'Aborting ...'
    print *
    stop

  end subroutine df0_unknown


!> Print which distribution is in use, with its constants, so that the run log
!! records the shape the particles were drawn from.

  subroutine df0_report

    implicit none

    print *
    print *, 'Initial distribution: ',trim(dftype)

    select case (dftype)
    case ('gauss')
       print '(a,es10.3,a,es10.3)', '   sigma_Q = ',sp,'   sigma_J = ',sr
    case ('bimodal')
       print '(a,f5.2,a,f5.2)', '   b1 = ',bim_b1,'   b2 = ',bim_b2
       print '(a,f6.3,a,f6.3)', '   Ja = ',bim_Ja,'   sa = ',bim_sa
       print '(a,f6.3,a,f6.3)', '   Jb = ',bim_Jb,'   sb = ',bim_sb
       print '(a,f5.2)',        '   wb = ',bim_wb
    case ('spiral')
       print '(a,f6.3,a,f6.3)', '   J0 = ',spi_J0,'   s0 = ',spi_s0
       print '(a,f7.2,a,f5.2)', '   beta = ',spi_beta,'   sigma_Q = ',spi_sq
    end select

    print *

  end subroutine df0_report


end module distribution
