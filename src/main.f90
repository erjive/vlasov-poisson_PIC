! ************************
! ***   MAIN PROGRAM   ***
! ************************

program VP_PIC

! Include modules

  use parameters
  use arrays
  use utils
  use hdf5_io
  use raw_io
  use paramfile


! Declare variables.

  implicit none

  integer i,j,k,l       ! Counters
  integer :: isub       ! Sub-step counter for the Yoshida composition
  real(8) :: dsub       ! Sub-step size for the Yoshida composition
  integer :: nstage     ! Number of leapfrog sub-steps in the composition
  real(8) :: wcomp(15)  ! Composition coefficients (symmetric, sum = 1)

! --- Yoshida (1990) symmetric compositions of leapfrog sub-steps -------
!
! One step = LF(w_1 dt) o LF(w_2 dt) o ... o LF(w_n dt), symmetric, with
! sum(w) = 1. Some w_i are NEGATIVE: those sub-steps run backwards in
! time, and that is exactly what cancels the lower-order error terms.
!
! 4th order, 3 stages: w1 = 1/(2-2^(1/3)), w0 = -2^(1/3)/(2-2^(1/3)).
  real(8), parameter :: y4a =  1.0d0/(2.0d0-2.0d0**(1.0d0/3.0d0))
  real(8), parameter :: y4b = -2.0d0**(1.0d0/3.0d0)/(2.0d0-2.0d0**(1.0d0/3.0d0))
!
! 6th order, 7 stages: Yoshida (1990) "Solution A". The central weight is
! fixed by the consistency condition w0 = 1 - 2*(w1+w2+w3).
  real(8), parameter :: y6c = -1.17767998417887d0
  real(8), parameter :: y6b =  0.235573213359357d0
  real(8), parameter :: y6a =  0.784513610477560d0
  real(8), parameter :: y6z =  1.0d0 - 2.0d0*(y6a+y6b+y6c)


  call read_parameters()

!  call test_consistency()

  call set_grid_size()

  call alloc_mem_set0()

  call construct_grid()

  call initial_data()

! ****************************************************
! ***   ACTION-ANGLE SETUP FOR THE EXACT ADVANCE   ***
! ****************************************************

! integrator="analytic" advances each particle with the closed-form
! solution Q3(t) = Q3(0) + omega(J3)*t, valid only where the radial action
! J3 is an exact constant of motion: static isochrone background, fixed
! L /= 0 and no self-interaction. The centrifugal softening must also be
! off (eps = 0), since the action-angle map is built from the unsoftened
! L^2/(2 r^2) term.

  if (integrator == 'analytic') then

     if (autointeraction .or. forcetype /= "bg" .or. BGtype /= "Isochrone" &
         .or. Lfix == 0.0d0 .or. eps /= 0.0d0) then
        print *
        print *, 'integrator="analytic" requires an integrable setup:'
        print *, '  autointeraction = .false.,  forcetype = "bg",'
        print *, '  BGtype = "Isochrone",  Lfix /= 0,  eps = 0.'
        print *, 'Aborting ...'
        print *
        stop
     end if

!    The analytic advance needs the initial (Q3,J3) of every particle, and
!    reduce_arrays resizes r_part/p_part/f without touching q0_part/j0_part.
     if (reduceparticles) then
        print *
        print *, 'integrator="analytic" is incompatible with reduceparticles=.true.'
        print *, '(reduce_arrays would resize r_part/p_part but not q0_part/j0_part).'
        print *, 'Aborting ...'
        print *
        stop
     end if

     call init_action_angle()

  end if

! ***************************
! ***   OUTPUT DIRECTORY  ***
! ***************************

! Create output directory and copy parameter file to it.

  call system('mkdir -p '//trim(directory))
  call system('cp '//trim(parameter_file)//' '//trim(directory))

  if (output_format=="hdf5") call open_hdf5_file()
  if (output_format=="raw")  call open_raw_file()

! Initialize time.

  t = 0.0d0

! **************************************
! ***   FIND DENSITY AND FLUX IN r   ***
! **************************************
  call density

! **************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE   ***
! **************************************************
  call grav_force()

! **************************************************************
! ***   FIND TOTAL NUMBER OF PARTICLES                       ***
! ***   AVERAGE KINETIC ENERGY, POTENTIAL AND TOTAL ENERGY   ***
! **************************************************************

!  call integrate
   call energy

! ************************
! *** INITIAL ANALYSIS ***
! ************************

   call analysish

!  ***************************
!  ***  INITIAL TIME STEP  ***
!  ***************************
  call set_timestep()

  print *
  print *, 'Time step fixed at size: ',dt

! *****************************
! ***   OUTUPUT TO SCREEN   ***
! *****************************

  print *
  print *,'------------------------------'
  print *,'|  Time step  |     Time     |'
  print *,'------------------------------'

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************
   if (output_format=="hdf5") then
      call save_data_hdf5(0)
   else if (output_format=="raw") then
      call save_data_raw(0)
   else
      call save_data()
   end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Euler method (forward differencing in time, first order).

     if (integrator=='euler') then

       r_part_p = r_part
       p_part_p = p_part

       r_part = r_part_p + p_part*dt
       p_part = p_part_p + force_part*dt

       call grav_force()
!   Second order leapfrog method

    else if (integrator == 'leapfrog') then

!     Leapfrog in 'kick-drift-kick' form: half kick, full drift with the
!     half-step momentum, force update, half kick. Second order and
!     symplectic, so the energy error stays bounded instead of drifting.
!
!     The kick and the drift are fused into one loop (and the final kick
!     into another) so each particle is read and written once per force
!     evaluation; the drift is done in place, so no copy of the old
!     positions is needed.

      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i=1,Npart
        p_part_h(i) = p_part(i)   + force_part(i)*dt*0.5D0
        r_part(i)   = r_part(i)   + p_part_h(i)  *dt
      end do
      !$OMP END PARALLEL DO

      call  grav_force()

      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i=1,Npart
        p_part(i)   = p_part_h(i) + force_part(i)*dt*0.5D0
      end do
      !$OMP END PARALLEL DO

!   Fourth order symplectic (Yoshida composition of three leapfrog steps)

    else if (integrator == 'yoshida4' .or. integrator == 'yoshida6') then

!     Table-driven symmetric composition of leapfrog sub-steps. Each
!     sub-step is a full kick-drift-kick of size w_i*dt, so the cost is
!     nstage force evaluations per step (nstage Poisson solves per step
!     when autointeraction=.true.). Symplectic at every order, unlike
!     the rk4 branch below, so no secular drift is reintroduced.
!
!     Raising the order is not automatically cheaper: reducing the error
!     by a factor R needs dt/R^(1/p) and therefore costs ~nstage*R^(1/p)
!     force evaluations, so the 7-stage 6th order only beats the 3-stage
!     4th order when very small errors are demanded.

      if (integrator == 'yoshida4') then
        nstage = 3
        wcomp(1:3) = (/ y4a, y4b, y4a /)
      else
        nstage = 7
        wcomp(1:7) = (/ y6a, y6b, y6c, y6z, y6c, y6b, y6a /)
      end if

      do isub = 1,nstage

        dsub = wcomp(isub)*dt

        !$OMP PARALLEL DO SCHEDULE(STATIC)
        do i=1,Npart
          p_part_h(i) = p_part(i)   + force_part(i)*dsub*0.5D0
          r_part(i)   = r_part(i)   + p_part_h(i)  *dsub
        end do
        !$OMP END PARALLEL DO

        call grav_force()

        !$OMP PARALLEL DO SCHEDULE(STATIC)
        do i=1,Npart
          p_part(i)   = p_part_h(i) + force_part(i)*dsub*0.5D0
        end do
        !$OMP END PARALLEL DO

      end do

!   Exact analytic advance (no self-interaction only)

    else if (integrator == 'analytic') then

!     At fixed L the radial motion in a static background is integrable:
!     J3 is conserved and Q3(t) = Q3(0) + omega(J3)*t. Each particle is
!     advanced to the absolute time t rather than step by step, so no
!     phase error accumulates and the diagnostics can be checked
!     independently of any integrator error.

      call advance_analytic(t)
      call grav_force()

!    Fourth order Runge-Kutta.

     else if (integrator=='rk4') then

        print *
        print *, 'Fourth order Runge-Kutta is deliberately NOT implemented:'
        print *, 'RK4 is not symplectic, so it reintroduces the secular drift'
        print *, 'in energy and J3 that floors h_k. Use integrator="yoshida4"'
        print *, '(4th order, symplectic) instead.'
        print *, 'Aborting ...'
        print *
        stop

!    Unknown integration method.

     else

        print *, 'Unknown integration method.'
        print *, 'Aborting ...'
        print *
        stop

     end if

!   At the origin impose symmetry condition f(r,p) = f(-r,-p)

!   "rmin == 0" is a run constant, so it is tested once here instead of
!   once per particle per step, and the remaining loop is parallel. When
!   rmin > 0 the whole pass is skipped outright (it could never fire).

    if (rmin == 0.0d0) then

      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do i=1,Npart
        if (r_part(i)<0.d0) then
          r_part(i) = -r_part(i)
          p_part(i) = -p_part(i)
        end if
      end do
      !$OMP END PARALLEL DO

    end if

!    **************************************
!    ***   FIND DENSITY AND FLUX IN r   ***
!    **************************************

     if (forcetype=="self") then

!    Save old value of rho and curr.
        rho_p  = rho
!        curr_p = curr

!    Integrate over phase space.
        !call density

     else

        if (mod(l+1,spatial_output).eq.0) then

!    Integrate over phase space.
!           call density

!    Save old value of rho and curr in order to calculate the continuity equation
!           rho_p  = rho
!           curr_p = curr

        end if

        if (mod(l,spatial_output).eq.0) then

!    Integrate over phase space.
           call density
           call energy
        end if

     end if


! ******************************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE ON THE PARTICLES  ***
! ******************************************************************
     !call grav_force()


!    **********************************
!    ***   SOLVE POISSON EQUATION   ***
!    **********************************

!    For the self gravitating case solve
!    the Poisson equation again.

!     if (forcetype=="self") then
!        !call poisson
!     end if



!    *****************************************
!    ***   CALCULATE CONTINUITY EQUATION   ***
!    *****************************************

!    The continuity equation has the form:
!
!    cont  =  0  =  d(rho)/dt + div(curr)  = d(rho)/dt + (1/r**2) d(r**2 curr)/dr
!
!                =  d(rho)/dt + d(curr)/dr + 2 curr / r
!
!    Notice that this should converge to zero.  The expression below is only
!    second order accurate.

!     if (mod(l,spatial_output).eq.0) then
!        do i=1,Nr-1
!           cont(i) = (rho(i) - rho_p(i))/dt &
!                + 0.25d0*(curr(i+1) + curr_p(i+1) - curr(i-1) - curr_p(i-1))/dr &
!                + (curr(i) + curr_p(i))/r(i)
!        end do
!     end if


!    ***************************
!    ***   ADAPT TIME STEP   ***
!    ***************************

!    For the self-gravitating case the force can change
!    with time, so one needs to adapt the time step.
!    Notice that the time step can go up and down in
!    response to the size of the force.

     !if (forcetype=="self") then
     !  call set_timestep()
     !end if


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

! field_output (independent of spatial_output) gates the expensive
! r_part/p_part/f snapshot -- the bulk of the disk footprint -- so
! hk1.tl/hk1_complex.tl (via analysish below) can be sampled finely in
! time without paying for an equally frequent field dump.
     if (mod(l,field_output).eq.0) then

       if (output_format=="hdf5") then
          call save_data_hdf5(l)
       else if (output_format=="raw") then
          call save_data_raw(l)
       else
          call save_data()
       end if

     end if

     if (mod(l,spatial_output).eq.0) then

        call analysish

     end if

!    *************************************************
!    ***   IF POSSIBLE REDUCE SIZE OF THE ARRAYS   ***
!    *************************************************

     if (reduceparticles .and. (mod(l,Nreduce).eq.0)) then
!     if (mod(l,time_output).eq.0) then
       call reduce_arrays

     end if

!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

     if (mod(l,time_output).eq.0) then
        write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',l,'   | ',t,'  | '
     end if

  end do

  print *,'------------------------------'


! ***************
! ***   END   ***
! ***************
  print *, 'Maximum radii of particles = ', maxval(r_part)

  if (output_format=="hdf5") call close_hdf5_file()
  if (output_format=="raw")  call close_raw_file()

  call deallocate_mem()

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *
  print *

end program VP_PIC
