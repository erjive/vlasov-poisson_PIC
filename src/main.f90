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


! Declare variables.

  implicit none

  integer i,j,k,l       ! Counters
  integer :: isub       ! Sub-step counter for the Yoshida composition
  real(8) :: dsub       ! Sub-step size for the Yoshida composition
! Yoshida (1990) 4th-order symplectic composition coefficients:
! one step = LF(w1*dt) o LF(w0*dt) o LF(w1*dt), with w1 = 1/(2-2^(1/3))
! and w0 = -2^(1/3)/(2-2^(1/3)). Note w0 < 0: the middle sub-step goes
! BACKWARDS in time, which is what cancels the O(dt^2) error term.
  real(8), parameter :: w1y =  1.0d0/(2.0d0-2.0d0**(1.0d0/3.0d0))
  real(8), parameter :: w0y = -2.0d0**(1.0d0/3.0d0)/(2.0d0-2.0d0**(1.0d0/3.0d0))


  call read_initial_param()

!  call test_consistency()

  call set_grid_size()

  call alloc_mem_set0()

  call construct_grid()

  call initial_data()

! ****************************************************
! ***   ACTION-ANGLE SETUP FOR THE EXACT ADVANCE   ***
! ****************************************************

! integrator="analytic" advances each particle with the closed-form
! solution Q3(t)=Q3(0)+omega(J3)*t, which only exists while J3 is exactly
! conserved: static background, fixed L /= 0, no self-interaction, and no
! centrifugal softening (eps must be 0, or the dynamics would not match
! the unsoftened action-angle map -- the very inconsistency documented in
! BUGS_TODO.md).

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

!    reduce_arrays reallocates r_part/p_part/f with a smaller Npart but
!    knows nothing about q0_part/j0_part, so the per-particle arrays would
!    silently get out of step with each other.
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
  call system('cp input_parameters '//trim(directory))

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

!    Save old time step.

      r_part_p = r_part

      p_part_p = p_part

 
!    Euler method (forward differencing in time, first order).

     if (integrator=='euler') then

       r_part = r_part_p + p_part*dt
       p_part = p_part_p + force_part*dt

       call grav_force()
!   Second order leapfrog method

    else if (integrator == 'leapfrog') then

!     Leapfrog integration 'kick-drift-kick' form

      p_part_h = p_part_p + force_part*dt*0.5D0
      r_part   = r_part_p + p_part_h * dt

      call  grav_force()

      p_part   = p_part_h + force_part*dt*0.5D0

!   Fourth order symplectic (Yoshida composition of three leapfrog steps)

    else if (integrator == 'yoshida4') then

!     Three kick-drift-kick sub-steps with sizes w1*dt, w0*dt, w1*dt.
!     Phase error drops from O(dt^2) to O(dt^4) at the cost of 3 force
!     evaluations per step (3 Poisson solves per step when
!     autointeraction=.true., where that path already dominates).
!     Symplectic, unlike the rk4 branch below, so there is still no
!     secular drift in the energy or in J3.

      do isub = 1,3

        if (isub == 2) then
          dsub = w0y*dt
        else
          dsub = w1y*dt
        end if

        p_part_h = p_part   + force_part*dsub*0.5D0
        r_part   = r_part   + p_part_h  *dsub

        call grav_force()

        p_part   = p_part_h + force_part*dsub*0.5D0

      end do

!   Exact analytic advance (no self-interaction only)

    else if (integrator == 'analytic') then

!     At fixed L the radial motion in a static background is integrable:
!     J3 is exactly conserved and Q3(t) = Q3(0) + omega(J3)*t. Advance to
!     the ABSOLUTE time t (not incrementally), so there is no accumulated
!     phase error of any kind. Validation path: isolates the quadrature,
!     analysish and the normalisations from all integrator error.

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

    do i=1,Npart

      if (rmin == 0 .and. r_part(i)<0.d0) then

        r_part(i) = -r_part(i)
        p_part(i) = -p_part(i)

      end if

    end do

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
