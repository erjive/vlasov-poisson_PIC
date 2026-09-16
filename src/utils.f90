! ===========================================================================
! utils.f90
! ===========================================================================
!> Utility file


module utils

  use parameters
  use arrays

 contains



! Sanity check.
  subroutine test_consistency

  if (rmin<0.0D0) then
     print *
     print *, 'rmin must be greater than or equal to zero'
     print *, 'Aborting ...'
     print *
     stop
  end if

  if (rmax<=rmin) then
     print *
     print *, 'rmin must be smaller than rmax'
     print *, 'Aborting ...'
     print *
     stop
  end if

!  if (pmax<=0.d0) then
!     print *
!     print *, 'pmax must be greater then zero'
!     print *, 'Aborting ...'
!     print *
!     stop
!  end if

!  if (pmax<=pmin) then
!     print *
!     print *, 'pmin must be smaller than pmax'
!     print *, 'Aborting ...'
!     print *
!     stop
!  end if

  end subroutine test_consistency



  !> Set all the parameters that are not set in the input file
  subroutine set_grid_size
! ***************************
! ***   FIND GRID SIZES   ***
! ***************************

    ! Find out number of grid points in r direction.
    Nr = int((rmax-rmin)/dr)
    if (rmin==0.0d0) Nr=Nr+1

    ! Find out number of grid points in p direction.
    !Np = 2*int(pmax/dp)

    print *
    print *, 'Number of points in r direction ',Nr
    !print *, 'Number of points in p direction ',Np

  end subroutine set_grid_size

  !> Allocate all memory
  subroutine alloc_mem_set0

! Find the number of ghost zones

  if (rmin>0.d0) then
    ghost = 0
  else if (rmin == 0.d0) then
    if (spatialorder == "two") then
        ghost = 2
    else if (spatialorder == "four") then
        ghost = 3
    end if
  end if

! Find the number of particles

  Npart = Nrc*Npc
  print *, "Number of computational particles=",Npart
! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Position and momentum of particles

  allocate(r_part      (1:Npart))
  allocate(r_part_p    (1:Npart))
  allocate(p_part      (1:Npart))
  allocate(p_part_p    (1:Npart))
  allocate(p_part_h   (1:Npart))
  allocate(pot_part    (1:Npart))
  allocate(force_part  (1:Npart))

  r_part     = 0.d0
  r_part_p   = 0.d0
  p_part     = 0.d0
  p_part_p   = 0.d0
  p_part_h   = 0.d0
  pot_part   = 0.d0
  force_part = 0.d0

! Coordinates, force and potential.

! construct_grid fills r(i) for i=0,Nr when rmin>0 (ghost=0), but for
! i=1-ghost,Nr when rmin=0 (ghost=2 or 3) -- the rmin>0 case still
! writes index 0, one element below the "1-ghost=1" lower bound that
! formula alone would give, so allocate the wider of the two ranges
! explicitly instead of just "r(1-ghost:Nr)" (an out-of-bounds write
! to r(0) for rmin>0, silently corrupting the heap until something
! else's deallocate happens to detect it).

  if (rmin>0.d0) then
     allocate(r(0:Nr))
  else
     allocate(r(1-ghost:Nr))
  end if
  r = 0.0d0

  if (autointeraction) then
    allocate(force  (1-ghost:Nr))
    allocate(pot    (1-ghost:Nr))
    allocate(dev_pot(1-ghost:Nr))
    force   = 0.0d0
    pot     = 0.0d0
    dev_pot = 0.0d0
  end if


! Density function f and sources.

  allocate(f  (1:Npart))
  !allocate(f_p(1:Npart))

  f   = 0.0d0
  !f_p = 0.0d0

! Integrates density, current and continuity equation.

  allocate(rho   (1-ghost:Nr))
  allocate(rho_p (1-ghost:Nr))
  allocate(avg_rho (1-ghost:Nr))
  allocate(curr  (1-ghost:Nr))
  allocate(curr_p(1-ghost:Nr))
  allocate(cont  (1-ghost:Nr))

  rho    = 0.0d0
  rho_p  = 0.0d0
  avg_rho= 0.0d0
  curr   = 0.0d0
  curr_p = 0.0d0
  cont   = 0.0d0

  end subroutine alloc_mem_set0


  !> Free all the memory in the allocated arrays.
  !!
  !! Mirrors alloc_mem_set0: every array is released under the same
  !! condition that allocated it, since deallocating an unallocated
  !! array is a runtime error.
  subroutine deallocate_mem

! q0_part/j0_part exist only for integrator="analytic".
  if (integrator == 'analytic') then
    deallocate(q0_part)
    deallocate(j0_part)
  end if

  deallocate(r_part)
  deallocate(r_part_p)
  deallocate(p_part)
  deallocate(p_part_p)
  deallocate(p_part_h)
  deallocate(pot_part)
  deallocate(force_part)

  deallocate(r)

! The grid fields force, pot and dev_pot exist only in the
! self-gravitating case. "res" is declared in arrays.f90 but never
! allocated, so it is not released here.

  if (autointeraction) then
    deallocate(force)
    deallocate(pot)
    deallocate(dev_pot)
  end if

! Density function f

  deallocate(f)
!  deallocate(f_p)

  deallocate(rho)
  deallocate(rho_p)
  deallocate(avg_rho)
  deallocate(curr)
  deallocate(curr_p)
  deallocate(cont)

  print *, "Memory deallocated"

  end subroutine deallocate_mem


subroutine construct_grid
! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Position in direction r.  In case that rmin=0
! we make sure to stagger the origin and we add
! two ghost points to the left of the origin for 
! second order integration, and three points for 
! fourth order integration.

  integer i

  if (rmin>0.d0) then
     do i=0,Nr
        r(i) = rmin + dble(i)*dr
     end do
  else

  do i=1-ghost,Nr
    r(i) = (dble(i)-0.5d0)*dr
  end do
      
  end if

! Position in direction p. Notice that the grid in the
! p direction will always cover the region (-pmax,pmax).

!  do j=0,Np
!     p(j) = pmin + dble(j) *dp
!  end do


end subroutine construct_grid

  !> Build a cell list for the current particles.
  !!
  !! Groups particle indices (1:Npart) by the radial grid cell
  !! (spacing dr) nearest their position r_part(j), so that a
  !! deposit loop over grid points i=1,Nr only needs to scan
  !! particles in a few nearby cells instead of all Npart particles
  !! -- turning the O(Nr*Npart) brute-force search in density() and
  !! avg_density() into ~O(Nr+Npart).
  !!
  !! On output, the particles assigned to grid cell c (1<=c<=Nr) are
  !! particle_order(cell_start(c):cell_start(c+1)-1). Both arrays
  !! are allocated here; the caller must deallocate them.
  !!
  !! The grid is uniform with spacing dr for both grid conventions
  !! used in construct_grid (rmin>0 and the staggered rmin=0 case),
  !! so r(k) = r(1) + (k-1)*dr always holds, and the nearest grid
  !! index to a position x is nint((x-r(1))/dr) + 1, independent of
  !! which convention built the grid.  Particles that fall outside
  !! the physical range are clamped into the boundary cell: harmless,
  !! since the caller still applies the exact distance cutoff and
  !! will simply reject them.
  !> Invert one (Q3,J3) pair back to (r,p_r) at fixed L=Lfix.
  !!
  !! Same Kepler-like inversion used by initial_data's "aa_quad" state:
  !! J3 fixes the energy, the energy fixes the turning points, and the
  !! angle Q3 is mapped to the eccentric-anomaly-like variable eta by
  !! Newton-Raphson on Q3 = eta - ecc*sin(eta). Factored out here so the
  !! analytic integrator and the initial data share one implementation.
  subroutine invert_QJ_to_rp(Qv,Jv,rv,pv)

    implicit none

    real(8), intent(in)  :: Qv,Jv
    real(8), intent(out) :: rv,pv

    real(8) :: Eg,er1,er2,s1,s2,ecc,etaNR,gNR,gpNR,argaux,sg,pv2,smallpi
    integer :: it

    smallpi = acos(-1.0d0)

    Eg = -1.d0/(2.d0*(Jv+0.5d0*(Lfix+sqrt(Lfix**2+4.d0)))**2)

    er1 = dsqrt((1.d0+Eg*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*Eg*(2.d0+2.d0*Eg+Lfix**2)))/(2.d0*Eg**2))
    er2 = dsqrt((1.d0+Eg*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*Eg*(2.d0+2.d0*Eg+Lfix**2)))/(2.d0*Eg**2))
    s1 = 1.d0 + sqrt(1.d0+er1**2)
    s2 = 1.d0 + sqrt(1.d0+er2**2)

    ecc = sqrt((-2.d0*Eg)**3)*sqrt(-Lfix**2-2.d0*Eg-2.d0-0.5D0/Eg)/(-2.d0*Eg)

    etaNR = modulo(Qv,2.0d0*smallpi)
    do it=1,50
      gNR  = etaNR - ecc*sin(etaNR) - modulo(Qv,2.0d0*smallpi)
      gpNR = 1.d0 - ecc*cos(etaNR)
      etaNR = etaNR - gNR/gpNR
      if (abs(gNR) < 1.0d-14) exit
    end do

    argaux = cos(etaNR)
    sg = (s1+s2-argaux*(s2-s1))/2.0d0
    rv = sqrt(max((sg-1.d0)**2-1.d0,0.0d0))

    pv2 = 2.d0*(Eg + 1.d0/(1.d0+dsqrt(1.d0+rv**2)) - 0.5d0*Lfix**2/max(rv**2,1.0d-12))
    pv2 = sqrt(max(pv2,0.0d0))
    if (modulo(etaNR,2.0d0*smallpi) > smallpi) then
      pv = -pv2
    else
      pv =  pv2
    end if

  end subroutine invert_QJ_to_rp


  !> Store each particle's initial (Q3,J3), for integrator="analytic".
  !!
  !! Uses the same forward map as analysish.f90, so it works from any
  !! initial state, not just the ones built in (Q3,J3) to begin with.
  subroutine init_action_angle

    implicit none

    integer :: i
    real(8) :: en,er1,er2,s1,s2,ss,argaux,eta,smallpi

    smallpi = acos(-1.0d0)

    allocate(q0_part(1:Npart))
    allocate(j0_part(1:Npart))

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(en,er1,er2,s1,s2,ss,argaux,eta)
    do i=1,Npart
      en = -1.0/(1.0D0+dsqrt(1.0D0+r_part(i)**2)) + 0.5d0*Lfix**2/(r_part(i)**2) + 0.5D0*p_part(i)**2
      er1 = dsqrt((1.d0+en*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*en*(2.d0+2.d0*en+Lfix**2)))/(2.d0*en**2))
      er2 = dsqrt((1.d0+en*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*en*(2.d0+2.d0*en+Lfix**2)))/(2.d0*en**2))
      s1 = 1.d0 + dsqrt(1.d0+er1**2)
      s2 = 1.d0 + dsqrt(1.d0+er2**2)
      ss = 1.d0 + dsqrt(1.d0+r_part(i)**2)
      argaux = (s1+s2-2.0d0*ss)/(s2-s1)
      if (p_part(i)>=0.d0) then
        eta = dacos(sign(min(abs(argaux),1.0D0),argaux))
      else
        eta = dacos(-sign(min(abs(argaux),1.0D0),argaux))+smallpi
      end if
      q0_part(i) = eta - dsqrt((-2.d0*en)**3)*sqrt(-Lfix**2-2.d0*en-2.d0-0.5D0/en)/(-2.d0*en)*dsin(eta)
      j0_part(i) = 1.d0/dsqrt(-2.d0*en)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))
    end do
    !$OMP END PARALLEL DO

  end subroutine init_action_angle


  !> Advance every particle ANALYTICALLY to absolute time "tnow".
  !!
  !! Without self-interaction the radial motion at fixed L is integrable:
  !! J3 is conserved and the angle advances linearly,
  !! Q3(t) = Q3(0) + omega(J3)*t with omega(J) = 1/(J+c)^3 and
  !! c = (L + sqrt(L^2+4))/2 for the isochrone. The particle state is then
  !! recovered by inverting (Q3,J3) back to (r,p_r).
  !!
  !! This is exact at any t, not an approximation refined by smaller dt,
  !! so it carries no phase error and serves as the reference against
  !! which the symplectic integrators are measured.
  subroutine advance_analytic(tnow)

    implicit none

    real(8), intent(in) :: tnow
    integer :: i
    real(8) :: om,Qt

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(om,Qt)
    do i=1,Npart
      om = 1.0d0/(j0_part(i)+0.5d0*(Lfix+sqrt(Lfix**2+4.d0)))**3
      Qt = q0_part(i) + om*tnow
      call invert_QJ_to_rp(Qt,j0_part(i),r_part(i),p_part(i))
    end do
    !$OMP END PARALLEL DO

  end subroutine advance_analytic


  subroutine build_cell_list(cell_start,particle_order)

    use omp_lib

    implicit none

    integer, allocatable, intent(out) :: cell_start(:)
    integer, allocatable, intent(out) :: particle_order(:)

    integer :: j,c,th,nth
    integer, allocatable :: cell_count(:),ic(:)
    integer, allocatable :: local_count(:,:),local_offset(:,:)

    allocate(cell_start(1:Nr+1))
    allocate(particle_order(1:Npart))
    allocate(cell_count(1:Nr))
    allocate(ic(1:Npart))

! Counting sort that groups the particles by radial cell, done in two
! parallel passes over per-thread histograms:
!
!   Pass 1: each particle's cell index ic(j) is computed and tallied
!   into local_count(:,th), the private column of the thread that owns
!   it, so no atomics are needed.
!
!   Between passes (serial, O(Nr*nth)): the columns are reduced into the
!   global cell_start, and each thread gets its own write offset per
!   cell, local_offset(c,th) = cell_start(c) + what threads 0..th-1 put
!   in cell c. This splits every cell's slot range into disjoint
!   per-thread sub-ranges.
!
!   Pass 2: each thread writes its own particles into its own sub-range.
!   Both loops use SCHEDULE(STATIC) with no chunk size, which OpenMP
!   guarantees to partition identically for the same bounds and team
!   size, so a particle is handled by the thread that counted it.
!
! Memory layout: the histograms are indexed (cell,thread) rather than
! (thread,cell). Fortran is column-major, so each thread then owns one
! contiguous block of Nr elements; with the transposed layout, entries
! of different threads for the same cell would sit a few elements apart
! and share cache lines, causing false sharing.
    nth = omp_get_max_threads()
    allocate(local_count(1:Nr,0:nth-1))
    allocate(local_offset(1:Nr,0:nth-1))
    local_count = 0

    !$OMP PARALLEL PRIVATE(j,th) SHARED(ic,local_count)
    th = omp_get_thread_num()
    !$OMP DO SCHEDULE(STATIC)
    do j=1,Npart
      ic(j) = nint((r_part(j)-r(1))/dr) + 1
      ic(j) = max(1,min(Nr,ic(j)))
      local_count(ic(j),th) = local_count(ic(j),th) + 1
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    do c=1,Nr
      cell_count(c) = sum(local_count(c,:))
    end do

    cell_start(1) = 1
    do c=1,Nr
      cell_start(c+1) = cell_start(c) + cell_count(c)
    end do

    do c=1,Nr
      local_offset(c,0) = cell_start(c)
      do th=1,nth-1
        local_offset(c,th) = local_offset(c,th-1) + local_count(c,th-1)
      end do
    end do

    !$OMP PARALLEL PRIVATE(j,th,c) SHARED(ic,local_offset,particle_order)
    th = omp_get_thread_num()
    !$OMP DO SCHEDULE(STATIC)
    do j=1,Npart
      c = ic(j)
      particle_order(local_offset(c,th)) = j
      local_offset(c,th) = local_offset(c,th) + 1
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(cell_count,ic,local_count,local_offset)

  end subroutine build_cell_list

  !> Set the time step.
  !! Here we find the time step using information from the
  !! Courant factor and the maximum value of the momentum.
  subroutine set_timestep

! *********************
! ***   TIME STEP   ***
! *********************

! Here we find the time step using information from the
! Courant factor and the maximum value of the momentum.

! Make sure that the Courant condition is satisfied in
! the r direction.

  integer i
  real(8) dtr,dtp       ! Auxiliary variables.
!  real(8) pmax_aux,Fmax_aux
  
!  pmax_aux = MAXVAL(abs(p_part))
!  dtr = courant*dr/pmax_aux
  dtr = courant*dr/pmax
! Second bound, from the acceleration: the time a particle needs to
! cross the radial resolution scale drc under the largest force present,
!
!   dtp = courant * sqrt(2*drc/Fmax),
!
! the usual criterion for symplectic N-body integrators. It tracks the
! curvature of the force field, since leapfrog applied to an oscillator
! is stable for dt*omega <~ 2 with omega^2 ~ dF/dr.
!
! A force can come either from the background or from self-gravity, so
! both cases enable the criterion.

  if (BGtype /= "null" .or. autointeraction) then
    Fmax = 0.0d0

    do i=1,Npart
       Fmax = max(Fmax,abs(force_part(i)))
    end do

    if (Fmax>0.0d0) then
       dtp = courant*sqrt(2.0d0*drc/Fmax)
       dt  = min(dtr,dtp)
    else
       dt = dtr
    end if
  else
    dt = dtr
  end if

  end subroutine set_timestep



  !> Save all the data to the corresponding files
  subroutine save_data

    character(100) :: filename     !< Name of output file

!   Save distribution function f.
    filename = 'vlasov_fdist'
    call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,f)

!   Save rho, curr and cont (multiplied by r**2).
!
!   NOTE: r, rho, avg_rho and curr are allocated with ghost points
!   (bounds 1-ghost:Nr, see alloc_mem_set0), but save1Ddata's dummy
!   arguments are explicit-shape (1:Nr). Passing the whole arrays (or
!   an expression built from them) here associates them by sequence,
!   silently shifting the data by "ghost" points: the output would
!   start with the unphysical ghost points and drop the last "ghost"
!   physical points near r=rmax. Slicing to (1:Nr) selects exactly
!   the physical points and matches the dummy's shape.
    filename = 'vlasov_density'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*rho(1:Nr))
    filename = 'vlasov_avg_density'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*avg_rho(1:Nr))
    filename = 'vlasov_curr'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*curr(1:Nr))
    filename = 'vlasov_energy'
    call save0Ddata(directory,filename,t,total_energy)
    filename = 'vlasov_k_phi_e'
    call save_energy(directory,filename,t,kinetic,potential,total_energy)
!    filename = 'vlasov_cont'
!    call save1Ddata(directory,filename,Nr,Np,t,r,r**2*cont)


!   Save force and potential.
    if (autointeraction) then

       filename = 'vlasov_force'
       call save1Ddata(directory,filename,Nr,t,r(1:Nr),force(1:Nr))

       filename = 'vlasov_potential'
       call save1Ddata(directory,filename,Nr,t,r(1:Nr),pot(1:Nr))

       filename = 'vlasov_potential_r0'
       call save0Ddata(directory,filename,t,pot(1))

       filename = 'vlasov_force_r0'
       call save0Ddata(directory,filename,t,force(1))
    end if

  end subroutine save_data


  subroutine save0Ddata(directory,filename,t,var)

! **********************
! ***   SAVE0DDATA   ***
! **********************

! This subroutine saves 1D data to files.

  implicit none

  real(8) t,var

  character(20) directory,filename,filestatus


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
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(2ES16.8)") t,var


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save0Ddata



  subroutine save1Ddata(directory,filename,Nr,t,r,var)

! **********************
! ***   SAVE1DDATA   ***
! **********************

! This subroutine saves 1D data to files.

  implicit none

  integer i
  integer Nr

  real(8) t

  real(8), dimension(1:Nr) :: r,var

  character(20) directory,filename,filestatus


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
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=1,Nr
     if (dabs(var(i)).gt.1.0D-50) then
        write(101,"(2ES16.8)") r(i),var(i)

     else
        write(101,"(2ES16.8)") r(i),0.0d0
     end if
  end do

! Leave two blank spaces before next time.
! The reason to leave two spaces is that 'gnuplot' asks
! for two spaces to distinguish different records.

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save1Ddata



  subroutine save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,var)

! **********************
! ***   SAVE2DDATA   ***
! **********************

! This subroutine saves 2D data to files.

  implicit none

  integer i
  integer Npart

  real(8) t

  real(8), dimension(1:Npart) :: r_part
  real(8), dimension(1:Npart) :: p_part

  real(8), dimension(1:Npart) :: var

  character(20) directory,filename,filestatus


! ***************************
! ***   OPEN DATA FILES   ***
! ***************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open files.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus)

  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus,position='append')
  end if


! ************************
! ***   SAVE 2D DATA   ***
! ************************

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=1,Npart
    if (dabs(var(i)).gt.1.0D-50) then
      write(101,"(3ES16.8)") r_part(i),p_part(i),var(i)
    else
      write(101,"(3ES16.8)") r_part(i),p_part(i),0.0D0
    end if
  end do

  write (101,*)
  write (101,*)


! ****************************
! ***   CLOSE DATA FILES   ***
! ****************************

  close(101)

! ***************
! ***   END   ***
! ***************

  end subroutine save2Ddata_particles


  subroutine save_energy(directory,filename,t,kinetic,potential,energy)

! **********************
! ***   SAVE ENERGY  ***
! **********************

! This subroutine saves the energy (kinectic,potential,virial) data to files.

  implicit none

  real(8) t,kinetic,potential,energy

  character(20) directory,filename,filestatus


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
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(5ES16.8)") t,kinetic,potential,energy


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_energy



  subroutine save_dens_current(directory,filename,Nr,Np,t,r,density,current,error)

! ************************************************
! ***   SAVE DENSITY, CURRENT AND ERROR DATA   ***
! ************************************************

! This subroutine saves the density, current and error in the
! continuity equation.

  implicit none

  integer i
  integer Nr,Np

  real(8) t

  real(8), dimension(0:Nr) :: r,density,current,error

  character(20) directory,filename,filestatus


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
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  density = density + 1.0D-50
  current = current + 1.0D-50
  error   = error   + 1.0D-50

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(101,"(4ES16.8)") r(i),density(i),current(i),error(i)
  end do

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_dens_current


  subroutine save_force_pot(directory,filename,Nr,Np,t,r,force,pot)

! *****************************************
! ***   SAVE FORCE AND POTENTIAL DATA   ***
! *****************************************

! This subroutine saves the force and potential data for the
! self gravitating case

  implicit none

  integer i
  integer Nr,Np

  real(8) t

  real(8), dimension(0:Nr) :: r,force,pot

  character(20) directory,filename,filestatus


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
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  force = force + 1.0D-50
  pot   = pot   + 1.0D-50

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(101,"(3ES16.8)") r(i),force(i),pot(i)
  end do

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_force_pot


subroutine reduce_arrays
  use parameters
  use arrays

 
  integer i,j
  integer counter,Npart_aux
  real(8), allocatable, dimension (:) :: r_aux,p_aux,f_aux
!  real(8), dimension (1:Npart) :: r_aux,p_aux,f_aux


! Count How many particles are still in the grid

  counter = 0

  do i=1,Npart
    if (r_part(i)<= rmax) then
      counter = counter + 1
    end if 
  end do

! Reduce the size of arrays if we have less than 90% 
! of the original number of particles.

  if (dble(counter) < 0.9D0* Npart) then 
  allocate(r_aux(1:Npart))
  allocate(p_aux(1:Npart))
  allocate(f_aux(1:Npart))

! Copy the position and momentum of particles

  r_aux = r_part
  p_aux = p_part
  f_aux = f

! Deallocate arrays

  deallocate(r_part)
  deallocate(r_part_p)
  deallocate(p_part)
  deallocate(p_part_p)
  deallocate(p_part_h)
  deallocate(pot_part)
  deallocate(force_part)
  deallocate(f)
!  deallocate(f_p)

! Allocate arrays 

  Npart_aux = Npart
  Npart     = counter

  allocate(r_part(1:Npart))
  allocate(r_part_p(1:Npart))
  allocate(p_part(1:Npart))
  allocate(p_part_p(1:Npart))
  allocate(p_part_h(1:Npart))
  allocate(pot_part(1:Npart))
  allocate(force_part(1:Npart))
  allocate(f(1:Npart))
!  allocate(f_p(1:Npart))

  j = 1
  do i=1,Npart_aux
    if (r_aux(i)<=rmax) then
      r_part(j) = r_aux(i) 
      p_part(j) = p_aux(i)
      f(j)      = f_aux(i)
      j = j+1
    end if

  end do
  print *, "In the grid are still", Npart, "computational particles"
  end if
end subroutine reduce_arrays

end module utils

