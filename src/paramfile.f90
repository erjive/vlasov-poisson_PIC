! ===========================================================================
! paramfile.f90
! ===========================================================================
!> Reads the run parameters from a "name = value" file, plus optional
!! overrides given on the command line:
!!
!!     ./VP_PIC [file] [name=value] [name=value] ...
!!
!! The file defaults to "input_parameters" in the current directory. Every
!! parameter keeps the default declared in parameters.f90 unless the file or
!! the command line sets it, so a run file only needs to list what differs
!! from those defaults, in any order. Names are case-insensitive, "!" and "#"
!! start a comment, and blank lines are ignored.
!!
!! Overrides are applied after the file, so a parameter sweep can reuse one
!! base file and vary only the parameters that change from run to run.

module paramfile

  use parameters

  implicit none

  private
  public :: read_parameters, parameter_file

  character(256) :: parameter_file = 'input_parameters'  !< File actually read

!> Every name the parser accepts. Used to reject unknown names and to
!! suggest the intended one when a name is misspelled, so it must list
!! exactly the names handled by assign_param below.
  integer, parameter :: NPARAM = 41
  character(16), parameter :: pname(NPARAM) = [ character(16) :: &
      'dr', 'Nrc', 'Npc', 'courant', 'Nt',                       &
      'rmin', 'rmax', 'rminc', 'rmaxc', 'pminc', 'pmaxc', 'pmax',&
      'Lfix', 'eps', 'reduceparticles', 'Nreduce',               &
      'time_output', 'spatial_output', 'field_output',           &
      'directory', 'output_format',                              &
      'a0', 'r0', 'p0', 'sr', 'sp', 'state',                     &
      'j1', 'j2', 'sj1', 'sj2', 'sq1', 'sq2',                    &
      'r1', 'r2',                                                &
      'bsplineorder', 'integrator', 'spatialorder',              &
      'forcetype', 'BGtype', 'autointeraction' ]

 contains

! ***************************************
! ***   READ FILE AND COMMAND LINE    ***
! ***************************************

  subroutine read_parameters

    implicit none

    integer :: nargs,iarg,u,ios,lineno,eq
    character(512) :: line,arg
    character(300) :: origin      ! file name plus line number, long enough for a full path
    logical :: thereis

    nargs = command_argument_count()

    if (nargs >= 1) then
       call get_command_argument(1,line)
       if (line(1:1) == '-') call usage_and_stop(trim(line))
!      A first argument holding "=" is an override, not a file name, so the
!      default file is kept and every argument is treated as an override.
       if (index(line,'=') == 0) then
          parameter_file = line
          iarg = 2
       else
          iarg = 1
       end if
    else
       iarg = 1
    end if

    inquire(file=trim(parameter_file),exist=thereis)

    if (.not.thereis) then
       print *
       print *, 'Parameter file not found: ',trim(parameter_file)
       call usage_and_stop('')
    end if

    open(newunit=u,file=trim(parameter_file),status='old',action='read')

    lineno = 0

    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       lineno = lineno + 1

       call strip_comment(line)
       if (len_trim(line) == 0) cycle

       eq = index(line,'=')
       write(origin,'(a,i0)') trim(parameter_file)//', line ',lineno

       if (eq == 0) then
          print *
          print *, trim(origin),': expected "name = value", found:'
          print *, '   ',trim(adjustl(line))
          call stop_run
       end if

       call assign_param(line(1:eq-1),line(eq+1:),origin)
    end do

    close(u)

!   Command line overrides, applied last so they win over the file.

    do while (iarg <= nargs)
       call get_command_argument(iarg,arg)
       if (arg(1:1) == '-') call usage_and_stop(trim(arg))
       eq = index(arg,'=')
       if (eq == 0) then
          print *
          print *, 'Command line arguments after the file must be name=value, found:'
          print *, '   ',trim(arg)
          call usage_and_stop('')
       end if
       call assign_param(arg(1:eq-1),arg(eq+1:),'command line')
       iarg = iarg + 1
    end do

    call report(nargs)
    call validate

  end subroutine read_parameters


! *********************************
! ***   ASSIGN ONE PARAMETER    ***
! *********************************

!> Store one "name = value" pair. "origin" names the file and line (or the
!! command line) the pair came from, so any error can point at it.

  subroutine assign_param(rawname,rawvalue,origin)

    implicit none

    character(*), intent(in) :: rawname,rawvalue,origin

    character(64)  :: name
    character(256) :: value

    name  = adjustl(rawname)
    value = adjustl(rawvalue)
    call strip_quotes(value)

    if (len_trim(value) == 0) then
       print *
       print *, trim(origin),': parameter "',trim(name),'" has no value.'
       call stop_run
    end if

    select case (lower(name))

!   Grid.
    case ('dr')              ; call get_real(value,dr,name,origin)
    case ('nrc')             ; call get_int (value,Nrc,name,origin)
    case ('npc')             ; call get_int (value,Npc,name,origin)
    case ('rmin')            ; call get_real(value,rmin,name,origin)
    case ('rmax')            ; call get_real(value,rmax,name,origin)
    case ('rminc')           ; call get_real(value,rminc,name,origin)
    case ('rmaxc')           ; call get_real(value,rmaxc,name,origin)
    case ('pminc')           ; call get_real(value,pminc,name,origin)
    case ('pmaxc')           ; call get_real(value,pmaxc,name,origin)
    case ('pmax')            ; call get_real(value,pmax,name,origin)

!   Time.
    case ('courant')         ; call get_real(value,courant,name,origin)
    case ('nt')              ; call get_int (value,Nt,name,origin)

!   Physical setup.
    case ('lfix')            ; call get_real(value,Lfix,name,origin)
    case ('eps')             ; call get_real(value,eps,name,origin)

!   Particle bookkeeping.
    case ('reduceparticles') ; call get_log (value,reduceparticles,name,origin)
    case ('nreduce')         ; call get_int (value,Nreduce,name,origin)

!   Output.
    case ('time_output')     ; call get_int (value,time_output,name,origin)
    case ('spatial_output')  ; call get_int (value,spatial_output,name,origin)
    case ('field_output')    ; call get_int (value,field_output,name,origin)
    case ('directory')       ; call get_str (value,directory,name,origin)
    case ('output_format')   ; call get_str (value,output_format,name,origin)

!   Initial distribution.
    case ('a0')              ; call get_real(value,a0,name,origin)
    case ('r0')              ; call get_real(value,r0,name,origin)
    case ('p0')              ; call get_real(value,p0,name,origin)
    case ('sr')              ; call get_real(value,sr,name,origin)
    case ('sp')              ; call get_real(value,sp,name,origin)
    case ('state')           ; call get_str (value,state,name,origin)

!   Test functions for the h_k modes.
    case ('j1')              ; call get_real(value,j1,name,origin)
    case ('j2')              ; call get_real(value,j2,name,origin)
    case ('sj1')             ; call get_real(value,sj1,name,origin)
    case ('sj2')             ; call get_real(value,sj2,name,origin)
    case ('sq1')             ; call get_real(value,sq1,name,origin)
    case ('sq2')             ; call get_real(value,sq2,name,origin)

!   Radial window of the averaged density.
    case ('r1')              ; call get_real(value,r1,name,origin)
    case ('r2')              ; call get_real(value,r2,name,origin)

!   Numerical methods.
    case ('bsplineorder')    ; call get_int (value,bsplineorder,name,origin)
    case ('integrator')      ; call get_str (value,integrator,name,origin)
    case ('spatialorder')    ; call get_str (value,spatialorder,name,origin)
    case ('forcetype')       ; call get_str (value,forcetype,name,origin)
    case ('bgtype')          ; call get_str (value,BGtype,name,origin)
    case ('autointeraction') ; call get_log (value,autointeraction,name,origin)

    case default
       call unknown_name(name,origin)

    end select

  end subroutine assign_param


!> Report an unrecognised name, with the closest known name as a hint.

  subroutine unknown_name(name,origin)

    implicit none

    character(*), intent(in) :: name,origin
    integer :: i,d,dbest,ibest

    dbest = huge(1)
    ibest = 1

    do i=1,NPARAM
       d = edit_distance(lower(name),lower(pname(i)))
       if (d < dbest) then
          dbest = d
          ibest = i
       end if
    end do

    print *
    print *, trim(origin),': unknown parameter "',trim(name),'".'

!   A small edit distance means a typo rather than a different name.
    if (dbest <= 3) print *, '   Did you mean "',trim(pname(ibest)),'"?'

    print *, '   Run with --help for the full list.'
    call stop_run

  end subroutine unknown_name


! *****************************
! ***   VALUE CONVERSION    ***
! *****************************

  subroutine get_real(value,var,name,origin)

    implicit none

    character(*), intent(in)  :: value,name,origin
    real(8),      intent(out) :: var
    integer :: ios

    read(value,*,iostat=ios) var
    if (ios /= 0) call bad_value(value,name,origin,'a real number')

  end subroutine get_real


  subroutine get_int(value,var,name,origin)

    implicit none

    character(*), intent(in)  :: value,name,origin
    integer,      intent(out) :: var
    integer :: ios

    read(value,*,iostat=ios) var
    if (ios /= 0) call bad_value(value,name,origin,'an integer')

  end subroutine get_int


!> Accepts the Fortran spellings plus the ones a user is likely to type.

  subroutine get_log(value,var,name,origin)

    implicit none

    character(*), intent(in)  :: value,name,origin
    logical,      intent(out) :: var

    select case (lower(value))
    case ('.true.','true','t','yes','y','on','1')
       var = .true.
    case ('.false.','false','f','no','n','off','0')
       var = .false.
    case default
       call bad_value(value,name,origin,'true or false')
    end select

  end subroutine get_log


!> Character parameters are fixed length, so a value that does not fit is
!! rejected instead of being silently truncated.

  subroutine get_str(value,var,name,origin)

    implicit none

    character(*), intent(in)  :: value,name,origin
    character(*), intent(out) :: var

    if (len_trim(value) > len(var)) then
       print *
       print *, trim(origin),': value of "',trim(name),'" is too long'
       print '(a,i0,a)', '    (maximum ',len(var),' characters).'
       call stop_run
    end if

    var = trim(value)

  end subroutine get_str


  subroutine bad_value(value,name,origin,expected)

    implicit none

    character(*), intent(in) :: value,name,origin,expected

    print *
    print *, trim(origin),': "',trim(name),'" expects ',trim(expected),', found "',trim(value),'".'
    call stop_run

  end subroutine bad_value


! ************************
! ***   VALIDATION     ***
! ************************

  subroutine validate

    implicit none

    call check_option(output_format,'output_format','ascii hdf5 raw')
    call check_option(state,'state','gaussian aa aa_halton aa_quad aa_random checkpoint')
    call check_option(integrator,'integrator','euler leapfrog yoshida4 yoshida6 analytic rk4')
    call check_option(spatialorder,'spatialorder','two four')
    call check_option(forcetype,'forcetype','bg self')
    call check_option(BGtype,'BGtype','null sphere Isochrone iso isotrun nfw burkert')

    if (rmin < 0.0d0)  call fail('rmin must be greater than or equal to zero.')
    if (rmax <= rmin)  call fail('rmax must be greater than rmin.')
    if (rmaxc <= rminc) call fail('rmaxc must be greater than rminc.')
    if (pmaxc <= pminc) call fail('pmaxc must be greater than pminc.')
    if (dr <= 0.0d0)   call fail('dr must be positive.')
    if (Nrc <= 0 .or. Npc <= 0) call fail('Nrc and Npc must be positive.')
    if (courant <= 0.0d0) call fail('courant must be positive.')
    if (spatial_output <= 0 .or. time_output <= 0 .or. field_output <= 0) &
       call fail('time_output, spatial_output and field_output must be positive.')

!   A particle snapshot also stores rho and the energies, which are only
!   recomputed every spatial_output steps. Unless field_output is a
!   multiple of spatial_output, a snapshot would carry grid quantities
!   belonging to an earlier time.
    if (mod(field_output,spatial_output) /= 0) then
       print *
       print *, 'field_output must be a multiple of spatial_output,'
       print *, 'because the grid quantities a particle snapshot stores are'
       print *, 'only recomputed every spatial_output steps.'
       print '(a,i0,a,i0)', '    field_output = ',field_output, &
                            ',  spatial_output = ',spatial_output
       call stop_run
    end if

!   eps softens the centrifugal term of the Hamiltonian that moves the
!   particles, while the action-angle diagnostics reconstruct (E,J3,Q3)
!   from the unsoftened one. With eps /= 0 the particles are therefore
!   evolved and analysed under different Hamiltonians, J3 is no longer
!   conserved and h_k floors out. For L0 /= 0 the barrier already keeps r
!   away from the origin, so the softening is only useful as L0 -> 0.
    if (eps /= 0.0d0 .and. Lfix /= 0.0d0) then
       print *
       print *, 'WARNING: eps /= 0 with Lfix /= 0.'
       print *, 'The dynamics will use a softened centrifugal term while'
       print *, 'analysish reconstructs (E,J3,Q3) unsoftened: J3 will drift'
       print *, 'and h_k will floor out. Set eps = 0 unless you know why.'
       print *, 'eps =',eps
       print *
    end if

  end subroutine validate


!> Check a parameter against a blank-separated list of valid values, and
!! print that list if it does not match.

  subroutine check_option(value,name,options)

    implicit none

    character(*), intent(in) :: value,name,options

    if (index(' '//lower(options)//' ',' '//trim(lower(value))//' ') == 0) then
       print *
       print *, 'Unknown ',trim(name),': "',trim(value),'".'
       print *, '   Valid values: ',trim(options)
       call stop_run
    end if

  end subroutine check_option


  subroutine fail(message)

    implicit none

    character(*), intent(in) :: message

    print *
    print *, trim(message)
    call stop_run

  end subroutine fail


! ***********************
! ***   REPORTING     ***
! ***********************

!> Echo what was read, so the run log records the file and any override.

  subroutine report(nargs)

    implicit none

    integer, intent(in) :: nargs
    integer :: i,eq
    character(512) :: arg
    logical :: any

    print *
    print *, 'Parameters read from: ',trim(parameter_file)

    any = .false.

    do i=1,nargs
       call get_command_argument(i,arg)
       eq = index(arg,'=')
       if (eq == 0) cycle
       if (.not.any) then
          print *, 'Command line overrides:'
          any = .true.
       end if
       print *, '   ',trim(arg)
    end do

  end subroutine report


  subroutine usage_and_stop(arg)

    implicit none

    character(*), intent(in) :: arg
    integer :: i,j

    if (arg == '--help' .or. arg == '-h') then
       print *
       print *, 'Usage: VP_PIC [file] [name=value] ...'
       print *
       print *, '  file        parameter file (default "input_parameters").'
       print *, '  name=value  override a parameter read from the file.'
       print *
       print *, 'Parameter names:'
       do i=1,NPARAM,4
          print '(4(3x,a16))', (pname(j),j=i,min(i+3,NPARAM))
       end do
       print *
       print *, 'Anything not listed in the file keeps the default in parameters.f90.'
       print *
       stop
    end if

    print *
    print *, 'Usage: VP_PIC [file] [name=value] ...'
    if (len_trim(arg) > 0) print *, 'Unrecognised option: ',trim(arg)
    print *, 'Run with --help for the list of parameter names.'
    print *
    stop

  end subroutine usage_and_stop


  subroutine stop_run

    implicit none

    print *, 'Aborting ...'
    print *
    stop

  end subroutine stop_run


! **********************
! ***   UTILITIES    ***
! **********************

!> Drop everything from the first comment marker onwards.

  subroutine strip_comment(line)

    implicit none

    character(*), intent(inout) :: line
    integer :: i

    i = scan(line,'!#')
    if (i > 0) line = line(1:i-1)

  end subroutine strip_comment


!> Remove one pair of surrounding quotes, so quoted values are accepted too.

  subroutine strip_quotes(value)

    implicit none

    character(*), intent(inout) :: value
    integer :: n

    n = len_trim(value)
    if (n < 2) return

    if ((value(1:1) == '"'  .and. value(n:n) == '"' ) .or. &
        (value(1:1) == "'"  .and. value(n:n) == "'" )) then
       value = value(2:n-1)
    end if

  end subroutine strip_quotes


  function lower(s) result(out)

    implicit none

    character(*), intent(in) :: s
    character(len(s)) :: out
    integer :: i,c

    out = s

    do i=1,len(s)
       c = iachar(s(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c+32)
    end do

  end function lower


!> Levenshtein distance, used only to suggest a name after a typo. Two rows
!! of the matrix are enough, since each row depends only on the previous.

  function edit_distance(a,b) result(d)

    implicit none

    character(*), intent(in) :: a,b
    integer :: d,la,lb,i,j,cost
    integer, allocatable :: prev(:),cur(:)

    la = len_trim(a)
    lb = len_trim(b)

    allocate(prev(0:lb),cur(0:lb))

    do j=0,lb
       prev(j) = j
    end do

    do i=1,la
       cur(0) = i
       do j=1,lb
          cost = 1
          if (a(i:i) == b(j:j)) cost = 0
          cur(j) = min(cur(j-1)+1, prev(j)+1, prev(j-1)+cost)
       end do
       prev = cur
    end do

    d = prev(lb)

  end function edit_distance


end module paramfile
