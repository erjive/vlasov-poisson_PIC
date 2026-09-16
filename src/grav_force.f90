
! **************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE   ***
! **************************************************

! Background field case.  Here we assume that the background
! gravitational field corresponds to the case of a constant
! density star of total mass 1 and radius 1, which has a
! gravitational potential "pot" given by:
!
! pot  =  1/2 ( r**2 - 3 )     r <  1
!
! pot  =  - 1 / r              r >= 1  
!
! for which the force is (force = - dpot/dr):
!
! force  = - r                 r <  1
! 
! force  = - 1 / r**2          r >= 1

subroutine grav_force



  use parameters
  use arrays
  use utils
  implicit none
  integer i
  real(8) :: smallpi
  real(8) :: sq,den       ! per-particle temporaries for the fused loops
                          ! (NB: "r2" would clash with parameters::r2)
  character(100) :: filename
  smallpi = acos(-1.0d0)

! Self-gravitating case.  In this case we need to
! solve the Poisson equation.

  if (autointeraction) then

     if (rmin/=0.0d0) then
        print *
        print *, 'For the self-gravitating case you should have rmin=0.'
        print *, 'Aborting ...'
        print *
        stop
     else
        call poisson_rk
     end if

  end if 

  if (forcetype=="bg") then

     if (BGtype == "null") then
     
        pot = pot 
        force = force
        pot_part = pot_part
        force_part = force_part

     else if (BGtype == "sphere") then

       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do i=1,Npart
            if (r_part(i)<1.d0) then
               pot_part(i)   =  0.5d0*(r_part(i)**2 - 3.d0)
               force_part(i) = - r_part(i)
            else
               pot_part(i)   = - 1.d0/r_part(i)
               force_part(i) = - 1.d0/r_part(i)**2
            end if
       end do
       !$OMP END PARALLEL DO

     else if (BGtype == "Isochrone") then

!      Fused into a single pass over the particles, with sqrt(1+r^2)
!      evaluated ONCE. Written as whole-array expressions, this was two
!      passes with the same sqrt recomputed 2-3 times per particle (and
!      two more passes in the angular-momentum block below). Those passes
!      were serial and, once analysish stopped dominating the run cost,
!      became the bottleneck: the measured OpenMP parallel fraction had
!      dropped to p~0.11. See BUGS_TODO.md.

       if (autointeraction) then

         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sq)
         do i=1,Npart
           sq = sqrt(1.0D0+r_part(i)**2)
           pot_part(i)   = pot_part(i)   - 1.0D0/(1.0D0+sq)
           force_part(i) = force_part(i) - r_part(i)/(sq*(1.0D0+sq)**2)
         end do
         !$OMP END PARALLEL DO

         pot = pot + (-1.0D0/(1.0D0+sqrt(1.0D0+r**2)))
         force = force + (-r/(sqrt(1.D0+r**2)*(1.D0+sqrt(1.D0+r**2))**2))

       else

         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sq)
         do i=1,Npart
           sq = sqrt(1.0D0+r_part(i)**2)
           pot_part(i)   = -1.0D0/(1.0D0+sq)
           force_part(i) = -r_part(i)/sq*pot_part(i)**2
         end do
         !$OMP END PARALLEL DO

       end if

     else if (BGtype == "iso") then

        pot   = 3.0d0*log(r)
        force =-3.0d0/r

     else if (BGtype == "isotrun") then

        pot   = (10.0d0/6.0d0)*( 2.0d0*atan(r)/r + log(r**2+1) )
        force = -(10.0d0/3.0d0)*( r-atan(r) )/r**2

     else if (BGtype == "nfw") then

        pot   = -16.0d0*log( 1.0d0+r )/r
        force = -16.0d0*( log(1.0d0+r)-r/(1.0d0+r) )/r**2

     else if (BGtype == "burkert") then

        pot =    ( 10.0d0/(3.0d0*r) )*( 2.0d0*(1.0d0+r)*atan(r) -2.0d0*(1.0d0+r)*log(1.0d0+r) -(1.0d0-r)*log(1.0d0+r**2) )
        force = -( 10.0d0/(3.0d0*r*r) )*( log( (1.0d0+r**2)*(1.0d0+r)**2 ) - 2.0d0*atan(r) )

     else

       print *
       print *, 'Unknown type of gravitational force'
       print *, 'Aborting ...'
       print *
       stop

     end if 



  end if
  

! *******************************
! ***   ADD ANGULAR MOMENTUM  ***
! *******************************
     if(Lfix /= 0.0d0) then

!       Also fused: one pass, with the denominator formed once instead of
!       twice per particle.
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
        do i=1,Npart
          den = r_part(i)**2 + eps*eps
          pot_part(i)   = pot_part(i)   + 0.5d0*Lfix**2/den
          force_part(i) = force_part(i) + Lfix**2*r_part(i)/den**2
        end do
        !$OMP END PARALLEL DO

        !pot   = pot + 0.5d0*Lfix**2/(r**2 + eps*eps)
        !force = force + Lfix**2*r/(r**2 + eps*eps)**2

     end if


  !filename = 'vlasov_potpart'
  !call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,pot_part)

  !filename = 'vlasov_comparepotential'
  !call save1Ddata(directory,filename,Nr,t,r,pot)

end subroutine grav_force
