! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays

    implicit none

    integer :: i,j
    real(8) :: smallpi,f_max
    real(8) :: raux,paux
    real(8) :: rand3(3)
    real(8) :: gaussian_fixedL
    real(8) :: w,x,y,z
    real(8) :: energy 

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
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = gaussian_fixedL(a0,r0,p0,Lfix,raux,paux,sr,sp)
        end do
      end do
      !$OMP END PARALLEL DO

      
          f = f*drc*dpc*8.D0*smallpi**2*Lfix
    
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


