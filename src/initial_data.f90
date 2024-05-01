! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays
    use utils

    implicit none

    logical :: accepted
    integer :: i,j
    real(8) :: smallpi,f_max
    real(8) :: raux,paux
    real(8) :: rand3(3)
    real(8) :: gaussian_fixedL
    real(8) :: w,x,y,z
    real(8) :: energy 

!   Auxiliary variables for a distribution function 
!   that depends on action-angle varialbes
    real(8) :: Jr, Qr, s, s1, s2, er1, er2, eta,argaux

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
          raux = rminc+(dble(i)+0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = gaussian_fixedL(1.0D0,r0,p0,Lfix,raux,paux,sr,sp)
        end do
      end do
      !$OMP END PARALLEL DO

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc
    
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

    else if(state .eq."aa") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux     

          energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
          er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          argaux = (s1+s2-2.0*s)/(s2-s1)

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - sqrt((-2.d0*energy)**3)*sqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*sin(eta)
          Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(Lfix+sqrt(Lfix**2+4.d0))

          f((i-1)*Npc+j) = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

          if ((f((i-1)*Npc+j) /= f((i-1)*Npc+j) )) then

            f((i-1)*Npc+j) = 0.D0
            r_part((i-1)*Npc+j) = 10000.D0
          end if

        end do
      end do
      !$OMP END PARALLEL DO
      
      f_max = maxval(f)
      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
           if ((r0-0.00)*f_max<=f((i-1)*Npc+j) .and. f((i-1)*Npc+j)<= r0*f_max ) then
              f((i-1)*Npc+j)=0.0D0
           else if (f((i-1)*Npc+j)<= (r0-0.00)*f_max ) then
              r_part((i-1)*Npc+j) = 100000.D0
           end if
!          !print *, Qr,Jr
        end do
      end do 
      !$OMP END PARALLEL DO
      call reduce_arrays

      print *, a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))
      f = a0/(drc*dpc*8.0*smallpi**2*Lfix*sum(f))*f
      !f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)*8.0*smallpi**2*Lfix*drc*dpc

    else if(state .eq."aa_random") then

      Npart = Nrc*Npc
      f_max = dexp(-dsin(0.5d0*0.0)**2/sp**2)*dexp(-sr**2/sr**2)*sr**2

      i = 1
      do while(i<= Npart)
!      do i=1,Npart
!        accepted = .false.
!        do while (.not. accepted)
! The subroutine random_number creates a pseudo-random number in the interval (0,1],
! we would like to include the zero, so here I just create a triad of three uniform 
! random numbers in the interval [0,1]

        call random_number(rand3)

        x = 1.d0-rand3(1)
        y = 1.d0-rand3(2)
        z = 1.d0-rand3(3)
  
! Random position and momentum of the particle.

        raux = (rmaxc-rminc)*x + rminc
        paux = (pmaxc-pminc)*y + pminc
        z = f_max*z 

        energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*Lfix**2/(raux**2) + 0.5D0*paux**2
        if (energy < 0.0D0) then

          er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
          s1 = 1.d0 + sqrt(1.d0+er1**2)
          s2 = 1.d0 + sqrt(1.d0+er2**2)
          s  = 1.d0 + sqrt(1.d0+raux**2)
          !eta = dacos(sign(min(abs(2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)),1.0),2.d0/(s1-s2)*(s-(s1+s2)*0.5d0)))
          argaux = (s1+s2-2.0*s)/(s2-s1)

          if (paux>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
          end if

          Qr = eta - dsqrt((-2.d0*energy)**3)*dsqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*dsin(eta)
          Jr = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))
           w = dexp(-dsin(0.5d0*Qr)**2/sp**2)*dexp(-Jr**2/sr**2)*Jr**2

          if (z <= w ) then
        
            r_part(i) = raux
            p_part(i) = paux
            f(i)      = w
            !accepted = .true.
            i = i+1
          end if
        end if
        !end do
      end do
      



!      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
!      !$OMP END PARALLEL DO
      print *, a0/sum(f)
      f = a0/sum(f)*f
      !f = f*drc*dpc*8.D0*smallpi**2
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


