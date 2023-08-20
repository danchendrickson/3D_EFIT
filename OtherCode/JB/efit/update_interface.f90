
Subroutine interface_velocity(lsecxStart,lsecyStart,  &
     &    lsecxEnd,lsecyEnd,zlevel,lx,ly,numz,dtodsp,dden,rden,  &
     &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
                    
  Implicit none

  Integer, intent(in) :: lsecxStart,lsecyStart,lsecxEnd,lsecyEnd,zlevel
  Integer, intent(in) :: lx,ly,numz,xlb,xub,ylb,yub,zlb,zub
  Double Precision, intent(in) :: dtodsp,dden,rden
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)  &
            &    :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)  &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz

  Integer :: ix,iy,tix,tiy,iz

  Integer :: sxi,syi,szi

  if (xlb == 1) then
     sxi = lsecxStart + 1
  else
     sxi = lsecxStart
  end if
  
  !   go across in x direction
  !   update interior points
  do ix = sxi,lsecxEnd-1
     do tiy = lsecyStart+1,lsecyEnd
           vx(ix,tiy,zlevel) = vx(ix,tiy,zlevel) + dtodsp*( &
                &	Txx(ix+1,tiy,zlevel)-Txx(ix,tiy,zlevel) + &
                &  Txy(ix,tiy,zlevel)-Txy(ix,tiy-1,zlevel) + &
                &  Txz(ix,tiy,zlevel)-Txz(ix,tiy,zlevel-1))
     end do
  end do

  !  update boundary points
  if (xlb == 1) then
     do tiy = lsecyStart,lsecyEnd
           vx(lsecxStart,tiy,zlevel) = vx(lsecxStart,tiy,zlevel) + &
                &   dtodsp*(2*Txx(lsecxStart+1,tiy,zlevel))
     end do
  end if
  if (xub == 1) then
     do tiy = lsecyStart,lsecyEnd
           vx(lsecxEnd,tiy,zlevel) = vx(lsecxEnd,tiy,zlevel) + &
                &                   dtodsp*(-2*Txx(lsecxEnd,tiy,zlevel))
     end do
  end if

  if (ylb == 1) then
     syi = lsecyStart + 1
  else
     syi = lsecyStart
  end if
  !   go across in y direction
  !   update interior points
  do tix = lsecxStart+1,lsecxEnd
     do iy = syi,lsecyEnd-1
           vy(tix,iy,zlevel) = vy(tix,iy,zlevel) + dtodsp*( &
                &	Txy(tix,iy,zlevel)-Txy(tix-1,iy,zlevel) + &
                &	Tyy(tix,iy+1,zlevel)-Tyy(tix,iy,zlevel) + &
                &	Tyz(tix,iy,zlevel)-Tyz(tix,iy,zlevel-1))
     end do
  end do
  !  update boundary points
  if (ylb == 1) then
     do tix = lsecxStart,lsecxEnd
           vy(tix,lsecyStart,zlevel) = vy(tix,lsecyStart,zlevel) + &
                &   dtodsp*(2*Tyy(tix,lsecyStart+1,zlevel))
     end do
  end if
  if (yub == 1) then
     do tix = lsecxStart,lsecxEnd
           vy(tix,lsecyEnd,zlevel) = vy(tix,lsecyEnd,zlevel) + &
                &  dtodsp*(-2*Tyy(tix,lsecyEnd,zlevel)) 
     end do
  end if

  !   go down in z direction
  !   update interior points
!  do tix = lsecxStart+1,lsecxEnd
!     do tiy = lsecyStart+1,lsecyEnd
!        vz(tix,tiy,zlevel-1) = vz(tix,tiy,zlevel-1) + dtodsp*( &
!             &	Txz(tix,tiy,zlevel-1)-Txz(tix-1,tiy,zlevel-1) + &
!             &	Tyz(tix,tiy,zlevel-1)-Tyz(tix,tiy-1,zlevel-1) + &
!             &	Tzz(tix,tiy,zlevel)-Tzz(tix,tiy,zlevel-1))
!     end do
!  end do
  do tix = lsecxStart+1,lsecxEnd
     do tiy = lsecyStart+1,lsecyEnd
           vz(tix,tiy,zlevel-1) = vz(tix,tiy,zlevel-1) +        &
                &   dtodsp*dden*(2/(dden+rden))*( &
                &	Txz(tix,tiy,zlevel-1)-Txz(tix-1,tiy,zlevel-1) + &
                &	Tyz(tix,tiy,zlevel-1)-Tyz(tix,tiy-1,zlevel-1) + &
                &	Tzz(tix,tiy,zlevel)-Tzz(tix,tiy,zlevel-1))
     end do
  end do
end Subroutine interface_velocity



Subroutine interface_stress(lsecxStart,lsecyStart,  &
     &    lsecxEnd,lsecyEnd,zlevel,lx,ly,numz,  &
     &    dtods,l2m,lambda,mu,mu2,xlb,xub,ylb,yub,zlb,zub, &
     &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
                   
  Implicit none

  Integer, intent(in) :: lsecxStart,lsecyStart
  Integer, intent(in) :: lsecxEnd,lsecyEnd,zlevel,lx,ly,numz
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
  Double Precision, intent(in) :: dtods,l2m,lambda,mu,mu2
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(in) :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)      &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz

  Integer :: ix,iy,tix,tiy,iz

  Integer :: sxi,syi,szi

  do tix = lsecxStart+1,lsecxEnd
     do tiy = lsecyStart+1,lsecyEnd
        Txx(tix,tiy,zlevel) = Txx(tix,tiy,zlevel) + dtods*( &
             &	l2m*(vx(tix,tiy,zlevel)-vx(tix-1,tiy,zlevel)) + &
             &	lambda*(vy(tix,tiy,zlevel)-vy(tix,tiy-1,zlevel) + &
             &	vz(tix,tiy,zlevel)-vz(tix,tiy,zlevel-1)))
        Tyy(tix,tiy,zlevel) = Tyy(tix,tiy,zlevel) + dtods*( &
             &	l2m*(vy(tix,tiy,zlevel)-vy(tix,tiy-1,zlevel)) + &
             &	lambda*(vx(tix,tiy,zlevel)-vx(tix-1,tiy,zlevel) + &
             &	vz(tix,tiy,zlevel)-vz(tix,tiy,zlevel-1)))
        Tzz(tix,tiy,zlevel) = Tzz(tix,tiy,zlevel) + dtods*( &
             &	l2m*(vz(tix,tiy,zlevel)-vz(tix,tiy,zlevel-1)) + &
             &	lambda*(vx(tix,tiy,zlevel)-vx(tix-1,tiy,zlevel) + &
             &	vy(tix,tiy,zlevel)-vy(tix,tiy-1,zlevel)))
     end do
  end do
  
  if (xlb == 1) then
     sxi = lsecxStart+1
  else
     sxi = lsecxStart
  end if
  if (ylb == 1) then
     syi = lsecyStart+1
  else
     syi = lsecyStart
  end if

  do ix = sxi,lsecxEnd-1
     do iy = syi,lsecyEnd-1
        Txy(ix,iy,zlevel) = Txy(ix,iy,zlevel) + dtods*mu*( & 
             &                vx(ix,iy+1,zlevel)-vx(ix,iy,zlevel) + &
             &                vy(ix+1,iy,zlevel)-vy(ix,iy,zlevel))
     end do
  end do

!  do ix = sxi,lsecxEnd-1
!     do iy = lsecyStart,lsecyEnd
!           Txz(ix,iy,zlevel-1) = Txz(ix,iy,zlevel-1) +         &
!                &                dtods*(4/((2/mu)+(2/mu2)))*( & 
!                &                vx(ix,iy,zlevel)-vx(ix,iy,zlevel-1) + &
!                &                vz(ix+1,iy,zlevel-1)-vz(ix,iy,zlevel-1))
!     end do
!  end do
!  do ix = lsecxStart,lsecxEnd
!     do iy = syi,lsecyEnd-1
!           Tyz(ix,iy,zlevel-1) = Tyz(ix,iy,zlevel-1) +     &
!                &                dtods*(4/((2/mu)+(2/mu2)))*( & 
!                &                vy(ix,iy,zlevel)-vy(ix,iy,zlevel-1) + &
!                &                vz(ix,iy+1,zlevel-1)-vz(ix,iy,zlevel-1))
!     end do
!  end do
  do ix = sxi,lsecxEnd-1
     do iy = lsecyStart,lsecyEnd
        Txz(ix,iy,zlevel-1) = Txz(ix,iy,zlevel-1) +     &
             &                dtods*(4/((2/mu)+(2/mu2)))*( & 
             &                vx(ix,iy,zlevel)-vx(ix,iy,zlevel-1) + &
             &                vz(ix+1,iy,zlevel-1)-vz(ix,iy,zlevel-1))
     end do
  end do

  do ix = lsecxStart,lsecxEnd
     do iy = syi,lsecyEnd-1
        Tyz(ix,iy,zlevel-1) = Tyz(ix,iy,zlevel-1) +     &
             &                dtods*(4/((2/mu)+(2/mu2)))*( & 
             &                vy(ix,iy,zlevel)-vy(ix,iy,zlevel-1) + &
             &                vz(ix,iy+1,zlevel-1)-vz(ix,iy,zlevel-1))
    end do
  end do
end Subroutine interface_stress

  