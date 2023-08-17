Subroutine section_velocity(lsecxStart,lsecyStart,lseczStart,  &
     &    lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz,dtodsp,  &
     &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
                    
  Implicit none

  Integer, intent(in) :: lsecxStart,lsecyStart,lseczStart
  Integer, intent(in) :: lsecxEnd,lsecyEnd,lseczEnd
  Integer, intent(in) :: lx,ly,numz,xlb,xub,ylb,yub,zlb,zub
  Double Precision, intent(in) :: dtodsp
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout) :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)      &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz

  Integer :: ix,iy,iz,tix,tiy,tiz

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
        do tiz = lseczStart+1,lseczEnd
           vx(ix,tiy,tiz) = vx(ix,tiy,tiz) + dtodsp*( &
                &	Txx(ix+1,tiy,tiz)-Txx(ix,tiy,tiz) + &
                &  Txy(ix,tiy,tiz)-Txy(ix,tiy-1,tiz) + &
                &  Txz(ix,tiy,tiz)-Txz(ix,tiy,tiz-1))
        end do
     end do
  end do

  if (zlb == 1) then
     szi = lseczStart
  else
     szi = lseczStart+1
  end if

  !  update boundary points
  if (xlb == 1) then
     do tiy = lsecyStart,lsecyEnd
        do tiz = szi,lseczEnd
           vx(lsecxStart,tiy,tiz) = vx(lsecxStart,tiy,tiz) + &
                &   dtodsp*(2*Txx(lsecxStart+1,tiy,tiz))
        end do
     end do
  end if
  if (xub == 1) then
     do tiy = lsecyStart,lsecyEnd
        do tiz = szi,lseczEnd
           vx(lsecxEnd,tiy,tiz) = vx(lsecxEnd,tiy,tiz) + &
                &                   dtodsp*(-2*Txx(lsecxEnd,tiy,tiz))
        end do
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
        do tiz = lseczStart+1,lseczEnd
           vy(tix,iy,tiz) = vy(tix,iy,tiz) + dtodsp*( &
                &	Txy(tix,iy,tiz)-Txy(tix-1,iy,tiz) + &
                &	Tyy(tix,iy+1,tiz)-Tyy(tix,iy,tiz) + &
                &	Tyz(tix,iy,tiz)-Tyz(tix,iy,tiz-1))
        end do
     end do
  end do
  !  update boundary points
  if (ylb == 1) then
     do tix = lsecxStart,lsecxEnd
        do tiz = szi,lseczEnd
           vy(tix,lsecyStart,tiz) = vy(tix,lsecyStart,tiz) + &
                &   dtodsp*(2*Tyy(tix,lsecyStart+1,tiz))
        end do
     end do
  end if
  if (yub == 1) then
     do tix = lsecxStart,lsecxEnd
        do tiz = szi,lseczEnd
           vy(tix,lsecyEnd,tiz) = vy(tix,lsecyEnd,tiz) + &
                &  dtodsp*(-2*Tyy(tix,lsecyEnd,tiz)) 
        end do
     end do
  end if

  if (zlb == 1) then
     szi = lseczStart + 1
  else
     szi = lseczStart
  end if
  !   go down in z direction
  !   update interior points
  do tix = lsecxStart+1,lsecxEnd
     do tiy = lsecyStart+1,lsecyEnd
        do iz = szi,lseczEnd-1
           vz(tix,tiy,iz) = vz(tix,tiy,iz) + dtodsp*( &
                &	Txz(tix,tiy,iz)-Txz(tix-1,tiy,iz) + &
                &	Tyz(tix,tiy,iz)-Tyz(tix,tiy-1,iz) + &
                &	Tzz(tix,tiy,iz+1)-Tzz(tix,tiy,iz))
        end do
     end do
  end do
  !  update boundary points
  if (zlb == 1) then
     do tix = lsecxStart,lsecxEnd
        do tiy = lsecyStart,lsecyEnd
           vz(tix,tiy,lseczStart) = vz(tix,tiy,lseczStart) + &
                &     dtodsp*(2*Tzz(tix,tiy,lseczStart+1))
        end do
     end do
  end if
  if (zub == 1) then
     do tix = lsecxStart,lsecxEnd
        do tiy = lsecyStart,lsecyEnd
           vz(tix,tiy,lseczEnd) = vz(tix,tiy,lseczEnd) + &
                &                   dtodsp*(-2*Tzz(tix,tiy,lseczEnd)) 
        end do
     end do
  end if
end Subroutine section_velocity



Subroutine section_stress(lsecxStart,lsecyStart,lseczStart,  &
     &    lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz,  &
     &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,zub, &
     &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
                   
  Implicit none

  Integer, intent(in) :: lsecxStart,lsecyStart,lseczStart
  Integer, intent(in) :: lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
  Double Precision, intent(in) :: dtods,l2m,lambda,mu
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(in) :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)      &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz

  Integer :: ix,iy,iz,tix,tiy,tiz

  Integer :: sxi,syi,szi

  do tix = lsecxStart+1,lsecxEnd
     do tiy = lsecyStart+1,lsecyEnd
        do tiz = lseczStart+1,lseczEnd
           Txx(tix,tiy,tiz) = Txx(tix,tiy,tiz) + dtods*( &
                &	l2m*(vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz)) + &
                &	lambda*(vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz) + &
                &	vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)))
           Tyy(tix,tiy,tiz) = Tyy(tix,tiy,tiz) + dtods*( &
                &	l2m*(vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz)) + &
                &	lambda*(vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz) + &
                &	vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)))
           Tzz(tix,tiy,tiz) = Tzz(tix,tiy,tiz) + dtods*( &
                &	l2m*(vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)) + &
                &	lambda*(vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz) + &
                &	vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz)))
        end do
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
  if (zlb == 1) then
     szi = lseczStart+1
  else
     szi = lseczStart
  end if
  

  do ix = sxi,lsecxEnd-1
     do iy = syi,lsecyEnd-1
        do iz = lseczStart,lseczEnd
           Txy(ix,iy,iz) = Txy(ix,iy,iz) + dtods*mu*( & 
                &                vx(ix,iy+1,iz)-vx(ix,iy,iz) + &
                &                vy(ix+1,iy,iz)-vy(ix,iy,iz))
        end do
     end do
  end do
  do ix = sxi,lsecxEnd-1
     do iy = lsecyStart,lsecyEnd
        do iz = szi,lseczEnd-1
           Txz(ix,iy,iz) = Txz(ix,iy,iz) + dtods*mu*( & 
                &                vx(ix,iy,iz+1)-vx(ix,iy,iz) + &
                &                vz(ix+1,iy,iz)-vz(ix,iy,iz))
        end do
     end do
  end do
  do ix = lsecxStart,lsecxEnd
     do iy = syi,lsecyEnd-1
        do iz = szi,lseczEnd-1
           Tyz(ix,iy,iz) = Tyz(ix,iy,iz) + dtods*mu*( & 
                &                vy(ix,iy,iz+1)-vy(ix,iy,iz) + &
                &                vz(ix,iy+1,iz)-vz(ix,iy,iz))
        end do
     end do
  end do

end Subroutine section_stress

