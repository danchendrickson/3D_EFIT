Subroutine find_bounds(lxS,lxE,lyS,lyE,rSz,rEz,starty,dens,bounds,Pid)

  Integer, Intent(in) :: lxS,lxE,lyS,lyE,rSz,rEz,starty,Pid
  Double Precision, Dimension(lxS:lxE+2,lyS:lyE+2,rSz:rEz),     & 
        &   Intent(inout) :: dens(:,:,:)
  Integer, Dimension(lxS:lxE,lyS:lyE,rSz:rEz), Intent(out) :: bounds(:,:,:)

  Integer :: x,y,z,Bx,By,Bz,repeat_flag

  Bz = 0
  Bx = 0
  By = 0
  repeat_flag = 1
  boounds = -4

!  write(*,*)'size:',Pid,lxE,lyE,rEz

  do while (repeat_flag == 1)
     repeat_flag = 0
     do x = 1,lxE
        do y = 1,lyE
           do z = rSz,rEz
              if (dens(x+1,y+1,z) <= 0.0) then
                 bounds(x,y,z) = -4
              elseif (dens(x+1,y+1,z) > 0.0) then
                 if (z == rEz) then
                    Bz = Bz +2
                 else if (dens(x+1,y+1,z+1) <= 0.0) then
                    Bz = Bz + 2
                 end if
!                 if ((y == 1) .and. (starty == 1)) then
!                    By = 10
!                 end if
                 if (dens(x+1,y,z) <= 0.0) then
                    By = By + 10
                 end if
                 if (dens(x+1,y+2,z) <= 0.0) then
                    By = By + 20
                 end if
                 if (dens(x,y+1,z) <= 0.0) then
                    Bx = Bx + 100
                 end if
                 if (dens(x+2,y+1,z) <= 0.0) then
                    Bx = Bx + 200
                 end if
                 bounds(x,y,z) = Bx + By + Bz
                 if ((Bx == 300) .or. (By == 30)) then
                    dens(x+1,y+1,z) = -222.0
                    repeat_flag = 1
                 end if
                 Bx = 0
                 By = 0
                 Bz = 0
              end if
           end do
        end do
     end do
!     bounds(:,:,rEz) = -4
     write(*,*)'Pid:',Pid,'Repeat_flag',repeat_flag
  end do

end Subroutine find_bounds
