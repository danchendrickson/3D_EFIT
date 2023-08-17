Subroutine map_transducer(startx,starty,endx,endy,numz,ly,lx,  &
    &   tposx1,tposx2,tposy1,tposy2,a,b,c,d)
    
    Integer, intent(in) :: startx,starty
    Integer, intent(in) :: endx,endy,numz,ly,lx
    Integer, intent(in) :: tposx1,tposx2,tposy1,tposy2
    Integer, intent(out) :: a,b,c,d
    
    a = 0
    b = 0
    c = 0
    d = 0
    
    ! if the left side is present
    if ((tposx1 >= startx) .AND. (tposx1 <= endx)) then
       ! if lower left corner is present
       if ((tposy1 >= starty) .AND. (tposy1 <= endy)) then
          a = tposx1 - (startx-1)
          c = tposy1 - (starty-1)
          ! if lower right corner is present
          if ((tposx2 >= startx) .AND. (tposx2 <= endx)) then
             b = tposx2 - (startx-1)
             ! else transducer goes to right edge
          else
             b = lx
          end if
          ! if upper left corner is present
          if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
             d  = tposy2 - (starty-1)
             ! else transducer goes to top edge
          else
             d = ly
          end if
          ! else if upper left corner is present
       else if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
          a = tposx1 - (startx-1)
          c = 1
          d = tposy2 - (starty-1)
          ! if upper right corner is present
          if ((tposx2 >= startx) .AND. (tposx2 <= endx)) then
             b = tposx2 - (startx-1)
             ! else transducer goes to right side
          else
             b = lx
          end if
          ! else if middle of left side (no corners)
       else if ((tposy1 < starty) .AND. (tposy2 > endy)) then
          a = tposx1 -(startx-1)
          c = 1
          d = ly
          ! if right side is present
          if ((tposx2 >= startx) .AND. (tposx2 <= endx)) then
             b = tposx2 - (startx-1)
             ! else transducer goes to right edge
          else
             b = lx
          end if
       end if
       ! else if right side is present
    else if ((tposx2 >= startx) .AND. (tposx2 <= endx)) then
       ! if lower right corner is present
       if ((tposy1 >= starty) .AND. (tposy1 <= endy)) then
          a = 1
          b = tposx2 - (startx-1)
          c = tposy1 - (starty-1)
          ! if upper right corner is present
          if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
             d = tposy2 - (starty-1)
             ! else transducer goes to top edge
          else
             d = ly
          end if
          ! else if upper right corner is present
       else if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
          a = 1
          b = tposx2 - (startx-1)
          c = 1
          d = tposy2 - (starty-1)
          ! else if middle of right side (no corners)
       else if ((tposy1 < starty) .AND. (tposy2 > endy)) then
          a = 1
          b = tposx2 - (startx-1)
          c = 1
          d = ly
       end if
       ! else if in between right and left sides
    else if ((tposx1 < startx) .AND. (tposx2 > endx)) then
       ! if bottom side is present
       if ((tposy1 >= starty) .AND. (tposy1 <= endy)) then
          a = 1
          b = lx
          c = tposy1 - (starty-1)
          ! if top side is present
          if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
             d = tposy2 - (starty-1)
             ! else transducer goes to top edge
          else
             d = ly
          end if
          ! else if top side is present
       else if ((tposy2 >= starty) .AND. (tposy2 <= endy)) then
          a = 1
          b = lx
          c = 1
          d = tposy2 - (starty-1)
          ! else if transducer covers entire region
       else if ((tposy1 < starty) .AND. (tposy2 > endy)) then
          a = 1
          b = lx
          c = 1
          d = ly
       end if
    end if
        
end Subroutine Map_transducer
                
                
                
