Subroutine map_region(startx,starty,endx,endy,lx,ly,   &
    &   rxStart,rxEnd,ryStart,ryEnd,lxStart,lxEnd,lyStart,lyEnd,flaw)

    Integer, intent(in) :: startx,starty
    Integer, intent(in) :: endx,endy,ly,lx
    Integer, intent(in) :: rxStart,rxEnd,ryStart,ryEnd
    Integer, intent(out) :: lxStart,lxEnd,lyStart,lyEnd,flaw

    lxStart = 0
    lxEnd = 0
    lyStart = 0
    lyEnd = 0
    flaw = 0

   ! if the left side is present
    if ((rxStart >= startx) .AND. (rxStart <= endx)) then
        ! if front left corner is present
        if ((ryStart >= starty) .AND. (ryStart <= endy)) then
            lxStart = rxStart - (startx-1)
            flaw = flaw + 10
            lyStart = ryStart - (starty-1)
            flaw = flaw + 1
            ! if front right corner is present
            if ((rxEnd >= startx) .AND. (rxEnd <= endx)) then
                lxEnd = rxEnd - (startx-1)
                flaw = flaw + 20
            ! else region goes to right edge
            else
                lxEnd = lx
            end if
            ! if back left corner is present
            if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
                lyEnd  = ryEnd - (starty-1)
                flaw = flaw + 2 
            ! else region goes to back edge
            else
                lyEnd = ly
            end if
        ! else if back left corner is present
        else if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
            lxStart = rxStart - (startx-1)
            flaw = flaw + 10
            lyStart = 1
            lyEnd = ryEnd - (starty-1)
            flaw = flaw + 2
            ! if back right corner is present
            if ((rxEnd >= startx) .AND. (rxEnd <= endx)) then
                lxEnd = rxEnd - (startx-1)
                flaw = flaw + 20
            ! else region goes to right side
            else
                lxEnd = lx
            end if
        ! else if middle of left side (no corners)
        else if ((ryStart < starty) .AND. (ryEnd > endy)) then
            lxStart = rxStart -(startx-1)
            flaw = flaw + 10
            lyStart = 1
            lyEnd = ly
            ! if right side is present
            if ((rxEnd >= startx) .AND. (rxEnd <= endx)) then
                lxEnd = rxEnd - (startx-1)
                flaw = flaw + 20
            ! else region goes to right edge
            else
                lxEnd = lx
            end if
        end if
    ! else if right side is present
    else if ((rxEnd >= startx) .AND. (rxEnd <= endx)) then
        ! if front right corner is present
        if ((ryStart >= starty) .AND. (ryStart <= endy)) then
            lxStart = 1
            lxEnd = rxEnd - (startx-1)
            flaw = flaw + 20
            lyStart = ryStart - (starty-1)
            flaw = flaw + 1
            ! if back right corner is present
            if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
                lyEnd = ryEnd - (starty-1)
                flaw = flaw + 2
            ! else region goes to back edge
            else
                lyEnd = ly
            end if
        ! else if back right corner is present
        else if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
            lxStart = 1
            lxEnd = rxEnd - (startx-1)
            flaw = flaw + 20
            lyStart = 1
            lyEnd = ryEnd - (starty-1)
            flaw = flaw + 2
        ! else if middle of right side (no corners)
        else if ((ryStart < starty) .AND. (ryEnd > endy)) then
            lxStart = 1
            lxEnd = rxEnd - (startx-1)
            flaw = flaw + 20
            lyStart = 1
            lyEnd = ly
        end if
    ! else if in between right and left sides
    else if ((rxStart < startx) .AND. (rxEnd > endx)) then
        ! if front side is present
        if ((ryStart >= starty) .AND. (ryStart <= endy)) then
            lxStart = 1
            lxEnd = lx
            lyStart = ryStart - (starty-1)
            flaw = flaw + 1
            ! if back side is present
            if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
                lyEnd = ryEnd - (starty-1)
                flaw = flaw + 2
            ! else region goes to back edge
            else
                lyEnd = ly
            end if
        ! else if back side is present
        else if ((ryEnd >= starty) .AND. (ryEnd <= endy)) then
            lxStart = 1
            lxEnd = lx
            lyStart = 1
            lyEnd = ryEnd - (starty-1)
            flaw = flaw + 2
        ! else if region covers entire area
        else if ((ryStart < starty) .AND. (ryEnd > endy)) then
            lxStart = 1
            lxEnd = lx
            lyStart = 1
            lyEnd = ly
            flaw = 4
        end if
     end if

   end Subroutine map_region
   
   Subroutine change_flaw(startx,starty,endx,endy,rxStart,ryStart,  &
        &   rxEnd,ryEnd,fxlb,fxub,fylb,fyub,flaw,Pid)

    Integer, intent(in) :: startx,starty,endx,endy
    Integer, intent(in) :: rxStart,rxEnd,ryStart,ryEnd
    Integer, intent(out) :: fxlb,fxub,fylb,fyub
    Integer, intent(inout) :: flaw,Pid

    fxlb = 0
    fxub = 0
    fylb = 0
    fyub = 0


! change the value of flaw according to if edge is at edge of cpu
     SELECT CASE (flaw)
        CASE (0,4)
!              write(*,*)'f1 Pid:',Pid
        CASE (1)
           fylb = 1
           if (ryStart == starty) then
              flaw = 4
           end if
        CASE (2)
           fyub = 1
           if (ryEnd == endy) then
              flaw = 4
           end if
        CASE (3)
           fylb = 1
           fyub = 1
           if (ryStart == starty) then
              if (ryEnd == endy) then
                 flaw = 4
              else
                 flaw = 2
              end if
           else if (ryEnd == endy) then
              flaw = 1
           end if
        CASE (10)
           fxlb = 1
           if (rxStart == startx) then
              flaw = 4
!              write(*,*)'f1 Pid:',Pid
           end if
        CASE (11)
           fxlb = 1
           fylb = 1
           if (rxStart == startx) then
              if (ryStart == starty) then
                 flaw = 4
!                 write(*,*)'f1 Pid:',Pid
              else
                 flaw = 1
              end if
           else if (ryStart == starty) then
              flaw = 10
           end if
        CASE (12)
           fxlb = 1
           fyub = 1
           if (rxStart == startx) then
              if (ryEnd == endy) then
                 flaw = 4
!                 write(*,*)'f1 Pid:',Pid
              else
                 flaw = 2
              end if
           else if (ryEnd == endy) then
              flaw = 10
           end if
        CASE (13)
           fxlb = 1
           fylb = 1
           fyub = 1
           if (rxStart == startx) then
              if (ryStart == starty) then
                 if (ryEnd == endy) then
                    flaw = 4
                 else
                    flaw = 2
                 end if
              else if (ryEnd == endy) then
                 flaw = 1
              end if
           else if (ryStart == starty) then
              if (ryEnd == endy) then
                 flaw = 10
              else
                 flaw = 12
              end if
           else if (ryEnd == endy) then
              flaw = 11
           end if
        CASE (20)
           fxub = 1
           if (rxEnd == endx) then
              flaw = 4
!              write(*,*)'f1 Pid:',Pid
           end if
        CASE (21)
           fxub = 1
           fylb = 1
           if (rxEnd == endx) then
              if (ryStart == starty) then
                 flaw = 4
!                 write(*,*)'f1 Pid:',Pid
              else
                 flaw = 1
              end if
           else if (ryStart == starty) then
              flaw = 20
           end if
        CASE (22)
           fxub = 1
           fyub = 1
           if (rxEnd == endx) then
              if (ryEnd == endy) then
                 flaw = 4
!                 write(*,*)'f1 Pid:',Pid
              else
                 flaw = 2
              end if
           else if (ryEnd == endy) then
              flaw = 20
           end if
        CASE (23)
           fxub = 1
           fylb = 1
           fyub = 1
           if (rxEnd == endx) then
              if (ryStart == starty) then
                 if (ryEnd == endy) then
                    flaw = 4
                 else 
                    flaw = 2
                 end if
              else if (ryEnd == endy) then
                 flaw = 1
              end if
           else if (ryStart == starty) then
              if (ryEnd == endy) then
                 flaw = 20
              else
                 flaw = 22
              end if
           else if (ryEnd == endy) then
              flaw = 21
           end if
        CASE (30)
           fxlb = 1
           fxub = 1
           if (rxStart == startx) then
              if (rxEnd == endx) then
                 flaw = 4
              else
                 flaw = 10
              end if
           else if (rxEnd == endx) then
              flaw = 20
           end if
        CASE (31)
           fxlb = 1
           fxub = 1
           fylb = 1
           if (rxStart == startx) then
              if (rxEnd == endx) then
                 if (ryStart == starty) then
                    flaw = 4
                 else
                    flaw = 1
                 end if
              else if (ryStart == starty) then
                 flaw = 20
              else
                 flaw = 21
              end if
           else if (rxEnd == endx) then
              if (ryStart == starty) then
                 flaw = 10
              else
                 flaw = 11
              end if
           else if (ryStart == starty) then
              flaw = 30
           end if
        CASE (32)
           fxlb = 1
           fxub = 1
           fyub = 1
           if (fxStart == startx) then
              if (rxEnd == endx) then
                 if (ryEnd == endy) then
                    flaw = 4
                 else
                    flaw = 2
                 end if
              else if (ryEnd == endy) then
                 flaw = 20
              else
                 flaw = 22
              end if
           else if (rxEnd == endx) then
              if (ryEnd == endy) then
                 flaw = 10
              else
                 flaw = 12
              end if
           else if (ryEnd == endy) then
              flaw = 30
           end if
        CASE (33)
           fxlb = 1
           fxub = 1
           fylb = 1
           fyub = 1
           if (rxStart == startx) then
              if (rxEnd == endx) then
                 if (ryStart == starty) then
                    if (ryEnd == endy) then
                       flaw = 4
                    else 
                       flaw = 2
                    end if
                 else if (ryEnd == endy) then
                    flaw = 1
                 else
                    flaw = 3
                 end if
              else if (ryStart == starty) then
                 if (ryEnd == endy) then
                    flaw = 20
                 else
                    flaw = 22
                 end if
              else if (ryEnd == endy) then
                 flaw = 21
              else
                 flaw = 23
              end if
           else if (rxEnd == endx) then
              if (ryStart == starty) then
                 if (ryEnd == endy) then
                    flaw = 10
                 else
                    flaw = 12
                 end if
              else if (ryEnd == endy) then
                 flaw = 11
              else
                 flaw = 13
              end if
           else if (ryStart == starty) then
              if (ryEnd == endy) then
                 flaw = 30
              else
                 flaw = 32
              end if
           else if (ryEnd == endy) then
              flaw = 31
           end if
        end SELECT

      end Subroutine change_flaw

   Subroutine map_T_region(starty,endy,ly,tryS,tryE,ltryS,ltryE) !,Sr,Er)

      Integer, intent(in) :: tryS
      Integer, intent(in) :: tryE
      Integer, intent(in) :: starty
      Integer, intent(in) :: endy
      Integer, intent(in) :: ly
      Integer, intent(out) :: ltryS
      Integer, intent(out) :: ltryE
!      Integer, intent(out) :: Sr,Er

      ltryS = 0
      ltryE = 0
!      Sr = 0
!      Er = 0

      ! if front is present
      if ((tryS >= starty) .AND. (tryS <= endy)) then
         ltryS = tryS - (starty - 1)
         Sr = 1
         ! if back is present
         if ((tryE >= starty) .AND. (tryE <= endy)) then
            ltryE = tryE - (starty - 1)
!            Er = 1
         ! region goes to back edge
         else
            ltryE = ly
         end if
      ! else if back is present
      else if ((tryE >= starty) .AND. (tryE <= endy)) then
         ltryE = tryE - (starty - 1)
!         Er = 1
         ltryS = 1
      ! else if in the middle
      else if ((tryS < starty) .AND. (tryE > endy)) then
         ltryS = 1
         ltryE = ly
      end if

    end Subroutine map_T_region
