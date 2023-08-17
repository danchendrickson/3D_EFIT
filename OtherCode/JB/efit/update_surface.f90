Subroutine update_sflaw_vel(flx,fly,numz,dtodsp,lxstart,lystart,    &
     &   lzstart,lxend,lyend,lzend,frSz,flxStart,flxEnd,  &
     &   flyStart,flyEnd,flaw,xlb,xub,ylb,yub,zlb,zub, &
     &   fxlb,fxub,fylb,fyub,  &
     &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)

  Implicit none

 
  Integer, intent(in) :: flx,fly,numz,frSz
  Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
  Integer, intent(in) :: fxlb,fxub,fylb,fyub
  Integer, intent(in) :: a,b,c,starty
  Double Precision, intent(in) :: dtodsp
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fvx,fvy,fvz
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
  Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds

  Integer :: ix,iy,iz,tix,tiy,tiz

  !flx - lx, local number in x-dir
  !fly - ly, local number in y-dir
  !lnz - numz or tfxE,thickness
  !numz - numz, global number in z-dir
  !gnx - numx, global number in x-dir
  !gny - numy, global number in y-dir
  !lsx - startx, global start of x on local
  !lsy - starty or tfyS, global start of y on local
  !lex - endx, global end of x on local
  !ley - endy or tfyE, global end of y on local
  !sx_index - s1_index
  !sy_index - s2_index
  !dtodsp - dtodsp, make sure to use the right density
  !frSz - rSz, bottom of flaw region
  !flxStart - lrxStart, local region start x
  !flxEnd - lrxEnd, local region end x
  !flyStart - lryStart, local region start x
  !flyEnd - lryEnd, local region end y
  !flaw - flaw, what part of flaw is present
  !numS - numS, surface
  

!  write(*,*) 'F1 updating flaw vel'

  ! This is called if the current CPU has a flawed region on it

  ! Update any part of the cpu not in flawed region normalfly.
  ! Basicalfly reverse mapping of the region.
  ! flaw: x lower +10, x upper +20, y lower +1, y upper +2
  !      possibles - 1,2,3,4,10,11,12,13,20,21,22,23,30,31,32,33
  
  SELECT CASE (flaw)
     CASE (1) ! lower y bound of flaw
        !bottom before flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(lxstart,lystart,frSz,lxEnd,flyStart,lzend,  &
             &    flx,fly,numz,dtodsp,xlb,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,lxend,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)

     CASE (2) ! upper y bound of flaw
        !bottom under flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,flyEnd,     &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom after flaw
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,frSz,  &
             &    flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after of flaw
        call section_velocity(lxstart,flyEnd,frSz,lxend,lyend,lzend,    &
             &    flx,fly,numz,dtodsp,xlb,xub,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (3) ! lower and upper y bounds of flaw
        !bottom before flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(lxstart,lystart,frSz,lxEnd,flyStart,      &
             &    lzend,flx,fly,numz,dtodsp,xlb,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,lxend,flyEnd,    &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom after flaw
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,frSz,  &
             &    flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after of flaw
        call section_velocity(lxstart,flyEnd,frSz,lxend,lyend,lzend,    &
             &    flx,fly,numz,dtodsp,xlb,xub,1,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (4) ! flaw completely covers
        !bottom under flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)

     CASE (10) ! lower x bound of flaw
        !bottom left of flaw
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,lystart,frSz,flxStart,lyend,      &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right under flaw
        call section_velocity(flxStart,lystart,lzstart,lxend,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)

     CASE (11) ! lower x & y bounds of flaw
        !lower left
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart,  &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_velocity(flxStart,lystart,lzstart,lxend,flyStart,  &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(flxStart,lystart,frSz,lxend,flyStart,     &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,flyStart,frSz,lxend,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,xlb,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(flxStart,flyStart,lzstart,lxend,lyend,    &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        
     CASE (12) ! lower x & upper y bounds of flaw
        !bottom left of flaw
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyEnd,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,lystart,frSz,flxStart,flyEnd,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right under flaw
        call section_velocity(flxStart,lystart,lzstart,lxend,flyEnd,    &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
!        write(*,*)'F1 vel flag'
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
!        write(*,*)'F2 vel flag'
        !upper left
        call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,    &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(flxStart,flyEnd,frSz,lxend,lyend,lzend,   &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (13) ! lower x & y & upper y bounds of flaw
        !lower left
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart,  &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_velocity(flxStart,lystart,lzstart,lxend,flyStart,  &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(flxStart,lystart,frSz,lxend,flyStart,     &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,flyStart,frSz,flxStart,lyend,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(flxStart,flyStart,lzstart,lxend,flyEnd,   &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !upper left
        call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,    &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(flxStart,flyEnd,frSz,lxend,lyend,     &
             &    lzend,flx,fly,numz,dtodsp,0,xub,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (20) ! upper x bound of flaw
        ! left under flaw
        call section_velocity(lxstart,lystart,lzstart,flxEnd,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,       &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,lystart,frSz,lxend,lyend,lzend,    &
             &    flx,fly,numz,dtodsp, 0,xub,ylb,yub,0,zub,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (21) ! upper x & lower y bounds of flaw
        !bottom left before of flaw
        call section_velocity(lxstart,lystart,lzstart,flxEnd,flyStart,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left before of flaw
        call section_velocity(lxstart,lystart,frSz,flxEnd,flyStart,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,flxEnd,lyend,    &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,flyStart,frSz,lxend,lyend,lzend,   &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (22) ! upper x & y bounds of flaw
        !bottom under flaw
        call section_velocity(lxstart,lystart,lzstart,flxEnd,flyEnd,    &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
 !       write(*,*)'vel flag'
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
 !       write(*,*)'F4 vel flag'
        !bottom right of flaw
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,lystart,frSz,lxend,flyEnd,lzend,   &
             &    flx,fly,numz,dtodsp,0,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(lxstart,flyEnd,lzstart,flxEnd,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(lxstart,flyEnd,frSz,flxEnd,lyend,lzend,   &
             &    flx,fly,numz,dtodsp,xlb,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (23) ! lower y & upper x & y bounds of flaw
        !bottom left before of flaw
        call section_velocity(lxstart,lystart,lzstart,flxEnd,flyStart,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left before of flaw
        call section_velocity(lxstart,lystart,frSz,flxEnd,flyStart,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,flxEnd,flyEnd,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,flyStart,lzstart,lxend,flyEnd,     &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,flyStart,frSz,lxend,flyEnd,        &
             &    lzend,flx,fly,numz,dtodsp,0,xub,0,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(lxstart,flyEnd,lzstart,flxEnd,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(lxstart,flyEnd,frSz,flxEnd,lyend,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,        &
             &    lzend,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (30) ! upper & lower x bounds of flaw
        !bottom left of flaw
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,lystart,frSz,flxStart,lyend,  &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right under flaw
        call section_velocity(flxStart,lystart,lzstart,flxEnd,lyend,    &
             &    frSz,flx,fly,numz,dtodsp,0,0,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,lystart,frSz,lxend,lyend,lzend,    &
             &    flx,fly,numz,dtodsp,0,xub,ylb,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (31) ! lower x & y & upper x bounds of flaw
        !lower left
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart,  &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_velocity(flxStart,lystart,lzstart,flxEnd,flyStart, &
             &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(flxStart,lystart,frSz,flxEnd,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,flyStart,frSz,flxStart,lyend,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(flxStart,flyStart,lzstart,flxEnd,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,flyStart,frSz,lxend,lyend,lzend,   &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (32) ! upper x & y & lower x bounds of flaw
        !bottom left of flaw
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyEnd,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,lystart,frSz,flxStart,flyEnd,     &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_velocity(flxStart,lystart,lzstart,flxEnd,flyEnd,   &
             &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !lbottom right of flaw
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,lystart,frSz,lxend,flyEnd,lzend,   &
             &    flx,fly,numz,dtodsp,0,xub,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper left
        call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,    &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(flxStart,flyEnd,lzstart,flxEnd,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(flxStart,flyEnd,frSz,flxEnd,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,0,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (33) ! upper & lower x & y bound of flaw
        !lower left
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart, &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_velocity(flxStart,lystart,lzstart,flxEnd,flyStart, &
             &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(flxStart,lystart,frSz,flxEnd,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,0,ylb,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,    &
             &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_velocity(lxstart,flyStart,lzstart,flxStart,flyEnd, &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,flyStart,frSz,flxStart,flyEnd,    &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(flxStart,flyStart,lzstart,flxEnd,flyEnd,  &
             &    frSz,flx,fly,numz,dtodsp,0,0,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
             &   flyStart,flyEnd,frSz,lzend,fvx,fvy,fvz,  &
             &   fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
        !bottom right of flaw
        call section_velocity(flxEnd,flyStart,lzstart,lxend,flyEnd,     &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_velocity(flxEnd,flyStart,frSz,lxend,flyEnd,lzend,  &
             &    flx,fly,numz,dtodsp,0,xub,0,0,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper left
        call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,    &
             &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(flxStart,flyEnd,lzstart,flxEnd,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_velocity(flxStart,flyEnd,frSz,flxEnd,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,0,0,0,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,lzend,  &
             &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     end SELECT
   end Subroutine update_sflaw_vel

!=========================================================================

   Subroutine update_sflaw_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
        &   frSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
        &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
        &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
        &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)

     Implicit none
     
     Integer, intent(in) :: flx,fly,numz,frSz
     Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
     Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
     Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
     Integer, intent(in) :: a,b,c
     Double Precision, intent(in) :: dtods,l2m,lambda,mu
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: fvx,fvy,fvz
     Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     
     Integer :: ix,iy,iz,tix,tiy,tiz
     
     !flx - lx, local number in x-dir
     !fly - ly or ltfyE, local number in y-dir
     !lnz - numz or tfzE, thickness
     !numz - numz, number in z-dir
     !sx_index - s1_index
     !sy_index - s2_index
     !dtods - dtods
     !l2m - l2m
     !lambda - lambda
     !mu - mu
     !frSz - rSz, bottom of flaw
     !flxStart - lrxStart, local region start x
     !flxEnd - lrxEnd, local region end x
     !flyStart - lryStart, local region start y
     !flyEnd - lryEnd, local region end y
     !flaw - flaw, which part of flaw is present
     !numS - numS, surface or not
     
     SELECT CASE (flaw)
     CASE (1)
        !bottom before flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
     
     CASE (2)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (3)
        !bottom before flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

      CASE (4)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)

     CASE (10)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)


     CASE (11)
        !lower left
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyStart,lzend,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(flxStart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)

     CASE (12)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyEnd,frSz,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !upper left
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(flxStart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (13)
        !lower left
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyStart,lzend,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(flxStart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !upper left
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(flxStart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (20)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (21)
        !bottom before flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (22)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_stress(flxEnd,flyEnd,lzstart,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (23)
        !bottom before flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_stress(flxEnd,flyEnd,lzstart,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (30)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (31)
        !lower left
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyStart,lzend,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(flxStart,lystart,frSz,  &
             &    flxEnd,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (32)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyEnd,frSz,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper left
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(flxStart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_stress(flxEnd,flyEnd,lzstart,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (33)
        !lower left
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyStart,lzend,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom before flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_stress(flxStart,lystart,frSz,  &
             &    flxEnd,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !surface
        call  update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
             &   flxStart,flxEnd,flyStart,flyEnd,frSz,lzend,  &
             &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper left
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(flxStart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_stress(flxEnd,flyEnd,lzstart,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     end SELECT

   end Subroutine update_sflaw_stress


   Subroutine update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
        &   flyStart,flyEnd,rSz,rEz,  &
        &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c,starty)

     Integer, intent(in) :: flx,fly,numz,rSz,rEz
     Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd
     Integer, intent(in) :: a,b,c,starty
     Double Precision, intent(in) :: dtodsp
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: vx,vy,vz
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds

     Integer :: ix,iy,iz,tix,tiy,tiz,syi

     if (starty == 1) then
        syi = flyStart
     else
        syi = flyStart+1
     end if

     !! Update velocities
     ! Update vx
     do ix = flxStart,flxEnd-1
        do tiy = syi,flyEnd
           do tiz = rSz+1,rEz
              SELECT CASE (bounds(ix-flxStart+1,tiy-flyStart+1,tiz-rSz+1))
              CASE (-4)
                 ! zero density
                 vx(ix,tiy,tiz) = 0.0
                 !                    write(*,*)' Flag 1'
              CASE (0,2,20,22)
                 ! interior of material
                 if (tiy /= flyStart) then
                    ! accounts for edge of CPU
!                    write(*,*)' Flag 2'        do ix = s1_index,lx-1
                    vx(ix,tiy,tiz) = vx(ix,tiy,tiz) + dtodsp* &
                         &     (Txx(ix+1,tiy,tiz)-Txx(ix,tiy,tiz) + &
                         &      Txy(ix,tiy,tiz)-Txy(ix,tiy-1,tiz) + &
                         &      Txz(ix,tiy,tiz)-Txz(ix,tiy,tiz-1))
                 end if
              CASE (100,101,102,110,111,112,120,121,122)
                 ! left boundary
!                 write(*,*)' Flag 3'
                 vx(ix,tiy,tiz) = vx(ix,tiy,tiz) +      &
                        &           dtodsp*(2*Txx(ix+1,tiy,tiz))
              CASE (200,201,202,210,211,212,220,221,222)
                 ! right boundary
!                 write(*,*)' Flag 4'
                 vx(ix,tiy,tiz) = vx(ix,tiy,tiz) +      &
                        &           dtodsp*(-2*Txx(ix,tiy,tiz))
              END SELECT
           end do
        end do
     end do
 !    write(*,*)'vx flag'
    ! Update vy
     do tix = flxStart+1,flxEnd
        do iy = flyStart,flyEnd-1
           do tiz = rSz+1,rEz
              SELECT CASE (bounds(tix-flxStart+1,iy-flyStart+1,tiz-rSz+1))
              CASE (-4)
                 ! zero density
                 vy(tix,iy,tiz) = 0.0
              CASE (0,2,200,202)
                 ! interior of material
!                 if ((tix /= 1) .AND. (tiz /= 1) .AND. (iy /= fly)) then
                    ! accounts for edge of CPU
                    vy(tix,iy,tiz) = vy(tix,iy,tiz) + dtodsp* &
                         &    (Txy(tix,iy,tiz)-Txy(tix-1,iy,tiz) + &
                         &     Tyy(tix,iy+1,tiz)-Tyy(tix,iy,tiz) + &
                         &     Tyz(tix,iy,tiz)-Tyz(tix,iy,tiz-1))
!                 end if
              CASE (10,11,12,110,111,112,210,211,212)
                 ! front boundary
                 vy(tix,iy,tiz) = vy(tix,iy,tiz) +      &
                        &           dtodsp*(2*Tyy(tix,iy+1,tiz))
              CASE (20,21,22,120,121,122,220,221,222)
                 ! back boundary
                 vy(tix,iy,tiz) = vy(tix,iy,tiz) +      &
                        &           dtodsp*(-2*Tyy(tix,iy,tiz)) 
              END SELECT
           end do
        end do
     end do
!     write(*,*)'vy flag'
     ! Update vz
     do tix = flxStart+1,flxEnd
        do tiy = syi,flyEnd
           do iz = rSz,rEz-1
              SELECT CASE (bounds(tix-flxStart+1,tiy-flyStart+1,iz-rSz+1))
              CASE (-4)
                 ! zero density
                 vz(tix,tiy,iz) = 0.0
              CASE (0,20,200,220)
                 if (tiy /= flyStart) then
                    ! accounts for edge of CPU
                    vz(tix,tiy,iz) = vz(tix,tiy,iz) + dtodsp* &
                         &          (Txz(tix,tiy,iz)-Txz(tix-1,tiy,iz) + &
                         &          Tyz(tix,tiy,iz)-Tyz(tix,tiy-1,iz) + &
                         &          Tzz(tix,tiy,iz+1)-Tzz(tix,tiy,iz))
                 end if
              CASE (1,11,21,101,111,121,201,211,221)
                 ! bottom boundary
                 vz(tix,tiy,iz) = vz(tix,tiy,iz) +  &
                        &           dtodsp*(2*Tzz(tix,tiy,iz+1))
              CASE (2,12,22,102,112,122,202,212,222)
                 ! top boundary
                 vz(tix,tiy,iz) = vz(tix,tiy,iz) +      &
                        &           dtodsp*(-2*Tzz(tix,tiy,iz))
              END SELECT
           end do
        end do
     end do
!     write(*,*)'vz flag'
     
   end Subroutine update_surf_vel


   Subroutine update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
        &   flxStart,flxEnd,flyStart,flyEnd,rSz,rEz,  &
        &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c)

     Integer, intent(in) :: flx,fly,numz,rSz,rEz
     Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd
     Integer, intent(in) :: a,b,c
     Double Precision, intent(in) :: dtods,l2m,lambda,mu
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: vx,vy,vz
     Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)     &
                &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds

     Integer :: ix,iy,iz,tix,tiy,tiz,bd

!     write(*,*)'flxStart:',flxStart,'flyStart:',flyStart,'rSz:',rSz,rEz
        do tix = flxStart+1,flxEnd
           do tiy = flyStart+1,flyEnd
              do tiz = rSz+1,rEz
                 bd = bounds(tix-flxStart+1,tiy-flyStart+1,tiz-rSz+1)
                 SELECT CASE (bd)
                 CASE (-4,1,10,11,12,21,22,100,101,102,110,111,112,     &
                  &   120,121,122,201,210,211,212,221,222,2,20,200,202,220)
                    ! zero density or lower boundary
                    Txx(tix,tiy,tiz) = 0.0
                    Tyy(tix,tiy,tiz) = 0.0
                    Tzz(tix,tiy,tiz) = 0.0

                 CASE (0)
                    ! interior of material
                    Txx(tix,tiy,tiz) = Txx(tix,tiy,tiz) + dtods*( l2m* &
                         &      (vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz)) + &
                         &      lambda* &
                         &      (vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz) + &
                         &      vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)))
                    Tyy(tix,tiy,tiz) = Tyy(tix,tiy,tiz) + dtods*( l2m* &
                         &      (vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz)) + &
                         &      lambda* &
                         &      (vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz) + &
                         &      vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)))
                    Tzz(tix,tiy,tiz) = Tzz(tix,tiy,tiz) + dtods*( l2m* &
                         &      (vz(tix,tiy,tiz)-vz(tix,tiy,tiz-1)) + &
                         &      lambda* &
                         &      (vx(tix,tiy,tiz)-vx(tix-1,tiy,tiz) + &
                         &      vy(tix,tiy,tiz)-vy(tix,tiy-1,tiz)))
                 end SELECT
              end do
           end do
        end do

        do ix = flxStart,flxEnd-1
           do iy = flyStart,flyEnd-1
              do iz = rSz+1,rEz
                 bd = bounds(ix-flxStart+1,iy-flyStart+1,iz-rSz+1)
                 SELECT CASE (bd)
                 CASE (-4,10,11,12,20,22,100,101,102,110,111,112,120,   &
                    &   121,122,200,201,202,210,211,212,220,221,222)
                    ! zero density or bc
                    Txy(ix,iy,iz) = 0.0
                    
                 CASE (0,1,2)
                    Txy(ix,iy,iz) = Txy(ix,iy,iz) + dtods*mu*( & 
                         &                vx(ix,iy+1,iz)-vx(ix,iy,iz) + &
                         &                vy(ix+1,iy,iz)-vy(ix,iy,iz))
                 end SELECT
              end do
           end do
        end do

        do ix = flxStart,flxEnd-1
           do iy = flyStart+1,flyEnd-1
              do iz = rSz,rEz
                 bd = bounds(ix-flxStart+1,iy-flyStart+1,iz-rSz+1)
                 SELECT CASE (bd)
                 CASE (-4,1,2,11,12,21,22,100,101,102,110,111,112,120,  &
                    &   121,122,200,201,202,210,211,212,220,221,222)
                    ! zero density
                    Txz(ix,iy,iz) = 0.0
                    
                 CASE (0,10,20)
                    Txz(ix,iy,iz) = Txz(ix,iy,iz) + dtods*mu*( & 
                         &                vx(ix,iy,iz+1)-vx(ix,iy,iz) + &
                         &                vz(ix+1,iy,iz)-vz(ix,iy,iz))
                 end SELECT
              end do
           end do
        end do

        do ix = flxStart+1,flxEnd-1
           do iy = flyStart,flyEnd-1
              do iz = rSz,rEz
                 bd = bounds(ix-flxStart+1,iy-flyStart+1,iz-rSz+1)
                 SELECT CASE (bd)
                 CASE (-4,1,2,10,11,12,20,21,22,101,102,110,111,112,    &
                        &   120,121,122,201,202,210,211,212,220,221,222)
                    ! zero density
                    Tyz(ix,iy,iz) = 0.0
                    
                 CASE (0,100,200)
                    Tyz(ix,iy,iz) = Tyz(ix,iy,iz) + dtods*mu*( & 
                         &                vy(ix,iy,iz+1)-vy(ix,iy,iz) + &
                         &                vz(ix,iy+1,iz)-vz(ix,iy,iz))
                 end SELECT
              end do
           end do
        end do

      end Subroutine update_surf_stress
