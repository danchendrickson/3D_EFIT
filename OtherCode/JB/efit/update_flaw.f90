Subroutine update_flaw_vel(flx,fly,numz,dtodsp,lxstart,lystart,lzstart, &
     &   lxend,lyend,lzend,frSz,flxStart,flxEnd,flyStart,flyEnd,  &
     &   flaw,xlb,xub,ylb,yub,zlb,zub, &
     &   fxlb,fxub,fylb,fyub,  &
     &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

  Implicit none

 
  Integer, intent(in) :: flx,fly,numz,frSz
  Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
  Integer, intent(in) :: fxlb,fxub,fylb,fyub
  Double Precision, intent(in) :: dtodsp
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fvx,fvy,fvz
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz

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
        call section_velocity(lxstart,lystart,frSz,lxEnd,flyStart,      &
             &    lzend,flx,fly,numz,dtodsp,xlb,xub,ylb,1,0,zub  &
             &    ,fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,lxend,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (2) ! upper y bound of flaw
        !bottom under flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,flyEnd,     &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,frSz,  &
             &    flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after of flaw
        call section_velocity(lxstart,flyEnd,frSz,lxend,lyend,lzend,    &
             &    flx,fly,numz,dtodsp,xlb,xub,1,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (3) ! lower and upper y bounds of flaw
        !bottom before flaw
        call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top before flaw
        call section_velocity(lxstart,lystart,frSz,lxEnd,flyStart,      &
             &    lzend,flx,fly,numz,dtodsp,xlb,xub,ylb,1,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom under flaw
        call section_velocity(lxstart,flyStart,lzstart,lxend,flyEnd,    &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,0,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
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
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (10) ! lower x bound of flaw
        !bottom left of flaw
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_velocity(lxstart,lystart,frSz,flxStart,lyend,      &
             &    lzend,flx,fly,numz,dtodsp,xlb,1,ylb,yub,0,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right under flaw
        call section_velocity(flxStart,lystart,lzstart,lxend,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    CASE (11) ! lower x & y bounds of flaw
       !lower left
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart, &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom before flaw
       call section_velocity(flxStart,lystart,lzstart,lxend,flyStart,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top before flaw
       call section_velocity(flxStart,lystart,frSz,lxend,flyStart,      &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom left of flaw
       call section_velocity(lxstart,flyStart,lzstart,flxStart,     &
            &    lyend,frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,flyStart,frSz,lxend,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(flxStart,flyStart,lzstart,lxend,lyend,     &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       
    CASE (12) ! lower x & upper y bounds of flaw
       !bottom left of flaw
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyEnd,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,lystart,frSz,flxStart,flyEnd,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,ylb,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right under flaw
       call section_velocity(flxStart,lystart,lzstart,lxend,flyEnd,     &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper left
       call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(flxStart,flyEnd,frSz,lxend,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,0,xub,1,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       
    CASE (13) ! lower x & y & upper y bounds of flaw
       !lower left
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart, &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom before flaw
       call section_velocity(flxStart,lystart,lzstart,lxend,flyStart,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top before flaw
       call section_velocity(flxStart,lystart,frSz,lxend,flyStart,      &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom left of flaw
       call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,flyStart,frSz,flxStart,lyend,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(flxStart,flyStart,lzstart,lxend,flyEnd,    &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper left
       call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(flxStart,flyEnd,frSz,lxend,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,0,xub,1,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (20) ! upper x bound of flaw
       ! left under flaw
       call section_velocity(lxstart,lystart,lzstart,flxEnd,lyend,      &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,        &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,lystart,frSz,lxend,lyend,lzend,     &
            &    flx,fly,numz,dtodsp,1,xub,ylb,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (21) ! upper x & lower y bounds of flaw
       !bottom left before of flaw
       call section_velocity(lxstart,lystart,lzstart,flxEnd,flyStart,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left before of flaw
       call section_velocity(lxstart,lystart,frSz,flxEnd,flyStart,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !lower right
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(lxstart,flyStart,lzstart,flxEnd,lyend,     &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,flyStart,frSz,lxend,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,1,xub,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       
    CASE (22) ! upper x & y bounds of flaw
       !bottom under flaw
       call section_velocity(lxstart,lystart,lzstart,flxEnd,flyEnd,     &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,lystart,frSz,lxend,flyEnd,lzend,    &
            &    flx,fly,numz,dtodsp,1,xub,ylb,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(lxstart,flyEnd,lzstart,flxEnd,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(lxstart,flyEnd,frSz,flxEnd,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,xlb,0,1,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper right
       call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,lzend,   &
            &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
 
     CASE (23) ! lower y & upper x & y bounds of flaw
       !bottom left before of flaw
       call section_velocity(lxstart,lystart,lzstart,flxEnd,flyStart,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left before of flaw
       call section_velocity(lxstart,lystart,frSz,flxEnd,flyStart,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !lower right
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(lxstart,flyStart,lzstart,flxEnd,flyEnd,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,flyStart,lzstart,lxend,flyEnd,      &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,flyStart,frSz,lxend,flyEnd,         &
            &    lzend,flx,fly,numz,dtodsp,1,xub,0,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(lxstart,flyEnd,lzstart,flxEnd,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(lxstart,flyEnd,frSz,flxEnd,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,xlb,0,1,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper right
       call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    CASE (30) ! upper & lower x bounds of flaw
       !bottom left of flaw
       call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,lystart,frSz,flxStart,lyend,       &
            &    lzend,flx,fly,numz,dtodsp, xlb,1,ylb,yub,0,zub,  &
            &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right under flaw
       call section_velocity(flxStart,lystart,lzstart,flxEnd,lyend,     &
            &    frSz,flx,fly,numz,dtodsp,0,0,ylb,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,frSz,   &
            &    flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,lystart,frSz,lxend,lyend,lzend,     &
            &    flx,fly,numz,dtodsp,1,xub,ylb,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (31) ! lower x & y & upper x bounds of flaw
       !lower left
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart, &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom before flaw
       call section_velocity(flxStart,lystart,lzstart,flxEnd,flyStart,  &
            &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top before flaw
       call section_velocity(flxStart,lystart,frSz,flxEnd,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,0,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !lower right
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom left of flaw
       call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,flyStart,frSz,flxStart,lyend,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(flxStart,flyStart,lzstart,flxEnd,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,flyStart,frSz,lxend,lyend,lzend,    &
            &    flx,fly,numz,dtodsp,1,xub,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (32) ! upper x & y & lower x bounds of flaw
       !bottom left of flaw
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyEnd,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,lystart,frSz,flxStart,flyEnd,      &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,ylb,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !under flaw
       call section_velocity(flxStart,lystart,lzstart,flxEnd,flyEnd,    &
            &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !lbottom right of flaw
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,       &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,lystart,frSz,lxend,flyEnd,lzend,    &
            &    flx,fly,numz,dtodsp,0,xub,ylb,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper left
       call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(flxStart,flyEnd,lzstart,flxEnd,lyend,      &
            &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(flxStart,flyEnd,frSz,flxEnd,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,0,0,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper right
       call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,lzend,   &
            &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)


     CASE (33) ! upper & lower x & y bound of flaw
       !lower left
       call section_velocity(lxstart,lystart,lzstart,flxStart,flyStart, &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom before flaw
       call section_velocity(flxStart,lystart,lzstart,flxEnd,flyStart,  &
            &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top before flaw
       call section_velocity(flxStart,lystart,frSz,flxEnd,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,0,ylb,1,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !lower right
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyStart,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom left of flaw
       call section_velocity(lxstart,flyStart,lzstart,flxStart,flyEnd,  &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top left of flaw
       call section_velocity(lxstart,flyStart,frSz,flxStart,flyEnd,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,1,0,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom under flaw
       call section_velocity(flxStart,flyStart,lzstart,flxEnd,flyEnd,   &
            &    frSz,flx,fly,numz,dtodsp,0,0,0,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom right of flaw
       call section_velocity(flxEnd,flyStart,lzstart,lxend,flyEnd,      &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top right of flaw
       call section_velocity(flxEnd,flyStart,frSz,lxend,flyEnd,     &
            &    lzend,flx,fly,numz,dtodsp,1,xub,0,0,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper left
       call section_velocity(lxstart,flyEnd,lzstart,flxStart,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !bottom after flaw
       call section_velocity(flxStart,flyEnd,lzstart,flxEnd,lyend,      &
            &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !top after flaw
       call section_velocity(flxStart,flyEnd,frSz,flxEnd,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,0,0,0,yub,0,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !upper right
       call section_velocity(flxEnd,flyEnd,lzstart,lxend,lyend,     &
            &    lzend,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    end SELECT
  end Subroutine update_flaw_vel

!========================================================================

Subroutine update_flaw_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
     &   frSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
     &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
     &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
     &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

  Implicit none

  Integer, intent(in) :: flx,fly,numz,frSz
  Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
  Double Precision, intent(in) :: dtods,l2m,lambda,mu
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fvx,fvy,fvz

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
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

      CASE (2)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,1,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,1,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

      CASE (4)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)


     CASE (10)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,1,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

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
             &    dtods,l2m,lambda,mu,0,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,1,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (12)
        !bottom left of flaw
        call section_stress(lxstart,lystart,lzstart, &
             &    flxStart,flyEnd,frSz,flx,fly,numz, &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,lystart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,1,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
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
             &    dtods,l2m,lambda,mu,0,xub,1,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,0,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom left of flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxStart,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top left of flaw
        call section_stress(lxstart,flyStart,frSz,  &
             &    flxStart,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,1,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
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
             &    dtods,l2m,lambda,mu,0,xub,1,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (20)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,ylb,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (22)
        !under flaw
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,1,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,1,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !lower right
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyStart,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(lxstart,flyStart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom after flaw
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top after flaw
        call section_stress(lxstart,flyEnd,frSz,  &
             &    flxEnd,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,1,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,1,ylb,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,ylb,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,0,0,ylb,1,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,1,0,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,0,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,1,ylb,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,lystart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,ylb,0,0,zub, &
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
             &    dtods,l2m,lambda,mu,0,0,1,yub,0,zub, &
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
             &    dtods,l2m,lambda,mu,0,0,ylb,1,0,zub, &
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
             &    dtods,l2m,lambda,mu,xlb,1,0,0,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !under flaw
        call section_stress(flxStart,flyStart,lzstart,  &
             &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !bottom right of flaw
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,0, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !top right of flaw
        call section_stress(flxEnd,flyStart,frSz,  &
             &    lxend,flyEnd,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,1,xub,0,0,0,zub, &
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
             &    dtods,l2m,lambda,mu,0,0,1,yub,0,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !upper right
        call section_stress(flxEnd,flyEnd,lzstart,  &
             &    lxend,lyend,lzend,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,zub, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     end SELECT

   end Subroutine update_flaw_stress
