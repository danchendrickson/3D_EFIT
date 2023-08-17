Subroutine update_mass_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,   &
     &   lxstart,lystart,lzstart,lxend,lyend,lzend,  &
     &   frSz,frEz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
     &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
     &   dexStart,dexEnd,deyStart,deyEnd,dxlb,dxub,dylb,dyub,delam,  &
     &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

  Implicit none

 
  Integer, intent(in) :: flx,fly,numz,frSz,frEz
  Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
  Integer, intent(in) :: fxlb,fxub,fylb,fyub
  Integer, intent(in) :: dexStart,dexEnd,deyStart,deyEnd
  Integer, intent(in) :: dxlb,dxub,dylb,dyub,delam
  Double Precision, intent(in) :: dtodsp,dtodsp2,dden,rden
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fvx,fvy,fvz
  Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)    &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz

  Integer :: ix,iy,iz,tix,tiy,tiz

  !flx - lx, local number in x-dir

!  write(*,*) 'F1 updating flaw vel'

  ! This is called if the current CPU has a flawed region on it

  ! Update any part of the cpu not in flawed region normalfly.
  ! Basicalfly reverse mapping of the region.
  ! flaw: x lower +10, x upper +20, y lower +1, y upper +2
  !      possibles - 1,2,3,4,10,11,12,13,20,21,22,23,30,31,32,33
  
!  write(*,*) 'flag 3 ',frEz

  SELECT CASE (flaw)
     CASE (1) ! lower y bound of flaw
        !plate before mass
        call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,flyStart,lzstart,lxend,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,xub,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,flyStart,lzstart,lxend,lyend,  &
                &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,flyStart,lxend,lyend,frSz+1,   &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,fxub,0,fyub,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,flyStart,frSz+1,lxEnd,lyend,   &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

     CASE (2) ! upper y bound of flaw
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,lystart,lzstart,lxend,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,xub,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,lystart,lzstart,lxend,flyEnd,  &
                &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,lystart,lxend,flyEnd,frSz+1, &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,fxub,fylb,0,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,lystart,frSz+1,lxend,flyEnd,   &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (3) ! lower and upper y bounds of flaw
        !plate before mass
        call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,0,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,flyStart,lzstart,lxend,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,flyStart,lzstart,lxend,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,flyStart,lxend,flyEnd,frSz+1,  &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,fxub,0,0,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,flyStart,frSz+1,lxEnd,flyEnd,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (4) ! flaw completely covers
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,lystart,lzstart,lxend,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,lystart,lzstart,lxend,lyend,   &
                &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,lystart,lxend,lyend,frSz+1,  &
                &   flx,fly,numz,dtodsp,dden,rden,fxlb,fxub,fylb,fyub,0,0,&
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass over all
           call section_velocity(lxstart,lystart,frSz+1,lxEnd,lyend,    &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

     CASE (10) ! lower x bound of flaw
        !plate left of mass
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,lystart,lzstart,lxend,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(flxStart,lystart,lzstart,lxend,lyend,  &
                &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,lystart,lxend,lyend,frSz+1, &
                &    flx,fly,numz,dtodsp,dden,rden,0,fxub,fylb,fyub,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,lystart,frSz+1,lxend,lyend,   &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

    CASE (11) ! lower x & y bounds of flaw
        !plate before mass
       call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !plate left of mass
       call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,flyStart,lzstart,lxend,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,xub,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !under mass
           call section_velocity(flxStart,flyStart,lzstart,lxend,lyend, &
                &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,flyStart,lxend,lyend,frSz+1,  &
                &   flx,fly,numz,dtodsp,dden,rden,0,fxub,0,fyub,0,0,  &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,flyStart,frSz+1,lxend,lyend,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

    CASE (12) ! lower x & upper y bounds of flaw
       !plate left of mass
       call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,lystart,lzstart,lxend,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,xub,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate right under mass
           call section_velocity(flxStart,lystart,lzstart,lxend,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,lystart,lxend,flyEnd,frSz+1, &
                &    flx,fly,numz,dtodsp,dden,rden,0,fxub,fylb,0,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,lystart,frSz+1,lxend,flyEnd,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate after mass
       call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
    CASE (13) ! lower x & y & upper y bounds of flaw
       !plate before mass
       call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !plate left of mass
       call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,flyStart,lzstart,lxend,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,xub,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !under mass
           call section_velocity(flxStart,flyStart,lzstart,lxend,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,flyStart,lxend,flyEnd,      &
                &    frSz+1,flx,fly,numz,dtodsp,dden,rden,0,fxub,0,0,0,0, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,flyStart,frSz+1,lxend,flyEnd, &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate after mass
       call section_velocity(flxStart,flyEnd,lzstart,lxend,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (20) ! upper x bound of flaw
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,lystart,lzstart,flxEnd,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,0,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,lystart,lzstart,flxEnd,lyend,  &
                &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,lystart,flxEnd,lyend,frSz+1, &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,0,fylb,fyub,0,0, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,lystart,frSz+1,flxEnd,lyend,   &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,  &
                &    zub,fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,frSz,   &
            &    flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (21) ! upper x & lower y bounds of flaw
       !plate before mass
       call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,flyStart,lzstart,flxEnd,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,0,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,flyStart,lzstart,flxEnd,lyend, &
                &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,flyStart,flxEnd,lyend,frSz+1,  &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,0,0,fyub,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,flyStart,frSz+1,flxEnd,lyend,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,frSz,  &
            &    flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    CASE (22) ! upper x & y bounds of flaw
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,lystart,lzstart,flxEnd,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,0,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,lystart,lzstart,flxEnd,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,lystart,flxEnd,flyEnd,frSz+1,  &
                &    flx,fly,numz,dtodsp,dden,rden,fxlb,0,fylb,0,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,lystart,frSz+1,flxEnd,flyEnd,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !plate after mass
       call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (23) ! lower y & upper x & y bounds of flaw
       !plate before mass
       call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   lxstart,flyStart,lzstart,flxEnd,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   xlb,0,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(lxstart,flyStart,lzstart,flxEnd,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,xlb,0,0,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(lxstart,flyStart,flxEnd,flyEnd,      &
                &    frSz+1,flx,fly,numz,dtodsp,dden,rden,fxlb,0,0,0,0,0, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(lxstart,flyStart,frSz+1,flxEnd,flyEnd, &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,flyStart,lzstart,lxend,flyEnd,  &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,0,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !plate after mass
       call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    CASE (30) ! upper & lower x bounds of flaw
        !plate left of mass
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,lystart,lzstart,flxEnd,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,0,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(flxStart,lystart,lzstart,flxEnd,lyend, &
                &    frSz,flx,fly,numz,dtodsp,0,0,ylb,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,lystart,flxEnd,lyend,frSz+1, &
                &    flx,fly,numz,dtodsp,dden,rden,0,0,fylb,fyub,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,lystart,frSz+1,flxEnd,lyend,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,lystart,lzstart,lxend,lyend,    &
            &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

    CASE (31) ! lower x & y & upper x bounds of flaw
       !plate before mass
       call section_velocity(lxstart,lystart,lzstart,lxend,flyStart,    &
            &    frSz,flx,fly,numz,dtodsp,xlb,xub,ylb,0,zlb,zub,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       !plate left of mass
       call section_velocity(lxstart,flyStart,lzstart,flxStart,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,xlb,0,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,flyStart,lzstart,flxEnd,lyend,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,0,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(flxStart,flyStart,lzstart,flxEnd,lyend,  &
                &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,flyStart,flxEnd,lyend,  &
                &    frSz+1,flx,fly,numz,dtodsp,dden,rden,0,0,0,fyub,0,0, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,flyStart,frSz+1,flxEnd,lyend, &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
       !plate right of mass
       call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,   &
            &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
            &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       
    CASE (32) ! upper x & y & lower x bounds of flaw
        !plate left of mass
        call section_velocity(lxstart,lystart,lzstart,flxStart,flyEnd,  &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,0,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,lystart,lzstart,flxEnd,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,0,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(flxStart,lystart,lzstart,flxEnd,flyEnd,  &
                &    frSz,flx,fly,numz,dtodsp,0,0,ylb,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,lystart,flxEnd,flyEnd,  &
                &    frSz+1,flx,fly,numz,dtodsp,dden,rden,0,0,fylb,0,0,0, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,lystart,frSz+1,flxEnd,flyEnd, &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_velocity(flxEnd,lystart,lzstart,lxend,flyEnd,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate after mass
        call section_velocity(lxstart,flyEnd,lzstart,lxend,lyend,frSz,  &
             &    flx,fly,numz,dtodsp,xlb,xub,0,yub,zlb,zub,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
     CASE (33) ! upper & lower x & y bound of flaw
        !plate left of mass
        call section_velocity(lxstart,lystart,lzstart,flxStart,lyend,   &
             &    frSz,flx,fly,numz,dtodsp,xlb,0,ylb,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate before mass
        call section_velocity(flxStart,lystart,lzstart,lxend,flyStart,  &
             &    frSz,flx,fly,numz,dtodsp,0,xub,ylb,0,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
                &   flxStart,flyStart,lzstart,flxEnd,flyEnd,lzend,  &
                &   frSz,frEz,dexStart,dexEnd,deyStart,deyEnd,delam,  &
                &   0,0,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        else
           !plate under mass
           call section_velocity(flxStart,flyStart,lzstart,flxEnd,flyEnd, &
                &    frSz,flx,fly,numz,dtodsp,0,0,0,0,zlb,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_velocity(flxStart,flyStart,flxEnd,flyEnd,     &
                &    frSz+1,flx,fly,numz,dtodsp,dden,rden,0,0,0,0,0,0,  &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_velocity(flxStart,flyStart,frSz+1,flxEnd,flyEnd,  &
                &    frEz,flx,fly,numz,dtodsp2,fxlb,fxub,fylb,fyub,0,zub, &
                &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_velocity(flxStart,flyEnd,lzstart,flxEnd,lyend,     &
             &    frSz,flx,fly,numz,dtodsp,0,0,0,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate right of mass
        call section_velocity(flxEnd,flyStart,lzstart,lxend,lyend,      &
             &    frSz,flx,fly,numz,dtodsp,0,xub,0,yub,zlb,1,  &
             &    fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        
    end SELECT
  end Subroutine update_mass_vel

!==========================================================================

Subroutine update_mass_stress(flx,fly,numz,dtods,l2m,l2m2,lambda,  &
     &   lambda2,mu,mu2,frSz,frEz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
     &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
     &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
     &   dexStart,dexEnd,deyStart,deyEnd,dxlb,dxub,dylb,dyub,delam,  &
     &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

  Implicit none

  Integer, intent(in) :: flx,fly,numz,frSz,frEz
  Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
  Integer, intent(in) :: dexStart,dexEnd,deyStart,deyEnd
  Integer, intent(in) :: dxlb,dxub,dylb,dyub,delam
  Double Precision, intent(in) :: dtods,l2m,lambda,mu,l2m2,lambda2,mu2
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
        !plate before mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,  &
                &   lambda,lambda2,mu,mu2,frSz,frEz,   &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,flyStart,lzstart,lxend,lyend,lzend,   &
                &   xlb,xub,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(lxstart,flyStart,lzstart,  &
                &    lxend,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,flyStart,  &
                &    lxend,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,fxub,0,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,flyStart,frSz+1,  &
                &    lxend,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

      CASE (2)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,  &
                &   lambda,lambda2,mu,mu2,frSz,frEz,        &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,lystart,lzstart,lxend,flyEnd,lzend,   &
                &   xlb,xub,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !under mass
           call section_stress(lxstart,lystart,lzstart,  &
                &    lxend,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,lystart,  &
                &    lxend,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,fxub,fylb,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,lystart,frSz+1,  &
                &    lxend,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (3)
        !plate before mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,flyStart,lzstart,lxend,flyEnd,lzend,   &
                &   xlb,xub,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(lxstart,flyStart,lzstart,  &
                &    lxend,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,0,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,flyStart,  &
                &    lxend,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,fxub,0,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,flyStart,frSz+1,  &
                &    lxend,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (4)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,     &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
                &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !under mass
           call section_stress(lxstart,lystart,lzstart,  &
                &    lxend,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,lystart,  &
                &    lxend,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,fxub,fylb,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,lystart,frSz+1,  &
                &    lxend,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

     CASE (10)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,lystart,lzstart,lxend,lyend,lzend,   &
                &   0,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,lystart,lzstart,  &
                &    lxend,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,lystart,  &
                &    lxend,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,fxub,fylb,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,lystart,frSz+1,  &
                &    lxend,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

     CASE (11)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate before mass
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,flyStart,lzstart,lxend,lyend,lzend,   &
                &   0,xub,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,flyStart,lzstart,  &
                &    lxend,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,flyStart,  &
                &    lxend,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,fxub,0,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,flyStart,frSz+1,  &
                &    lxend,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if

     CASE (12)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,lystart,lzstart,lxend,flyEnd,lzend,   &
                &   0,xub,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,lystart,lzstart,  &
                &    lxend,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,lystart,  &
                &    lxend,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,fxub,fylb,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,lystart,frSz+1,  &
                &    lxend,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (13)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate before mass
        call section_stress(flxStart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,flyStart,lzstart,lxend,flyEnd,lzend,   &
                &   0,xub,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,flyStart,lzstart,  &
                &    lxend,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,flyStart,  &
                &    lxend,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,fxub,0,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,flyStart,frSz+1,  &
                &    lxend,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (20)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,lystart,lzstart,flxEnd,lyend,lzend,   &
                &   xlb,0,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !under mass
           call section_stress(lxstart,lystart,lzstart,  &
                &    flxEnd,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,lystart,  &
                &    flxEnd,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,0,fylb,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,lystart,frSz+1,  &
                &    flxEnd,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (21)
        !plate before mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,flyStart,lzstart,flxEnd,lyend,lzend,   &
                &   xlb,0,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(lxstart,flyStart,lzstart,  &
                &    flxEnd,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,0,0,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,flyStart,  &
                &    flxEnd,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,0,0,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,flyStart,frSz+1,  &
                &    flxEnd,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (22)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,lystart,lzstart,flxEnd,flyEnd,lzend,   &
                &   xlb,0,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !under mass
           call section_stress(lxstart,lystart,lzstart,  &
                &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,lystart,  &
                &    flxEnd,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,0,fylb,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,lystart,frSz+1,  &
                &    flxEnd,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate after mass
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (23)
        !plate before mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    lxend,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &   
                &   lambda,lambda2,mu,mu2,frSz,frEz,dexStart,  &
                &   dexEnd,deyStart,deyEnd,delam, &
                &   lxstart,flyStart,lzstart,flxEnd,flyEnd,lzend,   &
                &   xlb,0,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(lxstart,flyStart,lzstart,  &
                &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,0,0,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(lxstart,flyStart,  &
                &    flxEnd,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,fxlb,0,0,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(lxstart,flyStart,frSz+1,  &
                &    flxEnd,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,flyStart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,0,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate after mass
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (30)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,  &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,lystart,lzstart,flxEnd,lyend,lzend,   &
                &   0,0,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,lystart,lzstart,  &
                &    flxEnd,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,0,ylb,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,lystart,  &
                &    flxEnd,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,0,fylb,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,lystart,frSz+1,  &
                &    flxEnd,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (31)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate before mass
        call section_stress(flxStart,lystart,lzstart,  &
             &    flyEnd,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,  &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,flyStart,lzstart,flxEnd,lyend,lzend,   &
                &   0,0,0,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,flyStart,lzstart,  &
                &    flxEnd,lyend,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,flyStart,  &
                &    flxEnd,lyend,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,0,0,fyub,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,flyStart,frSz+1,  &
                &    flxEnd,lyend,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (32)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,  &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,lystart,lzstart,flxEnd,flyEnd,lzend,   &
                &   0,0,ylb,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,lystart,lzstart,  &
                &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,lystart,  &
                &    flxEnd,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,0,fylb,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,lystart,frSz+1,  &
                &    flxEnd,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,flyEnd,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate after mass
        call section_stress(lxstart,flyEnd,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,xub,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     CASE (33)
        !plate left of mass
        call section_stress(lxstart,lystart,lzstart,  &
             &    flxStart,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,xlb,0,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate before mass
        call section_stress(flxStart,lystart,lzstart,  &
             &    flxEnd,flyStart,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,ylb,0,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        if (delam > 0) then
           call update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,    &
                &   lambda,lambda2,mu,mu2,frSz,frEz,  &
                &   dexStart,dexEnd,deyStart,deyEnd,delam, &
                &   flxStart,flyStart,lzstart,flxEnd,flyEnd,lzend,   &
                &   0,0,0,0,zlb,zub,fxlb,fxub,fylb,fyub,   &
                &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)        
        else
           !plate under mass
           call section_stress(flxStart,flyStart,lzstart,  &
                &    flxEnd,flyEnd,frSz,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,0,0,0,0,zlb,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !interface under mass
           call interface_stress(flxStart,flyStart,  &
                &    flxEnd,flyEnd,frSz+1,flx,fly,numz,  &
                &    dtods,l2m,lambda,mu,mu2,0,0,0,0,0,0, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
           !mass
           call section_stress(flxStart,flyStart,frSz+1,  &
                &    flxEnd,flyEnd,frEz,flx,fly,numz,  &
                &    dtods,l2m2,lambda2,mu2,fxlb,fxub,fylb,fyub,0,zub, &
                &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        end if
        !plate after mass
        call section_stress(flxStart,flyEnd,lzstart,  &
             &    flxEnd,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,0,0,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
        !plate right of mass
        call section_stress(flxEnd,lystart,lzstart,  &
             &    lxend,lyend,frSz,flx,fly,numz,  &
             &    dtods,l2m,lambda,mu,0,xub,ylb,yub,zlb,1, &
             &	  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)

     end SELECT

   end Subroutine update_mass_stress
