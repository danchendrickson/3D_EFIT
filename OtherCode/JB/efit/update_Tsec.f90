Subroutine T_velocity(lx,ly,numx,numz,startx,starty,endy,dtodsp, &
     &     ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,  &
     &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
     &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
     &     rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
     &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
     &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)

  Implicit none

 
  Integer, intent(in) :: lx,ly,numx,numz,startx,starty,endy
  Integer, intent(in) :: ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE
  Integer, intent(in) :: ltwyS,ltwyE,twyS,twyE,twzS,twzE
  Integer, intent(in) :: lttyS,lttyE,ttyS,ttyE,ttzS,ttzE
  Integer, intent(in) :: rSz,flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
  Integer, intent(in) :: numS,a,b,c
  Double Precision, intent(in) :: dtodsp
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout) :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)      &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
  Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds

  Integer :: wlb,wub,tlb,tub




  ! Have to be careful of what part of the stringer the CPU is on
  !  if (starty < twyS) then ! definitely have to update flage before web
  !     if (endy < twyS) then ! just flange, don't worry about web
  !     else if (endy > twyE) then ! have to do before, in and after web
  !     else ! haev to do before web and in web
  !     end if
  !  else if (starty > twyE) then ! only have to update flange after web
  !  else ! update in web definitely
  !     if (endy > twyE) then ! update in and after web
  !     else ! just in web
  !     end if
  !  end if

 
  if (starty < twyS) then ! definitely have to update flage before web
     if (endy < twyS) then ! just flange, don't worry about web
        if (flaw == 0) then
           ! Flange before the web
           call section_velocity(1,ltfyS,tfzS,lx,ltfyE,tfzE,lx,ly,numz,  &
                &    dtodsp,xlb,xub,ylb,0,zlb,1,vx,vy,vz,   &
                &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
              CASE (0)
                 call  update_flaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS,  &
                      &   lx,ltfyE,tfzE,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call  update_sflaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS, &
                      &   lx,ltfyE,tfzE,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,      &
                      &   a,b,c,starty)
              end SELECT
        end if
     else if (endy > twyE) then ! have to do before, in and after web
        if (flaw == 0) then
           ! Flange before the web
           call section_velocity(1,ltfyS,tfzS,lx,ltwyS,tfzE,lx,ly,numz,  &
                &    dtodsp,xlb,xub,ylb,0,zlb,1,vx,vy,vz,   &
                &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
              CASE (0)
                 call  update_flaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS,  &
                      &   lx,ltwyS,tfzE,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call  update_sflaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS, &
                      &   lx,ltwyS,tfzE,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,  &
                      &   a,b,c,starty)
              end SELECT                 
        end if                
        ! Flange under the web
        call section_velocity(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,0,zlb,0,vx,vy,vz,    &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Flange after the web
        call section_velocity(1,ltwyE,tfzS,lx,ltfyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,yub,zlb,1,vx,vy,vz,      &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     else ! have to do before web and in web
        if (flaw == 0) then
           ! Flange before the web
           call section_velocity(1,ltfyS,tfzS,lx,ltwyS,tfzE,lx,ly,numz,  &
                &    dtodsp,xlb,xub,ylb,0,zlb,1,vx,vy,vz,   &
                &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
           CASE (0)
              call  update_flaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS,     &
                   &   lx,ltwyS,tfzE,  &
                   &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                   &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                   &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
           CASE (1)
              call  update_sflaw_vel(lx,ly,numz,dtodsp,1,ltfyS,tfzS,    &
                   &   lx,ltwyS,tfzE,  &
                   &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
                   &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                   &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c,starty)
           end SELECT
        end if                
        ! Flange under the web
        call section_velocity(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,0,zlb,0,vx,vy,vz,    &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  else if (starty > twyE) then ! only have to update flange after web
     ! Flange after the web
     call section_velocity(1,ltfyS,tfzS,lx,ltfyE,tfzE,lx,ly,numz,dtodsp,  &
          &    xlb,xub,0,yub,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
  else ! update in web definitely
     if (endy > twyE) then ! update in and after web
        ! Flange under the web
        call section_velocity(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,0,zlb,0,vx,vy,vz,    &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Flange after the web
        call section_velocity(1,ltwyE,tfzS,lx,ltfyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,yub,zlb,1,vx,vy,vz,      &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     else ! just in web
        ! Flange under the web
        call section_velocity(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,0,zlb,0,vx,vy,vz,    &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  end if

  if (ltwyS > 0) then
     if (ltwyS == twyS-starty+1) then
        wlb = 1
     else
        wlb = 0
     end if
     if (ltwyE == twyE-starty+1) then
        wub = 1
     else
        wub = 0
     end if
     
     ! Web
     call section_velocity(1,ltwyS,twzS,lx,ltwyE,twzE,lx,ly,numz,dtodsp,  &
          &    xlb,xub,wlb,wub,0,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)

     if (lttyS > 0) then
        if (lttyS == ttyS-starty+1) then
           tlb = 1
        else
           tlb = 0
        end if
        if (lttyE == ttyE-starty+1) then
           tub = 1
        else
           tub = 0
        end if
        
        ! Top over the web
        call section_velocity(1,ltwyS,ttzS,lx,ltwyE,ttzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,0,0,1,vx,vy,vz,      &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Top before the web
        call section_velocity(1,lttyS,ttzS,lx,ltwyS,ttzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,tlb,0,1,1,vx,vy,vz,        &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Top after the web
        call section_velocity(1,ltwyE,ttzS,lx,lttyE,ttzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,0,tub,1,1,vx,vy,vz,        &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  else
     if (lttyS > 0) then
        if (lttyS == ttyS-starty+1) then
           tlb = 1
        else
           tlb = 0
        end if
        if (lttyE == ttyE-starty+1) then
           tub = 1
        else
           tub = 0
        end if
        
        ! Top before the web or after
        call section_velocity(1,lttyS,ttzS,lx,lttyE,ttzE,lx,ly,numz,  &
             &    dtodsp,xlb,xub,tlb,tub,1,1,vx,vy,vz,      &
             &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  end if

end Subroutine T_velocity


Subroutine T_stress(lx,ly,numx,numz,startx,starty,endy,dtods,l2m, &
     &     lambda,mu,ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,  &
     &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
     &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
     &     rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
     &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
     &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)

  Implicit none

 
  Integer, intent(in) :: lx,ly,numx,numz,startx,starty,endy
  Integer, intent(in) :: ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE
  Integer, intent(in) :: ltwyS,ltwyE,twyS,twyE,twzS,twzE
  Integer, intent(in) :: lttyS,lttyE,ttyS,ttyE,ttzS,ttzE
  Integer, intent(in) :: rSz,flxStart,flxEnd,flyStart,flyEnd,flaw
  Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
  Integer, intent(in) :: numS,a,b,c
  Double Precision, intent(in) :: dtods,l2m,lambda,mu
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout) :: vx,vy,vz
  Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)      &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
  Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds

  Integer :: wlb,wub,tlb,tub

  ! Have to be careful of what part of the stringer the CPU is on
  !  if (starty < twyS) then ! definitely have to update flage before web
  !     if (endy < twyS) then ! just flange, don't worry about web
  !     else if (endy > twyE) then ! have to do before, in and after web
  !     else ! haev to do before web and in web
  !     end if
  !  else if (starty > twyE) then ! only have to update flange after web
  !  else ! update in web definitely
  !     if (endy > twyE) then ! update in and after web
  !     else ! just in web
  !     end if
  !  end if

 
  if (starty < twyS) then ! definitely have to update flage before web
     if (endy < twyS) then ! just flange, don't worry about web  
        if (flaw == 0) then
           ! Flange before the web
           call  section_stress(1,ltfyS,tfzS,lx,ltfyE,tfzE,lx,ly,numz,  &
                &    dtods,l2m,lambda,mu, &
                &    xlb,xub,ylb,0,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
              CASE (0)
                 call update_flaw_stress(lx,ly,numz,dtods,l2m,lambda,mu,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltfyE,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,  &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call update_sflaw_stress(lx,ly,numz,dtods,l2m,lambda,mu, &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltfyE,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c)
              end SELECT
        end if
     else if (endy > twyE) then ! have to do before, in and after web
        if (flaw == 0) then
           ! Flange before the web
           call  section_stress(1,ltfyS,tfzS,lx,ltwyS,tfzE,lx,ly,numz,  &
                &    dtods,l2m,lambda,mu, &
                &    xlb,xub,ylb,0,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
              CASE (0)
                 call update_flaw_stress(lx,ly,numz,dtods,l2m,lambda,mu,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltwyS,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call update_sflaw_stress(lx,ly,numz,dtods,l2m,lambda,mu, &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltwyS,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c)
              end SELECT
        end if
        ! Flange under the web
        call  section_stress(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,0,zlb,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Flange after the web
        call  section_stress(1,ltwyE,tfzS,lx,ltfyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,yub,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     else ! have to do before and in web
        if (flaw == 0) then
           ! Flange before the web
           call  section_stress(1,ltfyS,tfzS,lx,ltwyS,tfzE,lx,ly,numz,  &
                &    dtods,l2m,lambda,mu, &
                &    xlb,xub,ylb,0,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           SELECT CASE (numS)
              CASE (0)
                 call update_flaw_stress(lx,ly,numz,dtods,l2m,lambda,mu,  &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltwyS,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call update_sflaw_stress(lx,ly,numz,dtods,l2m,lambda,mu, &
                      &   rSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
                      &   1,ltfyS,tfzS,lx,ltwyS,tfzE,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c)
              end SELECT
        end if
        ! Flange under the web
        call  section_stress(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,0,zlb,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Flange after the web
        call  section_stress(1,ltwyE,tfzS,lx,ltfyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,yub,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  else if (starty > twyE) then ! only have to update flangee after web
     ! Flange after the web
     call  section_stress(1,ltfyS,tfzS,lx,ltfyE,tfzE,lx,ly,numz,        &
          &    dtods,l2m,lambda,mu, &
          &    xlb,xub,0,yub,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
  else ! update in web definitely
     if (endy > twyE) then ! update in and after web
        ! Flange under the web
        call  section_stress(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,0,zlb,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Flange after the web
        call  section_stress(1,ltwyE,tfzS,lx,ltfyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,yub,zlb,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     else ! just in web
        ! Flange under the web
        call  section_stress(1,ltwyS,tfzS,lx,ltwyE,tfzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,0,zlb,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  end if

  if (ltwyS > 0) then
     if (ltwyS == twyS-starty+1) then
        wlb = 1
     else
        wlb = 0
     end if
     if (ltwyE == twyE-starty+1) then
        wub = 1
     else
        wub = 0
     end if
     
     ! Web
     call  section_stress(1,ltwyS,twzS,lx,ltwyE,twzE,lx,ly,numz,        &
          &    dtods,l2m,lambda,mu, &
          &    xlb,xub,wlb,wub,0,0,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     
     if (lttyS > 0) then
        if (lttyS == ttyS-starty+1) then
           tlb = 1
        else
           tlb = 0
        end if
        if (lttyE == ttyE-starty+1) then
           tub = 1
        else
           tub = 0
        end if
        
        ! Top over the web
        call  section_stress(1,ltwyS,ttzS,lx,ltwyE,ttzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,0,0,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Top before the web
        call  section_stress(1,lttyS,ttzS,lx,ltwyS,ttzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,tlb,0,1,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        ! Top after the web
        call  section_stress(1,ltwyE,ttzS,lx,lttyE,ttzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,0,tub,1,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  else
     if (lttyS > 0) then
        if (lttyS == ttyS-starty+1) then
           tlb = 1
        else
           tlb = 0
        end if
        if (lttyE == ttyE-starty+1) then
           tub = 1
        else
           tub = 0
        end if
        
        ! Top before or after web
        call  section_stress(1,lttyS,ttzS,lx,lttyE,ttzE,lx,ly,numz,     &
             &    dtods,l2m,lambda,mu, &
             &    xlb,xub,tlb,0,1,1,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
     end if
  end if

end Subroutine T_stress
