Program code
! Adds T structure
! Include bounds but not individual material properties 5-14-2008
! Adding mpi code parts 2-27-2008
! With input files
! Adding individual cell parameters in order to build up samples by cells
! Second try
! Conversion of efit2d.f90 to 3d JPB 11-2007
! Conversion matlab to f90 authors Jill Bingham and Chris Bording
! Original matlab code written by Kevin Rudd - NDE
! the college of william and mary - 2006
  USE exchange

  Implicit none

  include "mpif.h"
  include "silo.inc"

  Interface
     Subroutine tukeywin(w,N,alpha)
       Integer,Intent(in) :: N
       Double Precision ,Intent(in) :: alpha
       Double Precision ,Dimension(1:N),Intent(out) :: w(:) 
     end Subroutine tukeywin
     
     !Input parameter functions
     Subroutine read_model(alpha,plate,numR,delam,numS,plotevery,   &
        &   silo_on,siloSpace,fmax,SimulationTime,Pdim0,Pdim1)
       Double Precision, Intent(out) :: alpha
       Double Precision, Intent(out) :: fmax
       Double Precision, Intent(out) :: SimulationTime
       Integer, Intent(out) :: plate
       Integer, Intent(out) :: numR
       Integer, Intent(out) :: delam
       Integer, Intent(out) :: numS
       Integer, Intent(out) :: plotevery
       Integer, Intent(out) :: silo_on
       Integer, Intent(out) :: siloSpace
       Integer, Intent(out) :: Pdim0,Pdim1
     end Subroutine read_model

     Subroutine read_dmaterials(cmax,cmin,spacethickness,spacelength,   &
        &   spacewidth,dden,dcl,dcs,minz)
       Double Precision , Intent(out) :: cmax, cmin
       Double Precision , Intent(out) :: spacethickness        
       Double Precision , Intent(out) :: spacelength
       Double Precision , Intent(out) :: spacewidth
       Double Precision , Intent(out) :: dden, dcl, dcs
       Double Precision , Intent(out) :: minz
     end Subroutine read_dmaterials

     Subroutine read_T_region(tfyStart,tfyEnd,tfzStart,tfzEnd,twyStart, &
        &   twyEnd,twzStart,twzEnd,ttyStart,ttyEnd,ttzStart,ttzEnd)
       Double Precision, Intent(out) :: tfyStart
       Double Precision, Intent(out) :: tfyEnd
       Double Precision, Intent(out) :: tfzStart
       Double Precision, Intent(out) :: tfzEnd
       Double Precision, Intent(out) :: twyStart
       Double Precision, Intent(out) :: twyEnd
       Double Precision, Intent(out) :: twzStart
       Double Precision, Intent(out) :: twzEnd
       Double Precision, Intent(out) :: ttyStart
       Double Precision, Intent(out) :: ttyEnd
       Double Precision, Intent(out) :: ttzStart
       Double Precision, Intent(out) :: ttzEnd
     end Subroutine read_T_region
  
     Subroutine read_region(whichR,rxS,rxE,ryS,ryE,rzS,rzE,denr,clr,csr)
       Integer, Intent(in) :: whichR
       Double Precision, Intent(out) :: rxS
       Double Precision, Intent(out) :: rxE
       Double Precision, Intent(out) :: ryS
       Double Precision, Intent(out) :: ryE
       Double Precision, Intent(out) :: rzS
       Double Precision, Intent(out) :: rzE
       Double Precision, Intent(out) :: denr
       Double Precision, Intent(out) :: clr
       Double Precision, Intent(out) :: csr
     end Subroutine read_region

     Subroutine read_transducer(tposx,tposy,tposz,tthickness,tfreq,   &
        &   pcycles, catchtposx, catchtposy,catchtposz,     &
        &   catchtthickness,tmode) 
       Double Precision , Intent(out) :: tposx
       Double Precision , Intent(out) :: tposy
       Double Precision , Intent(out) :: tposz
       Double Precision , Intent(out) :: tthickness
       Double Precision , Intent(out) :: tfreq
       Double Precision , Intent(out) :: pcycles
       Double Precision , Intent(out) :: catchtposx
       Double Precision , Intent(out) :: catchtposy
       Double Precision , Intent(out) :: catchtposz
      Double Precision , Intent(out) :: catchtthickness
       Integer, Intent(out) :: tmode
     end Subroutine read_transducer

     Subroutine read_surface(rxStart,rxEnd,ryStart,ryEnd,Pid,S)
       Integer, Intent(in) :: rxStart,rxEnd,ryStart,ryEnd
       Integer, Intent(in) :: Pid
       Double Precision, dimension(rxStart:rxEnd,ryStart:ryEnd),    &
            &   Intent(out) :: S(:,:)
     end Subroutine read_surface
      
     Subroutine plot_pulse(dt,dfl,df,tkwin_vect)
       Integer , Intent(in) :: dfl
       Double Precision, Intent(in) :: dt       
       Double Precision, dimension(1:dfl), Intent(in) :: tkwin_vect      
       Double Precision, dimension(1:dfl), Intent(in) :: df
     end Subroutine plot_pulse
     
     Subroutine MPE_Decomp1d(n,Number_procs,Processor_ID,start,last,lsize)
       Integer, Intent(in) :: n
       Integer, Intent(in) :: Number_procs
       Integer, Intent(in) :: Processor_ID
       Integer, Intent(out) :: start
       Integer, Intent(out) :: last
       Integer, Intent(out) :: lsize  
     end Subroutine MPE_Decomp1d

     ! Map variable to the individual CPU functions
     Subroutine map_T_region(starty,endy,ly,tryS,tryE,ltryS,ltryE) !,Sr,Er)
      Integer, intent(in) :: tryS
      Integer, intent(in) :: tryE
      Integer, intent(in) :: starty
      Integer, intent(in) :: endy
      Integer, intent(in) :: ly
      Integer, intent(out) :: ltryS
      Integer, intent(out) :: ltryE
!      Integer, intent(out) :: Sr,Er
    end Subroutine map_T_region

     Subroutine  map_region(startx,starty,endx,endy,lx,ly,   &
          & rxStart,rxEnd,ryStart,ryEnd,lxStart,lxEnd,lyStart,lyEnd,flaw)  
       Integer, intent(in) :: startx,starty
       Integer, intent(in) :: endx,endy,ly,lx
       Integer, intent(in) :: rxStart,rxEnd,ryStart,ryEnd
       Integer, intent(out) :: lxStart,lxEnd,lyStart,lyEnd,flaw
     end Subroutine map_region

     Subroutine change_flaw(startx,starty,endx,endy,rxStart,ryStart,  &
          &   rxEnd,ryEnd,fxlb,fxub,fylb,fyub,flaw,Pid)
       Integer, intent(in) :: startx,starty,endx,endy
       Integer, intent(in) :: rxStart,rxEnd,ryStart,ryEnd
       Integer, intent(out) :: fxlb,fxub,fylb,fyub
       Integer, intent(inout) :: flaw,Pid
     end Subroutine change_flaw

     Subroutine  map_transducer(startx,starty,endx,endy,numz,ly,lx, &
          & tposx1,tposx2,tposy1,tposy2,a,b,c,d)     
       Integer, intent(in)  :: startx,starty
       Integer, intent(in)  :: endx,endy,numz,ly,lx
       Integer, intent(in)  :: tposx1,tposx2,tposy1,tposy2
       Integer, intent(out) :: a,b,c,d
     end Subroutine map_transducer

     Subroutine find_bounds(lxS,lxE,lyS,lyE,rSz,rEz,starty,dens,bounds,Pid)
       Integer, Intent(in) :: lxS,lxE,lyS,lyE,rSz,rEz,starty,Pid
       Double Precision, Dimension(lxS:lxE+2,lyS:lyE+2,rSz:rEz),    &
            &   Intent(inout) :: dens(:,:,:)
       Integer, Dimension(lxS:lxE,lyS:lyE,rSz:rEz), Intent(out) ::  &
            &   bounds(:,:,:)
     end Subroutine find_bounds

!!!!
     ! Update Velocity functions
     Subroutine update_flaw_vel(flx,fly,numz,dtodsp,lxstart,lystart,    &
          &     lzstart,lxend,lyend,lzend,  &
          &   frSz,flxStart,flxEnd,flyStart,flyEnd,flaw,xlb,    &
          &   xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,  &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Double Precision, intent(in) :: dtodsp
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
     end Subroutine update_flaw_vel

     Subroutine section_velocity(lsecxStart,lsecyStart,lseczStart,  &
          &    lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz,dtodsp,  &
          &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
       Integer, intent(in) :: lsecxStart,lsecyStart,lseczStart,     &
            &   lsecxEnd,lsecyEnd,lseczEnd
       Integer, intent(in) :: lx,ly,numz,xlb,xub,ylb,yub,zlb,zub
       Double Precision, intent(in) :: dtodsp
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &   :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            & :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     end Subroutine section_velocity

     Subroutine T_velocity(lx,ly,numx,numz,startx,starty,endy,dtodsp, &
          &     ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,   &
          &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
          &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
          &     rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
          &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
          &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)
       Integer, intent(in) :: lx,ly,numx,numz,startx,starty,endy
       Integer, intent(in) :: ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE
       Integer, intent(in) :: ltwyS,ltwyE,twyS,twyE,twzS,twzE
       Integer, intent(in) :: lttyS,lttyE,ttyS,ttyE,ttzS,ttzE
       Integer, intent(in) :: rSz,flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Integer, intent(in) :: numS,a,b,c
       Double Precision, intent(in) :: dtodsp
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine T_velocity

     Subroutine update_mass_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
          &   lxstart,lystart,lzstart,lxend,lyend,lzend,  &
          &   frSz,frEz,flxStart,flxEnd,flyStart,flyEnd,flaw,    &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
          &   dexStart,dexEnd,deyStart,deyEnd,dxlb,dxub,dylb,dyub,delam,  &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz,frEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Integer, intent(in) :: dexStart,dexEnd,deyStart,deyEnd
       Integer, intent(in) :: dxlb,dxub,dylb,dyub,delam
       Double Precision, intent(in) :: dtodsp,dtodsp2,dden,rden
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
     end Subroutine update_mass_vel

     Subroutine interface_velocity(lsecxStart,lsecyStart,  &
          &    lsecxEnd,lsecyEnd,zlevel,lx,ly,numz,dtodsp,dden,rden,  &
          &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
       Integer, intent(in) :: lsecxStart,lsecyStart,lsecxEnd,lsecyEnd
       Integer, intent(in) :: lx,ly,numz,xlb,xub,ylb,yub,zlb,zub,zlevel
       Double Precision, intent(in) :: dtodsp,dden,rden
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &   :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &   :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     end Subroutine interface_velocity

     Subroutine update_delam_vel(flx,fly,numz,dtodsp,dtodsp2,dden,rden,  &
          &   lxstart,lystart,lzstart,lxend,lyend,lzend,  &
          &   frSz,frEz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz,frEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
       Integer, intent(in) :: fxlb,fxub,fylb,fyub
       Double Precision, intent(in) :: dtodsp,dtodsp2,dden,rden
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
     end Subroutine update_delam_vel

     Subroutine update_sflaw_vel(flx,fly,numz,dtodsp,lxstart,lystart,  &
          &   lzstart,lxend,lyend,lzend,frSz,flxStart,flxEnd,flyStart, &
          &   flyEnd,flaw,xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,  &
          &  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
       Integer, intent(in) :: flx,fly,numz,frSz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
       Integer, intent(in) :: fxlb,fxub,fylb,fyub
       Integer, intent(in) :: a,b,c,starty
       Double Precision, intent(in) :: dtodsp
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            & :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine update_sflaw_vel

     Subroutine update_surf_vel(flx,fly,numz,dtodsp,flxStart,flxEnd,  &
          &   flyStart,flyEnd,rSz,rEz,  &
          &  fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c,starty)
       Integer, intent(in) :: flx,fly,numz,rSz,rEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd
       Integer, intent(in) :: a,b,c,starty
       Double Precision, intent(in) :: dtodsp
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            & :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine update_surf_vel

!!!!!
     ! Update Stress functions
     Subroutine update_flaw_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
          &   frSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
          &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,  &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Double Precision, intent(in) :: dtods,l2m,lambda,mu
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
     end Subroutine update_flaw_stress

     Subroutine section_stress(lsecxStart,lsecyStart,lseczStart,  &
          &    lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz,  &
          &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,zub, &
          &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
       Integer, intent(in) :: lsecxStart,lsecyStart,lseczStart
       Integer, intent(in) :: lsecxEnd,lsecyEnd,lseczEnd,lx,ly,numz
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
       Double Precision, intent(in) :: dtods,l2m,lambda,mu
       Double Precision, Dimension(1:lx,1:ly,1:numz),Intent(in) :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     end Subroutine section_stress

     Subroutine T_stress(lx,ly,numx,numz,startx,starty,endy,dtods, &
          &     l2m,lambda,mu,ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,  &
          &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
          &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
          &     rSz,flxStart,flxEnd,flyStart,flyEnd,flaw,  &
          &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
          &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)
       Integer, intent(in) :: lx,ly,numx,numz,startx,starty,endy
       Integer, intent(in) :: ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE
       Integer, intent(in) :: ltwyS,ltwyE,twyS,twyE,twzS,twzE
       Integer, intent(in) :: lttyS,lttyE,ttyS,ttyE,ttzS,ttzE
       Integer, intent(in) :: rSz,flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Integer, intent(in) :: numS,a,b,c
       Double Precision, intent(in) :: dtods,l2m,lambda,mu
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine T_stress

     Subroutine update_mass_stress(flx,fly,numz,dtods,l2m,l2m2,lambda,  &
          &   lambda2,mu,mu2,frSz,frEz,flxStart,flxEnd,flyStart,flyEnd, &
          &   flaw,lxstart,lystart,lzstart,lxend,lyend,lzend,   &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
          &   dexStart,dexEnd,deyStart,deyEnd,dxlb,dxub,dylb,dyub,delam,  &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz,frEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Integer, intent(in) :: dexStart,dexEnd,deyStart,deyEnd
       Integer, intent(in) :: dxlb,dxub,dylb,dyub,delam
       Double Precision, intent(in) :: dtods,l2m,lambda,mu,l2m2,lambda2,mu2
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
     end Subroutine update_mass_stress

     Subroutine interface_stress(lsecxStart,lsecyStart,  &
          &    lsecxEnd,lsecyEnd,zlevel,lx,ly,numz,  &
          &    dtods,l2m,lambda,mu,mu2,xlb,xub,ylb,yub,zlb,zub, &
          &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
       Integer, intent(in) :: lsecxStart,lsecyStart
       Integer, intent(in) :: lsecxEnd,lsecyEnd,zlevel,lx,ly,numz
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub
       Double Precision, intent(in) :: dtods,l2m,lambda,mu,mu2
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(in)        &
            &    :: vx,vy,vz
       Double Precision, Dimension(1:lx,1:ly,1:numz), Intent(inout)     &
            &    :: Txx,Tyy,Tzz,Txy,Txz,Tyz
     end Subroutine interface_stress

     Subroutine update_delam_stress(flx,fly,numz,dtods,l2m,l2m2,lambda,  &
          &   lambda2,mu,mu2,frSz,frEz,flxStart,flxEnd,flyStart,flyEnd, &
          &   flaw,lxstart,lystart,lzstart,lxend,lyend,lzend,   &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz)
       Integer, intent(in) :: flx,fly,numz,frSz,frEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Double Precision, intent(in) :: dtods,l2m,lambda,mu,l2m2,lambda2,mu2
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
     end Subroutine update_delam_stress

     Subroutine update_sflaw_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
          &   frSz,flxStart,flxEnd,flyStart,flyEnd,flaw, &
          &   lxstart,lystart,lzstart,lxend,lyend,lzend,   &
          &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
       Integer, intent(in) :: flx,fly,numz,frSz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd,flaw
       Integer, intent(in) :: lxstart,lystart,lzstart,lxend,lyend,lzend
       Integer, intent(in) :: xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub
       Integer, intent(in) :: a,b,c
       Double Precision, intent(in) :: dtods,l2m,lambda,mu
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine update_sflaw_stress

     Subroutine update_surf_stress(flx,fly,numz,dtods,l2m,lambda,mu,  &
          &   flxStart,flxEnd,flyStart,flyEnd,rSz,rEz,  &
          &   fvx,fvy,fvz,fTxx,fTyy,fTzz,fTxy,fTxz,fTyz,bounds,a,b,c)
       Integer, intent(in) :: flx,fly,numz,rSz,rEz
       Integer, intent(in) :: flxStart,flxEnd,flyStart,flyEnd
       Integer, intent(in) :: a,b,c
       Double Precision, intent(in) :: dtods,l2m,lambda,mu
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fvx,fvy,fvz
       Double Precision, Dimension(1:flx,1:fly,1:numz), Intent(inout)   &
            &    :: fTxx,fTyy,fTzz,fTxy,fTxz,fTyz
       Integer, Dimension(1:a,1:b,1:c), Intent(in) :: bounds
     end Subroutine update_surf_stress

  end Interface

  
  Double Precision, Parameter :: pi = 3.1415926535
!
! define MPI variables
!  
  Integer :: Pid, N_proc, ierr
!  Integer :: TstartA,nlocal,deficit
  Integer :: MPI_COMM_2D
!  Integer :: status
  Integer :: stride,stride1
  Integer :: period(2),Pdim_size(2),coords(2)
  Integer :: nbrleft,nbrright,nbrtop,nbrbottom
  Integer :: Pdim0, Pdim1


! Model Params 
  ! thickness of absorbing boundary layer (in simulation units ds)
  Double Precision :: abc 
  ! tells if the model is a simple plate or not
  Integer :: plate
  ! number of material regions to add
  Integer :: numR
  ! delamination
  Integer :: numD
  ! 1 if adding surface 0 if not
  Integer :: numS
  ! update the plot every <plotevery> timesteps
  Integer :: plotevery  
  ! maximum frequency in the simulation (used to determine time step)
  Double Precision :: fmax
  ! simulation runtime (seconds)
  Double Precision :: SimulationTime
  ! alpha - variable for tukeywin(Tapered-Cosine)
  Double Precision :: alpha

! Space Properties  
  ! cmax - Fastest Wave Speed (m/s)
  Double Precision :: cmax
  ! cmin - Slowest Wave Speed (m/s)
  Double Precision :: cmin
  ! spacethickness  - space thickness (meters) (z-dir)
  Double Precision :: spacethickness
  ! spacelength - space length (meters) (x-dir)
  Double Precision :: spacelength
  ! spacewidth - space width (meters) (y-dir)
  Double Precision :: spacewidth
  ! minz - minumin z height in space
  Double Precision :: minz
  
! Region Properties
!! Default
  ! dden - Default density
  Double Precision :: dden
  ! dcl - Default Longitudinal Wave Speed
  Double Precision :: dcl
  ! dcs - Default Shear Wave Speed
  Double Precision :: dcs

! T Region Info
  ! tfyStart - T flange start (meters) (y-dir)
  ! tfyEnd - T flange end (meters) (y-dir)
  Double Precision :: tfyStart, tfyEnd
  ! tfzStart - T flange start (meters) (z-dir)
  ! tfzEnd - T flange end (meters) (z-dir)
  Double Precision :: tfzStart, tfzEnd
  ! twyStart - T web start (meters) (y-dir)
  ! twyEnd - T web end (meters) (y-dir)
  Double Precision :: twyStart, twyEnd
  ! twzStart - T web start (meters) (z-dir)
  ! twzEnd - T web end (meters) (z-dir)
  Double Precision :: twzStart, twzEnd
  ! ttyStart - T top start (meters) (y-dir)
  ! ttyEnd - T top end (meters) (y-dir)
  Double Precision :: ttyStart, ttyEnd
  ! ttzStart - T top start (meters) (z-dir)
  ! ttzEnd - T top end (meters) (z-dir)
  Double Precision :: ttzStart, ttzEnd
  ! global variables
  Integer :: tfyS,tfyE,tfzS,tfzE
  Integer :: twyS,twyE,twzS,twzE
  Integer :: ttyS,ttyE,ttzS,ttzE
  ! local variables
  Integer :: ltfyS,ltfyE
  Integer :: ltwyS,ltwyE
  Integer :: lttyS,lttyE

!! Region Info
  ! rxStart - Region start (meters) (x-dir)
  Integer, Allocatable :: rxStart(:)
  ! rxEnd - Refion end (meters) (x-dir)
  Integer, Allocatable :: rxEnd(:)
  ! ryStart - Region start (meters) (y-dir)
  Integer, Allocatable :: ryStart(:)
  ! ryEnd - Refion end (meters) (y-dir)
  Integer, Allocatable :: ryEnd(:)
  ! rzStart - Region start (meters) (z-dir)
  Integer, Allocatable :: rzStart(:)
  ! rzEnd - Refion end (meters) (z-dir)
  Integer, Allocatable :: rzEnd(:)
  ! rden - Density of Region (kg/m^3)
  Double Precision, Allocatable :: rden(:)
  ! rcl - Region Longitudinal Wave Speed
  Double Precision, Allocatable :: rcl(:)
  ! rcs - Region Longitudinal Wave Speed
  Double Precision, Allocatable :: rcs(:)
  ! global region positions
  Integer :: rSx,rEx,rSy,rEy,rSz,rEz
  Integer :: fxlb,fxub,fylb,fyub
  Integer :: dxlb,dxub,dylb,dyub
  ! local region positions
  Integer :: lrxStart,lrxEnd,lryStart,lryEnd,flaw
  Integer :: dexStart,dexEnd,deyStart,deyEnd,delam
  ! Dummy integers
  Double Precision :: rxS,rxE,ryS,ryE,rzS,rzE,denr,clr,csr

  !dens - holds the surface information
  Double Precision, Allocatable :: dens(:,:,:)
  Integer, Allocatable :: bounds(:,:,:)
  Double Precision, Allocatable :: surface(:,:)
  Integer :: a,b,c


  ! Transducer Properties
  ! tposx - transducer position x (meters) - 
  !         from x=0 side of the plate length
  Double Precision :: tposx
  ! tposy - transducer position y (meters) - 
  !         from y=0 of the plate width
  Double Precision :: tposy
  ! tposz - transducer position z (meters) - 
  !         from z=0 of the space height
  Double Precision :: tposz
  ! tthickness - transducer diameter  (meters) - this is 2D
  Double Precision :: tthickness
  ! tfreq - frequency of transducer (Hertz) - must be less than fmax below
  Double Precision :: tfreq
  ! pcycles - number of cycles in the pulse
  Double Precision :: pcycles        
  ! tpulselength - transducer pulse length (seconds)
  Double Precision :: tpulselength 
  ! catchtposx - catch transducer position x (meters) - 
  !         from x=0 side of the plate length
  Double Precision :: catchtposx
  ! catchtposy - catch transducer position y (meters) - 
  !         from y=0 side of the plate width
  Double Precision :: catchtposy
  ! catchtposz - catch transducer position z (meters) - 
  !         from z=0 side of the space height
  Double Precision :: catchtposz
  ! catchtthickness - catch transducer diameter  (meters) - this is 2D
  Double Precision :: catchtthickness
  ! transducer mode - 1=vx 2=vy 3=vz
  Integer :: tmode
  ! global positions of the transducer
  Integer :: tposx1,tposx2     !  tposx1 - left edge
  Integer :: ctposx1,ctposx2   !  tposx2 - right edge
  Integer :: tposy1,tposy2     !  tposy1 - front edge
  Integer :: ctposy1,ctposy2   !  tposy2 - back edge
  ! local positions of the transducer
  Integer :: ltposx1,ltposx2
  Integer :: lctposx1,lctposx2
  Integer :: ltposy1,ltposy2
  Integer :: lctposy1,lctposy2

  
! Txx - Normal Stress Values in x direction
  Double Precision, Allocatable :: Txx(:,:,:)
! Tyy - Normal Stress Values in y direction
  Double Precision, Allocatable :: Tyy(:,:,:)
! Tzz - Normal Stress Values in z direction
  Double Precision, Allocatable :: Tzz(:,:,:)
! Txy - Sheer Stress Value 
  Double Precision, Allocatable :: Txy(:,:,:)
! Txz - Sheer Stress Value 
  Double Precision, Allocatable :: Txz(:,:,:)
! Tyz - Sheer Stress Value 
  Double Precision, Allocatable :: Tyz(:,:,:)
! vx -  Velocity values in the x direction
  Double Precision, Allocatable :: vx(:,:,:)
! vy - Velocity values in the y direction 
  Double Precision, Allocatable :: vy(:,:,:)
! vz - Velocity values in the z direction 
  Double Precision, Allocatable :: vz(:,:,:)
! disp - displacement 
  Double Precision, Allocatable :: disp(:,:,:)
! catchAline - Array to Store the Catch Aline 
  Double Precision, Allocatable :: catchAline(:)

! tkwin_vect -  tukeywin transform
  Double Precision, Allocatable :: tkwin_vect(:)
! df -  pulse form
  Double Precision, Allocatable :: df(:)

  Double Precision :: ds,dt,lambda,mu,lambda2,mu2
  Double Precision :: mean,Partial_mean
  Double Precision :: dtods,dtodsp,l2m,dtodsp2,l2m2
  Double Precision :: t1,t2
   
  Integer :: numt, t
  Integer :: numx,numy,numz,minnumz,tzs,catchtzs
  Integer :: dfl
  
  Integer :: startx,endx,starty,endy !,startz,endz 
  Integer :: s1_index,s2_index
  Integer :: i,j,k,iz,r !,ix,iy,tix,tiy,tiz,r,xe,xs,xs1,xs2
  Integer :: lx,ly,lz
  Integer :: xlb,xub,ylb,yub,zlb,zub

! Silo variables
  Integer :: silo_on
  Integer :: siloSpace,silo_lx,silo_ly
  Integer :: dbfile, lpname, ldiscp, ID, err, lfname, oldlen
  Character(40) :: pname,discp,fname
  Character(5) :: tt,pp,ii
  Double Precision, Allocatable :: xx(:),yy(:),zz(:)
  Integer, Allocatable :: dims(:),lmeshnames(:),meshtypes(:)
  Character(40), Allocatable :: meshnames(:)
  Integer, Allocatable :: ldispnames(:),vartypes(:)
  Character(40), Allocatable :: dispnames(:)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Pid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,N_proc,ierr) 
!  write(*,*) 'flag A, MPI_COMM_WORLD:',MPI_COMM_WORLD



!!!! Start defining the space by reading input files

  if (Pid .EQ. 0) then
     call read_model(alpha,plate,numR,numD,numS,plotevery,silo_on,  &
        &   siloSpace,fmax,SimulationTime,Pdim0,Pdim1)
     call read_transducer(tposx,tposy,tposz,tthickness,tfreq, pcycles,   &
        &   catchtposx, catchtposy,catchtposz, catchtthickness,tmode)
     call read_dmaterials(cmax,cmin,spacethickness,spacelength, &
        &   spacewidth,dden,dcl,dcs,minz)
     write(*,*)
  end if

  call MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(plate,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(numR,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(numD,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(plotevery,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(silo_on,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(siloSpace,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(numS,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(fmax,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(SimulationTime,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Pdim0,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Pdim1,1,MPI_INTEGER,0, &
       &  MPI_COMM_WORLD,ierr)

  if (numR > 1) then
     Allocate(rxStart(1:numR))
     Allocate(rxEnd(1:numR))
     Allocate(ryStart(1:numR))
     Allocate(ryEnd(1:numR))
     Allocate(rzStart(1:numR))
     Allocate(rzEnd(1:numR))
     Allocate(rden(1:numR))
     Allocate(rcl(1:numR))
     Allocate(rcs(1:numR))
  else
     Allocate(rxStart(1))
     Allocate(rxEnd(1))
     Allocate(ryStart(1))
     Allocate(ryEnd(1))
     Allocate(rzStart(1))
     Allocate(rzEnd(1))
     Allocate(rden(1))     
     Allocate(rcl(1))
     Allocate(rcs(1))     
  end if

  if (Pid .EQ. 0) then
     ! Compute number of nodes for the minumin thickness 
     minnumz = ceiling((10*fmax*minz)/cmin) + 1
     ! compute spatial step
     ds = minz/(minnumz-1)
     write(*,*)' minnumz: ',minnumz
     write(*,*)' ds: ',ds
     write(*,*)

     if (plate == 2) then
        call read_T_region(tfyStart,tfyEnd,tfzStart,tfzEnd,twyStart,    &
            &  twyEnd,twzStart,twzEnd,ttyStart,ttyEnd,ttzStart,ttzEnd)
        write(*,*)
        tfyS = ceiling(tfyStart/ds)+1
        tfyE = ceiling(tfyEnd/ds)+1
        tfzS = ceiling(tfzStart/ds)+1
        tfzE = ceiling(tfzEnd/ds)+1
        twyS = ceiling(twyStart/ds)+1
        twyE = ceiling(twyEnd/ds)+1
        twzS = ceiling(twzStart/ds)+1
        twzE = ceiling(twzEnd/ds)+1
        ttyS = ceiling(ttyStart/ds)+1
        ttyE = ceiling(ttyEnd/ds)+1
        ttzS = ceiling(ttzStart/ds)+1
        ttzE = ceiling(ttzEnd/ds)+1

        write(*,*) ' tfyS = ',tfyS
        write(*,*) ' tfyE = ',tfyE
        write(*,*) ' tfzS = ',tfzS
        write(*,*) ' tfzE = ',tfzE
        write(*,*)
        write(*,*) ' twyS = ',twyS
        write(*,*) ' twyE = ',twyE
        write(*,*) ' twzS = ',twzS
        write(*,*) ' twzE = ',twzE
        write(*,*)
        write(*,*) ' ttyS = ',ttyS
        write(*,*) ' ttyE = ',ttyE
        write(*,*) ' ttzS = ',ttzS
        write(*,*) ' ttzE = ',ttzE
        write(*,*)
     else
        tfyS = 0
        tfyE = 0
        tfzS = 0
        tfzE = 0
        twyS = 0
        twyE = 0
        twzS = 0
        twzE = 0
        ttyS = 0
        ttyE = 0
        ttzS = 0
        ttzE = 0
     end if
     
     if (numR > 0) then
        do i = 1,numR
           call read_region(i,rxS,rxE,ryS,ryE,rzS,rzE,denr,clr,csr)
           rxStart(i) = ceiling(rxS/ds)+1
           rxEnd(i) = ceiling(rxE/ds)+1
           ryStart(i) = ceiling(ryS/ds)+1
           ryEnd(i) = ceiling(ryE/ds)+1
           rzStart(i) = ceiling(rzS/ds)+1
           rzEnd(i) = ceiling(rzE/ds)+1
           rden(i) = denr
           rcl(i) = clr
           rcs(i) = csr
        end do
        write(*,*)' rxStart: ',rxStart
        write(*,*)' rxEnd: ',rxEnd
        write(*,*)' ryStart: ',ryStart
        write(*,*)' ryEnd: ',ryEnd
        write(*,*)' rzStart: ',rzStart
        write(*,*)' rzEnd: ',rzEnd
        write(*,*)' rden: ',rden
        write(*,*)' rcl: ',rcl
        write(*,*)' rcs: ',rcs
        write(*,*)
    else
        rxStart = 0
        rxEnd = 0
        ryStart = 0
        ryEnd = 0
        rzStart = 0
        rzEnd = 0
        rden = 0.0
        rcl = 0.0
        rcs = 0.0
     end if
  end if
    
!
! Broad cast input parameters values to all processors
!
  call MPI_Bcast(cmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(spacethickness,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(spacelength,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(spacewidth,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dden,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dcl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dcs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(minz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(tfyS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tfyE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tfzS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tfzE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(twyS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(twyE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(twzS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(twzE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ttyS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ttyE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ttzS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ttzE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  if (numR > 0) then
     call MPI_Bcast(rxStart,numR,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rxEnd,numR,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ryStart,numR,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ryEnd,numR,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rzStart,numR,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rzEnd,numR,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rden,numR,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rcl,numR,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rcs,numR,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     call MPI_Bcast(rxStart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rxEnd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ryStart,1,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ryEnd,1,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rzStart,1,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rzEnd,1,MPI_INTEGER,0, &
          &  MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rden,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rcl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(rcs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if

  call MPI_Bcast(tposx,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tposy,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tposz,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tthickness,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tfreq,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(pcycles,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(catchtposx,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(catchtposy,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(catchtposz,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(catchtthickness,1,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tmode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  Pdim_size(1) = Pdim0
  Pdim_size(2) = Pdim1

  period(1) = 0
  period(2) = 0

  call MPI_Dims_create(N_proc,2,Pdim_size,ierr)   
  call MPI_Cart_create(MPI_COMM_WORLD,2,Pdim_size,period,1,     &
        &   MPI_COMM_2D,ierr)
  call MPI_Comm_rank(MPI_COMM_2D,Pid,ierr);

!  neighbors to the left and the right of Pid
  call MPI_Cart_shift(MPI_COMM_2D,0,1,nbrleft,nbrright,ierr)
!  neighbors above and below   
  call MPI_Cart_shift(MPI_COMM_2D,1,1,nbrtop,nbrbottom,ierr) 
!  get mpi cartisian coordinates
  call MPI_Cart_get(MPI_COMM_2D,2,Pdim_size,period,coords,ierr)  

!
!Calc intial parameter values
!
! pulse length (seconds)
  tpulselength = pcycles*(1/tfreq)

! Compute number of nodes for the minumin thickness 
  minnumz = ceiling((10*fmax*minz)/cmin) + 1

! compute spatial step
  ds = minz/(minnumz-1)
  
! compute number of nodes in z-direction
  numz = ceiling(spacethickness/ds)+1

! compute number of nodes in x-direction
  numx = ceiling(spacelength/ds)+1

! compute number of nodes in y-direction
  numy = ceiling(spacewidth/ds)+1
   
! compute time step
  dt = 1/(2*cmax*sqrt(3/((ds)**2)))

! number of time steps
  numt=ceiling(SimulationTime/dt)
  
! compute default lame constant - lambda  lm=p(cp^2-2*cs^2) 
  lambda  = dden*(dcl**2-2*dcs**2)

! compute default lame constant - mu      mu=p*cs^2  
  mu = dden*(dcs*dcs)  
  
! Simulation Constant
  dtodsp = dt/(dden*ds)
  dtods = dt/ds
  l2m = lambda+2*mu

  if ((numR > 0) .and. (rden(1) > 0.0)) then
     lambda2  = rden(1)*(rcl(1)**2-2*rcs(1)**2)
     mu2 = rden(1)*(rcs(1)*rcs(1))  
     dtodsp2 = dt/(rden(1)*ds)
     l2m2 = lambda2+2*mu2
  else
     lambda2  = lambda
     mu2 = mu 
     dtodsp2 = dtodsp
     l2m2 = l2m
  end if

  
!  write(*,*) 'Pid:',Pid,' dtodsp: ',dtodsp,dtodsp2
!  write(*,*) 'Pid:',Pid,' lambda: ',lambda,lambda2
!  write(*,*) 'Pid:',Pid,' mu: ',mu,mu2

! drive function length in dt units
  dfl  = ceiling(tpulselength/dt)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (Pid == 0) then 
     write(*,*)' delta T dt = ',dt
     write(*,*)' number of time steps = ',numt
     write(*,*)
     write(*,*)' step size ds = ',ds 
     write(*,*)' number of grid points in x dir = ',numx
     write(*,*)' number of grid points in y dir = ',numy
     write(*,*)' number of grid points in z dir = ',numz
     write(*,*)' number of grid points in min z = ',minnumz
     write(*,*)
     write(*,*)' lambda default = ',lambda, lambda2
     write(*,*)' mu default = ',mu,mu2
     write(*,*)' dtods = ', dtods
     write(*,*)
     write(*,*)' number of dt steps the drive function is on ',dfl
     write(*,*)
   end if
    
! call domain decomposition function here
! determines the size for local num1 and local num2 values then 
! allocate arrays  including ghost boundaries
! print the results to screen for testing then  
! 
 
  call MPE_Decomp1d(numx,Pdim_size(1),coords(1),startx,endx,lx)
  call MPE_Decomp1d(numy,Pdim_size(2),coords(2),starty,endy,ly)

!  startz = 1
!  endz = numz


! startx,starty,endx,endy are the global index numbers
! lx, ly are the size of the local index(1:lx,1:ly)

  call MPI_TYPE_VECTOR(ly,1,lx,MPI_DOUBLE_PRECISION,stride,ierr)
  call MPI_TYPE_VECTOR(lx,1,1,MPI_DOUBLE_PRECISION,stride1,ierr)

  call MPI_TYPE_COMMIT(stride,ierr)
  call MPI_TYPE_COMMIT(stride1,ierr)
     

! Silo Spacing
!  silo_lx = lx !/siloSpace
!  silo_ly = ly !/siloSpace

!  if (Pid == 0) then
!     write(*,*)' silo_lx = ',silo_lx
!     write(*,*)' silo_ly = ',silo_ly
!     write(*,*)' silo_lz = ',numz
!     write(*,*)
!  end if

!
! allocate Txx Tyy Tzz Txy Txz Tyz vx vy vz and catchaline arrays
  Allocate(Txx(1:lx,1:ly,1:numz))
  Allocate(Tyy(1:lx,1:ly,1:numz))
  Allocate(Tzz(1:lx,1:ly,1:numz))
  Allocate(Txy(1:lx,1:ly,1:numz))
  Allocate(Txz(1:lx,1:ly,1:numz))
  Allocate(Tyz(1:lx,1:ly,1:numz))
   
  Allocate(vx(1:lx,1:ly,1:numz)) 
  Allocate(vy(1:lx,1:ly,1:numz))
  Allocate(vz(1:lx,1:ly,1:numz))
  Allocate(disp(1:lx,1:ly,1:numz))

  Allocate(catchAline(1:numt))  
  
! Allocate silo write arrays
  Allocate(xx(1:lx))
  Allocate(yy(1:ly))
  Allocate(zz(1:numz))
  Allocate(dims(1:3))
  Allocate(meshnames(1:N_proc))
  Allocate(lmeshnames(1:N_proc))
  Allocate(meshtypes(1:N_proc))
  Allocate(dispnames(1:N_proc))
  Allocate(ldispnames(1:N_proc))
  Allocate(vartypes(1:N_proc))
  
!
! make the pulse with a tukey window
!
  Allocate(tkwin_vect(1:dfl))
  Allocate(df(1:dfl))
    
!
!  Initialize Drive Function - df
  if (Pid == 0) then
     tkwin_vect = 0.0
     df = 0.0
       
     call tukeywin(tkwin_vect,dfl,alpha)  
       
     do i = 1,dfl
        df(i) = sin((i-1)*dt*tfreq*2*pi) * tkwin_vect(i)
     end do
  end if

  call MPI_Bcast(df,dfl,MPI_DOUBLE_PRECISION,0, &
       &  MPI_COMM_WORLD,ierr)
    
!
! Build density arrray
!-------------------------------------------------------------------------

  if (plate == 2) then
     call map_T_region(starty,endy,ly,tfyS,tfyE,ltfyS,ltfyE)
     call map_T_region(starty,endy,ly,twyS,twyE,ltwyS,ltwyE)
     call map_T_region(starty,endy,ly,ttyS,ttyE,lttyS,lttyE)
  else
      ltfyS = 0
      ltfyE = 0
      ltwyS = 0
      ltwyE = 0
      lttyS = 0
      lttyE = 0
  end if
     
  if (numR > 0) then
! Change starts and ends into ds units
     do r = 1,numR
        rSx = rxStart(r)
        rEx = rxEnd(r)
        rSy = ryStart(r)
        rEy = ryEnd(r)
        rSz = rzStart(r)
        rEz = rzEnd(r)
        
        if (r == 1) then
           ! Map flaw region dimensions to the local processor
           call map_region(startx,starty,endx,endy,lx,ly,   &
                &   rSx,rEx,rSy,rEy,lrxStart,lrxEnd,lryStart,lryEnd,  &
                &   flaw)
           call change_flaw(startx,starty,endx,endy,rSx,rSy,  &
                &   rEx,rEy,fxlb,fxub,fylb,fyub,flaw,Pid)
           !       write(*,*)'Pid:',Pid,' rSz:',rSz,' rEz:',rEz
           delam = 0
        else if (r == 2) then
           ! Map flaw region dimensions to the local processor
           call map_region(startx,starty,endx,endy,lx,ly,   &
                &   rSx,rEx,rSy,rEy,dexStart,dexEnd,deyStart,deyEnd,  &
                &   delam)
!                  write(*,*)'Pid:',Pid,' rSz:',rSz,' rEz:',rEz
           call change_flaw(startx+lrxStart-1,starty+lryStart-1,  &
                &   startx+lrxEnd-1,starty+lryEnd-1,rSx,rSy,  &
                &   rEx,rEy,dxlb,dxub,dylb,dyub,delam,Pid)

        end if
     end do
  else
     flaw = 0
     delam = 0
  end if
  

  ! Allocating and reading in the corrosion surface.
  if (numR == 1 .and. numS == 1) then
     if (flaw > 0) then
        a = lrxEnd-lrxStart+1
        b = lryEnd-lryStart+1
        c = rEz-rSz+1
        write(*,*)'Pid:',Pid,'global:',lrxStart,lryStart,rSz,'local:',a,b,c

        Allocate(dens(1:a+2,1:b+2,1:c))
        Allocate(surface(1:a+2,1:b+2))
        dens = -666.0
        call read_surface(1,a+2,1,b+2,Pid,surface(:,:))
        surface = surface-rSz
        ! build rest of density array
        write(*,*)'FLAG flaw:',Pid,flaw,a,b,c
        do i = 1,a+2
           do j = 1,b+2
              do k = 1,c
                 if (k <= surface(i,j)) then
                    dens(i,j,k) = dden
                 else
                    dens(i,j,c) = -555.0
                 end if
              end do
           end do
        end do
!        if (Pid == 2) then
!           dens(1:3,:,:) = dden
!           dens(:,b:b+2,:) = dden
!           dens(:,1,:) = -111.0
!        elseif (Pid == 4) then
!           dens(a:a+2,:,:) = dden
!           dens(:,b:b+2,:) = dden
!           dens(:,1,:) = -111.0
!        end if

        write(*,*)'Pid:',Pid,'surf:',shape(surface)
 
        Deallocate(surface)
        Allocate(bounds(1:a,1:b,1:c))

!lrxStart,lrxEnd,lryStart,lryEnd,rSz,rEz
        call find_bounds(1,a,1,b,1,c,starty,dens,bounds,Pid)
!        if (Pid == 2) then
!           write(*,*)
!           write(*,*),bounds
!        end if
!           do j = 40,50
!              do i = 1,c
!                 write(*,*)'Pid:',Pid,j,i,'bounds:',bounds(20,j,i)
!              end do
!           end do
!        end if
     end if
  else
     Allocate(dens(0:0,0:0,0:0))
     Allocate(bounds(0:0,0:0,0:0))
  end if
  

!-------------------------------------------------------------------------
! 
! Build the x's and y's of the transducers
!
!
! left edge and front edge of drive transducer 
!
  tposx1 = ceiling((tposx-(tthickness/2))/ds)
  tposy1 = ceiling((tposy-(tthickness/2))/ds)

!
! right edge and back edge of drive transducer
!
  tposx2 = tposx1 + floor(tthickness/ds)
  tposy2 = tposy1 + floor(tthickness/ds)

!
!  left edge and front edge of catch transducer 
!

  ctposx1 = ceiling((catchtposx - (catchtthickness/2))/ds)
  ctposy1 = ceiling((catchtposy - (catchtthickness/2))/ds)

!
! right edge and back edge of catch transducer
!
  ctposx2 = ctposx1 + floor(catchtthickness/ds)
  ctposy2 = ctposy1 + floor(catchtthickness/ds)
      
  call map_transducer(startx,starty,endx,endy,numz,ly,lx, &
       & tposx1,tposx2,tposy1,tposy2,ltposx1,ltposx2,ltposy1,ltposy2)

  call map_transducer(startx,starty,endx,endy,numz,ly,lx,ctposx1,  &
    &   ctposx2,ctposy1,ctposy2,lctposx1,lctposx2,lctposy1,lctposy2)
 
! find z level of pitch transducers position
! searching botton up  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if (plate == 1) then
!     tzs = numz
!  else if (plate == 2) then
!     tzs = tfzE
!  else if(plate == 3) then
!     tzs = numz
!  else
!      tzs = -5
!  end if

  tzs = ceiling(tposz/ds)
  
  if (Pid == 0) then
     write(*,*)' Pid = ',Pid,' tzs = ',tzs
  end if

! find z level of catch transducer position
! searching botton up  
!  if (plate == 1) then
!     catchtzs = numz
!  else if (plate == 2) then
!     catchtzs = tfzE
!  else if (plate == 3) then
!     catchtzs = rEz
!  else
!      catchtzs = -5
!  end if

  catchtzs = ceiling(catchtposz/ds)

  if (Pid == 0) then
     write(*,*)' Pid = ',Pid,' catchtzs = ',catchtzs
     write(*,*)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  if (Pid == 0) then
!
! plot the pulse signal and tukeywin filter
!
     call plot_pulse(dt,dfl,df,tkwin_vect)
!
! write velocity output header for plotmtv
!
! open output files
     open(112,FILE = "results/aline.mtv",STATUS = "UNKNOWN",    &
            &   ACCESS = "SEQUENTIAL")
     write(112,*)"$ DATA=CURVE2D"
     write(112,*)" "
     write(112,*)"% XMIN=0"
     write(112,*)"% XMAX=",SimulationTime
     write(112,*)"% YMIN=-4E-8"
     write(112,*)"% YMAX=4E-8"
     write(112,*)" "
      
  end if  

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   xlb = 0
   xub = 0
   ylb = 0 
   yub = 0
   zlb = 1
   zub = 1

  s1_index = 1
  s2_index = 1
  
  if (startx == 1) then
     s1_index = 2
     xlb = 1
!     write(*,*)' Pid:',Pid,'xlb=1'
  end if
  if (starty == 1) then
     s2_index = 2
     ylb = 1
!     write(*,*)' Pid:',Pid,'ylb=1'
  end if
  if (endx == numx) then
     xub = 1
!     write(*,*)' Pid:',Pid,'xub=1'
  end if
  if (endy == numy) then
     yub = 1
!     write(*,*)' Pid:',Pid,'yub=1'
  end if
  
! Check to see if transducer contacts sample
if ((tzs < 0) .AND. (ltposx1 > 0)) then
   write(*,*)' Transducer not in contact with sample'
else

! Initialize arrays to zero
   vx = 0.0
   vy = 0.0
   vz = 0.0
   Txx = 0.0
   Tyy = 0.0
   Tzz = 0.0
   Txy = 0.0
   Txz = 0.0
   Tyz = 0.0
   disp = 0.0
     
! Processor time
  t1 = MPI_WTIME()
  
!
! Time Loop
!=========================================================================
  do t= 1,numt
!  do t = 1,2000
     if (Pid == 0) then 
        write(*,*)' timestep t = ',t
     end if

!
! UPDATE VELOCITIES
! 
     if (plate == 1) then
        if (flaw > 0) then
           SELECT CASE (numS)
              CASE (0)  ! no surface
                 call  update_flaw_vel(lx,ly,numz,dtodsp,1,1,1,lx,ly,  &
                      &   numz,rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1) ! surface
!                 write(*,*)
!                 write(*,*)'F1 flag vel'
                 call  update_sflaw_vel(lx,ly,numz,dtodsp,1,1,1,lx,ly,  &
                      &   numz,rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw,  &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub, &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,  &
                      &   a,b,c,starty)
           end SELECT
        else
           call section_velocity(1,1,1,lx,ly,numz,lx,ly,numz,dtodsp,  &
                &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,      &
                &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        end if
     else if (plate == 2) then
        ! T stringer
!        write(*,*)'F1 flag vel'
        call  T_velocity(lx,ly,numx,numz,startx,starty,endy,dtodsp, &
             &     ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,  &
             &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
             &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
             &     rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw,  &
             &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
             &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)
     else if (plate == 3) then 
!        ! Mass loading in place of flaw
        if (flaw > 0) then
!           write(*,*)'Pid: ',Pid,' Flag 1 ',flaw,rEz
           call  update_mass_vel(lx,ly,numz,dtodsp,dtodsp2,dden,rden(1), &
                &   1,1,1,lx,ly,rSz,rSz,rEz,lrxStart,lrxEnd,  &
                &   lryStart,lryEnd,flaw,xlb,xub,ylb,yub,zlb,zub,  &
                &   fxlb,fxub,fylb,fyub,dexStart,dexEnd,deyStart,deyEnd, &
                &   dxlb,dxub,dylb,dyub,delam,  &
                &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
!           write(*,*)'Pid: ',Pid,' Flag 2 ',flaw
        else
           call section_velocity(1,1,1,lx,ly,rSz,lx,ly,numz,dtodsp,  &
                &    xlb,xub,ylb,yub,zlb,zub,vx,vy,vz,      &
                &    Txx,Tyy,Tzz,Txy,Txz,Tyz)
        end if
     else
        if (Pid == 0) then
           write(*,*)' Did nothing to velocities'
        end if
     end if
!------------------------------------------------------------------------

!
! DRIVE FUNCTION
!
 
! this is the top of the plate under the transmitting transducer
     if (t<dfl) then        
        if (ltposx1 > 0) then   
           do i = ltposx1,ltposx2 
              do j = ltposy1,ltposy2
                 if (tmode == 1) then
                    vx(i,j,tzs) =  vx(i,j,tzs) + dtods*(1/dden)*df(t)
                 else if (tmode == 2) then
                    vy(i,j,tzs) =  vy(i,j,tzs) + dtods*(1/dden)*df(t)
                 else
                    vz(i,j,tzs) =  vz(i,j,tzs) + dtods*(1/dden)*df(t)
                 end if
              end do
           end do
        end if
     end if

!
! Catch Aline at catch transducer
!
     Partial_mean = 0.0000
     if (lctposx1 > 0) then
        do i= lctposx1,lctposx2
           do j = lctposy1,lctposy2
              Partial_mean = Partial_mean + vz(i,j,catchtzs)
           end do
        end do
     end if
     
     mean = 0.000000
     
     call MPI_Allreduce(Partial_mean,mean,1,MPI_DOUBLE_PRECISION, &
          &   MPI_SUM,MPI_COMM_WORLD,ierr)
     
     catchAline(t) = mean/((floor(catchtthickness/ds))**2)
     
!
!   write out catchAline to velocity.mtv
!
     if (Pid == 0) then        
        write(112,*) t*dt , catchAline(t)
     end if

!     
! Absorbing Boundary Conditions : should go here
!

!
! MPI Communication Send and Receive VELOCITIES!!!
!
     do iz = 1,numz
        ! Update vx 
        call exchange_vx(vx,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
        ! Update vy
        call exchange_vy(vy,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
        ! Update vz
        call exchange_vz(vz,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
     end do
     
     
!
! UPDATE stresses!!!!!!
!
     if (plate == 1) then
        if (flaw > 0) then
           SELECT CASE (numS)
              CASE (0)
                 call update_flaw_stress(lx,ly,numz,dtods,l2m,lambda,mu,  &
                      &   rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw, &
                      &   1,1,1,lx,ly,numz,xlb,xub,ylb,yub,zlb,zub,   &
                      &   fxlb,fxub,fylb,fyub,  &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
              CASE (1)
                 call update_sflaw_stress(lx,ly,numz,dtods,l2m,lambda,mu, &
                      &   rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw, &
                      &   1,1,1,lx,ly,numz,   &
                      &   xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,   &
                      &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,bounds,a,b,c)
              end SELECT
        else
           call section_stress(1,1,1,lx,ly,numz,lx,ly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,zub, &
                &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        end if
!!!------------------------------------------------------------------------
     else if (plate == 2) then
        call T_stress(lx,ly,numx,numz,startx,starty,endy,dtods,l2m, &
                &     lambda,mu,ltfyS,ltfyE,tfyS,tfyE,tfzS,tfzE,  &
                &     ltwyS,ltwyE,twyS,twyE,twzS,twzE,   &
                &     lttyS,lttyE,ttyS,ttyE,ttzS,ttzE, &
                &      rSz,lrxStart,lrxEnd,lryStart,lryEnd,flaw,  &
                &     xlb,xub,ylb,yub,zlb,zub,fxlb,fxub,fylb,fyub,    &
                &     vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,numS,bounds,a,b,c)
     else if (plate == 3) then 
!        ! Mass loading in place of flaw
        if (flaw > 0) then
           call update_mass_stress(lx,ly,numz,dtods,l2m,l2m2,lambda,  &
                &   lambda2,mu,mu2,rSz,rEz,lrxStart,lrxEnd,lryStart, &
                &   lryEnd,flaw,1,1,1,lx,ly,rSz,xlb,xub,ylb,yub,zlb,zub,  &
                &   fxlb,fxub,fylb,fyub,dexStart,dexEnd,deyStart,deyEnd,  &
                &   dxlb,dxub,dylb,dyub,delam,  &
                &   vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        else
           call section_stress(1,1,1,lx,ly,rSz,lx,ly,numz,  &
                &    dtods,l2m,lambda,mu,xlb,xub,ylb,yub,zlb,zub, &
                &	  vx,vy,vz,Txx,Tyy,Tzz,Txy,Txz,Tyz)
        end if
     else
        if (Pid == 0) then
           write(*,*) 'Did nothing to stresses'
        end if
     end if
!--------------------------------------------------------------------------

     
!
! MPI Communication Send and Receive STRESSES!!!
!
     do iz = 1,numz
!
! Update Txx 
        call exchange_Txx(Txx,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
! Update Tyy 
        call exchange_Tyy(Tyy,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
! Update Tzz
        call exchange_Tzz(Tzz,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
! Update Txy 
        call exchange_Txy(Txy,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)     
! Update Txz
        call exchange_Txz(Txz,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D)
! Update Tyz
        call exchange_Tyz(Tyz,lx,ly,numz,iz,stride,stride1, &
             & nbrleft,nbrright,nbrbottom,nbrtop,MPI_COMM_2D) 
     end do



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!! Calls to write silo files
     if (silo_on == 1) then
        if ((MOD(t,plotevery) == 0 )) then

! Update displacement array
           do i = 1,lx
              do j = 1,ly
                 do k = 1,numz
                    disp(i,j,k) = sqrt(vx(i,j,k)**2 +   &
                         &     vy(i,j,k)**2 + vz(i,j,k)**2)
!              if (disp(i,j,k) >= 1) then
!                 write(*,*)'over 1, x,y,z:',i,j,k
!              end if
                 end do
              end do
           end do
!           write(*,*)' flag1'

           
!! without Subroutine writeSilo
           write(unit=tt,fmt='(I5)') t
           write(unit=pp,fmt='(I5)') Pid
           pname = 'results/velVector_'//TRIM(ADJUSTL(pp))//    &
                &   '_'//TRIM(ADJUSTL(tt))//'.silo'
           lpname = LEN_TRIM(pname)
           discp = "3d velocity mesh"
           ldiscp = LEN_TRIM(discp)
           
           ierr = dbcreate(pname,lpname,DB_NOCLOBBER,DB_LOCAL,      &
                    &   discp,ldiscp,DB_PDB,dbfile)
           if (dbfile == -1) then 
              write(*,*) ' could not create Silo file'
           endif
!           ierr = dbmkdir(dbfile,'mesh_data',9,ID)
!           ierr = dbsetdir(dbfile,'/mesh_data',10)

           do i = 1,lx
              xx(i) = (i+startx-1.0)*ds
           end do
           do j = 1,ly
              yy(j) = (j+starty-1.0)*ds
           end do
           do k = 1,numz
              zz(k) = k*1.0*ds
           end do
           dims(1) = lx
           dims(2) = ly
           dims(3) = numz

           err = dbputqm(dbfile,'rectmesh',8,'x',1,'y',1,'z',1,     &
                &   xx,yy,zz,dims,3,DB_DOUBLE,DB_COLLINEAR,DB_F77NULL,ID)
           err = dbputqv1(dbfile,'disp',4,'rectmesh',8,disp,dims,3,     &
                &   DB_F77NULL,0,DB_DOUBLE,DB_NODECENT,DB_F77NULL,ierr)
           ierr = dbclose(dbfile)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!

           if (Pid == 0) then

!! without Subroutine writeSilo
              write(unit=tt,fmt='(I5)') t
              fname = 'results/totalvelmesh_'//TRIM(ADJUSTL(tt))//'.silo'
              lfname = LEN_TRIM(fname)
              discp = 'Total vector velocities for simulation at t '//  &
                    &   TRIM(ADJUSTL(tt))
              ldiscp = LEN_TRIM(discp)
              write(*,*)' root filename = ',fname
           
              ! Create a new silo file
              ierr = dbcreate(fname,lfname,DB_NOCLOBBER,DB_LOCAL,discp, &
                        &   ldiscp,DB_PDB,dbfile)
              if(dbfile == -1) then
                 write(*,*)' Could not create Silo root file!'
              end if

              do i = 1,N_proc
                 write(unit=ii,fmt='(I5)') i-1
           !  meshnames(i) = ('velVector_'//ii//'_'//tt//'.silo:rectmesh')
                 meshnames(i) = ADJUSTL('velVector_'//TRIM(ADJUSTL(ii))// &
                    &   '_'//TRIM(ADJUSTL(tt))//'.silo:/rectmesh')
                 lmeshnames(i) = LEN_TRIM(meshnames(i))
                 meshtypes(i) = DB_QUAD_RECT
               !  dispnames(i) = ('velVector_'//ii//'_'//tt//'.silo:/disp')
                 dispnames(i) = ADJUSTL('velVector_'//TRIM(ADJUSTL(ii))// 
                    &   '_'//TRIM(ADJUSTL(tt))//'.silo:/disp')
                 ldispnames(i) = LEN_TRIM(dispnames(i))
                 vartypes(i) = DB_QUADVAR
              end do
              
!              do i = 1,N_proc
!                 write(*,*) meshnames(i),lmeshnames(i),dispnames(i),   &
                        &   ldispnames(i)
!              end do

           ! Set the maximum string length to 20 since 
           !    that's how long our strings are
              oldlen = dbget2dstrlen()
              err = dbset2dstrlen(40)

           ! Write the multimesh object
              err = dbputmmesh(dbfile,'rectmesh',8,N_proc,meshnames,    &
                        &   lmeshnames,meshtypes,DB_F77NULL,ierr)
              err = dbputmvar(dbfile,'disp',4,N_proc,dispnames,     &
                        &   ldispnames,vartypes,DB_F77NULL,ierr)
              !write(*,*) 'err disp: ',err

           ! Restore the previous value for maximum string length
              err = dbset2dstrlen(oldlen)
           
           ! Close the Silo file
              ierr = dbclose(dbfile)
!!
           end if
        end if !end silo on!
     end if


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do  !end time loop!


!=========================================================================
!

! Processor Time
  t2 = MPI_WTIME()
  if (Pid == 0) then
     write(*,*)' Run time: ',t2-t1,' seconds; ',(t2-t1)/60,' mins'
  end if

end if !end transducer contact if!

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_Type_free(stride,ierr)
  call MPI_Type_free(stride1,ierr)
  call MPI_Comm_free(MPI_COMM_2D,ierr)


! close plotmtv file.
  if (Pid == 0) then
     write(112,*)"$ END"
     close(112)   
  end if  

  213 format (42I3) 
 
  deallocate(rxStart)
  deallocate(rxEnd)
  deallocate(ryStart)
  deallocate(ryEnd)
  deallocate(rzStart)
  deallocate(rzEnd)
  
  deallocate(Txx)
  deallocate(Tyy)
  deallocate(Tzz)
  deallocate(Txy)
  deallocate(Txz)
  deallocate(Tyz)
 
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  
  deallocate(disp)
  deallocate(catchAline)
  
  deallocate(tkwin_vect)
  deallocate(df)

  deallocate(xx)
  deallocate(yy)
  deallocate(zz)
  deallocate(dims)
  deallocate(lmeshnames)
  deallocate(meshtypes)
  deallocate(meshnames)
  deallocate(ldispnames)
  deallocate(vartypes)
  deallocate(dispnames)

end Program code

