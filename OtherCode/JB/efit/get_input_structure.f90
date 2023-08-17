 Subroutine read_dmaterials(cmax,cmin,spacethickness,spacelength,   &
        &   spacewidth,dden,dcl,dcs,minz)

  Implicit none
 
  Double Precision, Intent(out) :: cmax
  Double Precision, Intent(out) :: cmin
  Double Precision, Intent(out) :: spacethickness 
  Double Precision, Intent(out) :: spacelength
  Double Precision, Intent(out) :: spacewidth
  Double Precision, Intent(out) :: dden
  Double Precision, Intent(out) :: dcl
  Double Precision, Intent(out) :: dcs
  Double Precision, Intent(out) :: minz
  
  open(200,FILE = "material_parameters",STATUS = "OLD",ACTION = "READ")
  read(200,*) cmax
  read(200,*) cmin
  read(200,*) spacethickness
  read(200,*) spacelength
  read(200,*) spacewidth  
  read(200,*) dden
  read(200,*) dcl
  read(200,*) dcs
  read(200,*) minz
  close(200)
  
  write(*,*)
  write(*,*) 'Material Properties for efit3d'
  write(*,*) 
  write(*,*) 'cmax = ',cmax
  write(*,*) 'cmin = ',cmin
  write(*,*) 'spacethickness = ',spacethickness
  write(*,*) 'spacelength = ',spacelength
  write(*,*) 'spacewidth = ',spacewidth
  write(*,*) 'dden = ',dden
  write(*,*) 'dcl = ',dcl
  write(*,*) 'dcs = ',dcs
  write(*,*) 'minz = ',minz

end Subroutine read_dmaterials

Subroutine read_region(whichR,rxS,rxE,ryS,ryE,rzS,rzE,denr,clr,csr)
  
  Implicit none
  
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
  
  Character(5) :: wR

  write(unit=wR,fmt='(I5)') whichR
  open(203,FILE = "region_"//(TRIM(ADJUSTL(wR))),STATUS = "OLD",    &
        &   ACTION = "READ")
  read(203,*) rxS
  read(203,*) rxE
  read(203,*) ryS
  read(203,*) ryE
  read(203,*) rzS
  read(203,*) rzE
  read(203,*) denr
  read(203,*) clr
  read(203,*) csr
  close(203)

!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'xStart = ',rxS
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'xEnd= ',rxE
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'yStart = ',ryS
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'yEnd= ',ryE
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'zStart = ',rzS
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'zEnd= ',rzE
!  write(*,*) 'r'//(TRIM(ADJUSTL(wR)))//'den = ',denr
  
end Subroutine read_region

Subroutine read_transducer(tposx,tposy,tposz,tthickness,tfreq,     &
    &   pcycles,catchtposx,catchtposy,catchtposz,catchtthickness,tmode)

  Implicit none

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

  open(201,FILE = "transducer_parameters",STATUS = "OLD",ACTION = "READ")
  read(201,*) tposx
  read(201,*) tposy
  read(201,*) tposz
  read(201,*) tthickness
  read(201,*) tfreq
  read(201,*) pcycles
  read(201,*) catchtposx
  read(201,*) catchtposy
  read(201,*) catchtposz
  read(201,*) catchtthickness
  read(201,*) tmode
  close(201)

  write(*,*)
  write(*,*) 'Transducer Properties for efit3d'
  write(*,*) 
  write(*,*) 'tposx = ',tposx
  write(*,*) 'tposy = ',tposy
  write(*,*) 'tposz = ',tposz
  write(*,*) 'tthickness = ',tthickness
  write(*,*) 'tfreq = ',tfreq
  write(*,*) 'pcycles = ',pcycles
  write(*,*) 'catchtposx = ',catchtposx 
  write(*,*) 'catchtposy = ',catchtposy
  write(*,*) 'catchtposz = ',catchtposz
  write(*,*) 'catchtthickness = ',catchtthickness
  write(*,*) 'tmode = ',tmode

end Subroutine read_transducer

Subroutine read_model(alpha,plate,numR,delam,numS,plotevery,    &
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

  open(202,FILE = "model_parameters",STATUS = "OLD",ACTION = "READ")
  read(202,*) alpha
  read(202,*) plate
  read(202,*) numR
  read(202,*) delam
  read(202,*) numS
  read(202,*) plotevery
  read(202,*) silo_on
  read(202,*) siloSpace
  read(202,*) fmax
  read(202,*) SimulationTime 
  read(202,*) Pdim0
  read(202,*) Pdim1
  close(202)

  write(*,*)
  write(*,*) 'Model parameters for efit3d'
  write(*,*) 
  write(*,*) 'alpha = ', alpha   
  write(*,*) 'plate = ', plate
  write(*,*) 'number of regions = ',numR   
  write(*,*) 'delamination = ',delam
  write(*,*) 'surface = ',numS
  write(*,*) 'plotevery = ',plotevery 
  write(*,*) 'silo_on = ' ,silo_on
  write(*,*) 'siloSpace = ' ,siloSpace
  write(*,*) 'fmax = ',fmax  
  write(*,*) 'SimulationTime = ',SimulationTime
  write(*,*) 'Processor Topology is (',Pdim0,' by ',Pdim1,')'
  write(*,*)
end Subroutine read_model


Subroutine read_surface(rxStart,rxEnd,ryStart,ryEnd,Pid,S)

  Implicit none
  
  Integer, Intent(in) :: rxStart,rxEnd,ryStart,ryEnd
  Integer, Intent(in) :: Pid
  Double Precision, dimension(rxStart:rxEnd,ryStart:ryEnd),Intent(out) :: S

  Character(5) :: pp

  write(unit=pp,fmt='(I5)') Pid
  open(205,FILE = "surf_"//TRIM(ADJUSTL(pp)),STATUS= "OLD",ACTION = "READ")
  read(205,*) S
  close(205)
  
!  write(*,*) shape(S)
!  write(*,*)'Pid:',Pid,'row1:',S(:,1)

!  if (Pid == 2) then
!     S = 12.0
!  else if (Pid == 4) then
!     S = 13.0
!  end if
!  S(10:103,10:103) = 14
!  S(40:70,40:70) = 6

!  if (Pid == 2) then
!     write(*,*) S
!  end if

end Subroutine read_surface

Subroutine read_T_region(tfyStart,tfyEnd,tfzStart,tfzEnd,twyStart,  &
        &   twyEnd,twzStart,twzEnd,ttyStart,ttyEnd,ttzStart,ttzEnd)
  
  Implicit none
  
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
  
  open(205,FILE = "region_T",STATUS = "OLD",ACTION = "READ")
  read(205,*) tfyStart
  read(205,*) tfyEnd
  read(205,*) tfzStart
  read(205,*) tfzEnd
  read(205,*) twyStart
  read(205,*) twyEnd
  read(205,*) twzStart
  read(205,*) twzEnd
  read(205,*) ttyStart
  read(205,*) ttyEnd
  read(205,*) ttzStart
  read(205,*) ttzEnd
  close(205)

  write(*,*) 'T flange yS,yE,zS,zE: ',tfyStart,tfyEnd,tfzStart,tfzEnd
  write(*,*) 'T web yS,yE,zS,zE: ', twyStart,twyEnd,twzStart,twzEnd
  write(*,*) 'T top yS,yE,zS,zE: ', ttyStart,ttyEnd,ttzStart,ttzEnd

end Subroutine read_T_region
