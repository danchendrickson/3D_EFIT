Module exchange

Contains

  Subroutine exchange_vx(vx,lx,ly,numz,iz,stride,stride1,       &
        &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE

    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: vx 
    
    Integer :: ierr, status(MPI_STATUS_SIZE)

!!  from right to left   
    
!        call MPI_Sendrecv(vx(1,1,iz),1,stride,nbrleft,1,  &
!             &    vx(lx,1,iz),1,stride,nbrright,1, &
!             &    comm2d,status,ierr)       

        call MPI_Send(vx(1,1,iz),1,stride,nbrleft,1,comm2d,ierr)        
        call MPI_Recv(vx(lx,1,iz),1,stride,nbrright,1,comm2d,status,ierr)
             
!! from bottom to top 
!        call MPI_Sendrecv(vx(1,ly,iz),1,stride1,nbrbottom,2,   &
!             &    vx(1,1,iz),1,stride1,nbrtop,2,   &
!             &    comm2d,status,ierr)

        call MPI_Send(vx(1,ly,iz),1,stride1,nbrbottom,2,comm2d,ierr)
        call MPI_Recv(vx(1,1,iz),1,stride1,nbrtop,2,comm2d,status,ierr) 

  end Subroutine exchange_vx

  Subroutine exchange_vy(vy,lx,ly,numz,iz,stride,stride1,       &
        &nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE

    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: vy 
    
    Integer :: ierr, status(MPI_STATUS_SIZE)

!!  from right to left   
! Update vy 
!! from left to right         
!        call MPI_Sendrecv(vy(lx,1,iz),1,stride,nbrright,3, &
!             &    vy(1,1,iz),1,stride,nbrleft,3, &
!             &    comm2d,status,ierr)
     
        call MPI_Send(vy(lx,1,iz),1,stride,nbrright,3,comm2d,ierr)
        call MPI_Recv(vy(1,1,iz),1,stride,nbrleft,3,comm2d,status,ierr)
             

!! from top to bottom  
!        call MPI_Sendrecv(vy(1,1,iz),1,stride1,nbrtop,4,   &
!             &    vy(1,ly,iz),1,stride1,nbrbottom,4,  &
!             &    comm2d,status,ierr)
      
        call MPI_Send(vy(1,1,iz),1,stride1,nbrtop,4,comm2d,ierr)
        call MPI_Recv(vy(1,ly,iz),1,stride1,nbrbottom,4,comm2d,status,ierr) 

  end Subroutine exchange_vy

  Subroutine exchange_vz(vz,lx,ly,numz,iz,stride,stride1,       &
        &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE

    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: vz
    
    Integer :: ierr, status(MPI_STATUS_SIZE)

!!  from right to left   
! Update vz 

! from left to right         
!        call MPI_Sendrecv(vz(lx,1,iz),1,stride,nbrright,5, &
!             &    vz(1,1,iz),1,stride,nbrleft,5, &
!             &    comm2d,status,ierr)

        call MPI_Send(vz(lx,1,iz),1,stride,nbrright,5,comm2d,ierr)        
        call MPI_Recv(vz(1,1,iz),1,stride,nbrleft,5,comm2d,status,ierr)
             
! from bottom to top 
!        call MPI_Sendrecv(vz(1,ly,iz),1,stride1,nbrbottom,6,   &
!             &    vz(1,1,iz),1,stride1,nbrtop,6,  &
!             &    comm2d,status,ierr)
      
        call MPI_Send(vz(1,ly,iz),1,stride1,nbrbottom,6,comm2d,ierr)
        call MPI_Recv(vz(1,1,iz),1,stride1,nbrtop,6,comm2d,status,ierr) 

  end Subroutine exchange_vz

  Subroutine exchange_Txx(Txx,lx,ly,numz,iz,stride,stride1,     &
        &nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Txx
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Txx 
! 
! from left to right         
!    call  MPI_Sendrecv(Txx(lx,1,iz),1,stride,nbrright,7, &
!         &    Txx(1,1,iz),1,stride,nbrleft,7, &
!         &    comm2d,status,ierr)

    call MPI_Send(Txx(lx,1,iz),1,stride,nbrright,7,comm2d,ierr)        
    call MPI_Recv(Txx(1,1,iz),1,stride,nbrleft,7,comm2d,status,ierr)
                 
! from bottom to top 
!    call  MPI_Sendrecv(Txx(1,ly,iz),1,stride1,nbrbottom,8,   &
!         &    Txx(1,1,iz),1,stride1,nbrtop,8,  &
!         &    comm2d,status,ierr)

    call MPI_Send(Txx(1,ly,iz),1,stride1,nbrbottom,8,comm2d,ierr)
    call MPI_Recv(Txx(1,1,iz),1,stride1,nbrtop,8,comm2d,status,ierr) 
    
  end Subroutine exchange_Txx

  Subroutine exchange_Tyy(Tyy,lx,ly,numz,iz,stride,stride1,     &
        &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Tyy
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Tyy
!
! from left to right         
!    call  MPI_Sendrecv(Tyy(lx,1,iz),1,stride,nbrright,9, &
!         &    Tyy(1,1,iz),1,stride,nbrleft,9, &
!         &    comm2d,status,ierr)

    call MPI_Send(Tyy(lx,1,iz),1,stride,nbrright,9,comm2d,ierr)        
    call MPI_Recv(Tyy(1,1,iz),1,stride,nbrleft,9,comm2d,status,ierr)
                 
! from bottom to top 
!        call  MPI_Sendrecv(Tyy(1,ly,iz),1,stride1,nbrbottom,10,   &
!             &    Tyy(1,1,iz),1,stride1,nbrtop,10,  &
!             &    comm2d,status,ierr)

    call MPI_Send(Tyy(1,ly,iz),1,stride1,nbrbottom,10,comm2d,ierr)
    call MPI_Recv(Tyy(1,1,iz),1,stride1,nbrtop,10,comm2d,status,ierr) 
    
  end Subroutine exchange_Tyy

  Subroutine exchange_Tzz(Tzz,lx,ly,numz,iz,stride,stride1,     &
        &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Tzz
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Tzz 
! 
! from left to right         
!    call  MPI_Sendrecv(Tzz(lx,1,iz),1,stride,nbrright,11, &
!         &    Tzz(1,1,iz),1,stride,nbrleft,11, &
!         &    comm2d,status,ierr)
    
    call MPI_Send(Tzz(lx,1,iz),1,stride,nbrright,11,comm2d,ierr)        
    call MPI_Recv(Tzz(1,1,iz),1,stride,nbrleft,11,comm2d,status,ierr)

! from bottom to top 
!    call  MPI_Sendrecv(Tzz(1,ly,iz),1,stride1,nbrbottom,12,   &
!             &    Tzz(1,1,iz),1,stride1,nbrtop,12,  &
!             &    comm2d,status,ierr)                 

    call MPI_Send(Tzz(1,ly,iz),1,stride1,nbrbottom,12,comm2d,ierr)
    call MPI_Recv(Tzz(1,1,iz),1,stride1,nbrtop,12,comm2d,status,ierr) 
    
  end Subroutine exchange_Tzz
  
  Subroutine exchange_Txy(Txy,lx,ly,numz,iz,stride,stride1,     &
        &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Txy
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Txy
! 
! from right to left     
!    call  MPI_Sendrecv(Txy(1,1,iz),1,stride,nbrleft,13,  &
!         &    Txy(lx,1,iz),1,stride,nbrright,13, &
!         &    comm2d,status,ierr)
    
    call MPI_Send(Txy(1,1,iz),1,stride,nbrleft,13,comm2d,ierr)        
    call MPI_Recv(Txy(lx,1,iz),1,stride,nbrright,13,comm2d,status,ierr)

! from top to bottom 
!        call  MPI_Sendrecv(Txy(1,1,iz),1,stride1,nbrtop,14,   &
!             &    Txy(1,ly,iz),1,stride1,nbrbottom,14,   &
!             &    comm2d,status,ierr)
              
    call MPI_Send(Txy(1,1,iz),1,stride1,nbrtop,14,comm2d,ierr)
    call MPI_Recv(Txy(1,ly,iz),1,stride1,nbrbottom,14,comm2d,status,ierr) 
    
  end Subroutine exchange_Txy

  Subroutine exchange_Txz(Txz,lx,ly,numz,iz,stride,stride1,     &
            &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Txz
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Txz
! 
! from right to left     
!        call  MPI_Sendrecv(Txz(1,1,iz),1,stride,nbrleft,15,  &
!             &    Txz(lx,1,iz),1,stride,nbrright,15, &
!             &    comm2d,status,ierr)

    call MPI_Send(Txz(1,1,iz),1,stride,nbrleft,15,comm2d,ierr)        
    call MPI_Recv(Txz(lx,1,iz),1,stride,nbrright,15,comm2d,status,ierr)

! from top to bottom 
!        call  MPI_Sendrecv(Txz(1,1,iz),1,stride1,nbrtop,16,   &
!             &    Txz(1,ly,iz),1,stride1,nbrbottom,16,   &
!             &    comm2d,status,ierr)

    call MPI_Send(Txz(1,1,iz),1,stride1,nbrtop,16,comm2d,ierr)
    call MPI_Recv(Txz(1,ly,iz),1,stride1,nbrbottom,16,comm2d,status,ierr) 
    
  end Subroutine exchange_Txz


  Subroutine exchange_Tyz(Tyz,lx,ly,numz,iz,stride,stride1,     &
            &   nbrleft,nbrright,nbrbottom,nbrtop,comm2d)

    Implicit NONE
    
    include "mpif.h"

    Integer, Intent(in) :: lx,ly,numz,iz
    Integer, Intent(in) :: stride,stride1
    Integer, Intent(in) :: nbrleft,nbrright,nbrbottom,nbrtop
    Integer, Intent(in) :: comm2d
    
    Double Precision,  dimension(1:lx,1:ly,1:numz), Intent(out) :: Tyz
    
    Integer :: ierr, status(MPI_STATUS_SIZE)
!
! Update Txy
! 
! from right to left         
!        call  MPI_Sendrecv(Tyz(1,1,iz),1,stride,nbrleft,17, &
!             &    Tyz(lx,1,iz),1,stride,nbrright,17, &
!             &    comm2d,status,ierr)

    call MPI_Send(Tyz(1,1,iz),1,stride,nbrleft,17,comm2d,ierr)        
    call MPI_Recv(Tyz(lx,1,iz),1,stride,nbrright,17,comm2d,status,ierr)

! from bottom to top 
!        call  MPI_Sendrecv(Tyz(1,1,iz),1,stride1,nbrtop,18,   &
!             &    Tyz(1,ly,iz),1,stride1,nbrbottom,18,  &
!             &    comm2d,status,ierr)
              
    call MPI_Send(Tyz(1,1,iz),1,stride1,nbrtop,18,comm2d,ierr)
    call MPI_Recv(Tyz(1,ly,iz),1,stride1,nbrbottom,18,comm2d,status,ierr) 
    
  end Subroutine exchange_Tyz

end Module exchange

