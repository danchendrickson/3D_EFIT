Subroutine MPE_Decomp1d(n,Number_procs,Processor_ID,start,last,lsize)

  Integer, Intent(in) :: n
  Integer, Intent(in) :: Number_procs
  Integer, Intent(in) :: Processor_ID
  Integer, Intent(out) :: start
  Integer, Intent(out) :: last
  Integer, Intent(out) :: lsize  

  Integer :: nlocal
  Integer :: deficit
  Integer :: startl
  Integer :: endl

  nlocal = n/Number_procs
  startl =  1 + Processor_ID*(nlocal)

  deficit = mod(n,Number_procs)

  startl = startl + min(Processor_ID,deficit)

  if (Processor_ID < deficit) then
     nlocal = nlocal + 1
  end if

  endl = startl + (nlocal -1)  



!         
!  endl and startl are set. Now determine local size
!  not this is fortran so it is 1 to N 
!  we want to know N-1 = N-1 we want (N+1) - 1 = N
!  add to for the size of the ghostboundries  
!  add (2) by default then modify if node is 
!  on boundary.
!
  if (Processor_ID > 0) then
     startl = startl -1
  end if
     
  lsize = (endl+1)-startl
  
!
!  'n' is the global number of points
!
  

!
! This is only true if you have 1xN or Nx1 processor topology 
!
!  if ((startl == 1) .AND. (endl == n)) then
!     lsize = n
!  end if

!!  return start and last values

  start = startl
  last  = endl
  
end Subroutine MPE_Decomp1d

