Subroutine tukeywin(w,N,alpha)

  Integer,Intent(in) :: N
  Double Precision,Intent(in) :: alpha
  Double Precision,Dimension(1:N),Intent(out) :: w(:) 
 
  Integer :: M
  Integer :: i
  Double Precision, Parameter :: pi = 3.1415926535 

  M = (N/2)* alpha
!  write(*,*)'M =',M
!  write(*,*)'N = ',N
!  write(*,*)'alpha = ',alpha

  do i = 1,N     
     if ((i < N-M) .and. (i > M)) then
        w(i) = 1
     else
        if (i <= M) then
           w(i) =  0.5 - 0.5*cos(((i-1)*pi/M))
        elseif (i >= N-M) then
           w(i)=  0.5 - 0.5*cos((pi*(N-i))/M)
        end if
     end if
  end do


  
end subroutine tukeywin
