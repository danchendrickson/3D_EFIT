Subroutine plot_pulse(dt,dfl,df,tkwin_vect)

  Integer , Intent(in) :: dfl
  Double Precision, Intent(in) :: dt
! tkwin_vect -  tukeywin transform
  Double Precision, dimension(1:dfl), Intent(in) :: tkwin_vect
! df -  pulse form
  Double Precision, dimension(1:dfl), Intent(in) :: df

  Integer :: i

! open output files
  
  open(300,FILE = "results/pulse.mtv",STATUS="UNKNOWN",ACCESS="SEQUENTIAL")
  write(300,*)"$ DATA=CURVE2D"
  write(300,*)" "
  write(300,*)"% XMIN=0"
  write(300,*)"% XMAX=10"
  write(300,*)"% YMIN=-2"
  write(300,*)"% YMAX=2"    
  write(300,*)" "
 
  do i = 1,dfl
     write(300,*) i*dt*1000000,  df(i)
  end do ! df -  pulse form
   write(300,*)" "
  do i = 1,dfl
     write(300,*) i*dt*1000000,  tkwin_vect(i)
  end do
 
  
  write(300,*)"$ END"
  close(300)   

end Subroutine plot_pulse