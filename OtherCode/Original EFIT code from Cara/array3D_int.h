class array3D_int
{
  private:
	int    *a;  
    int    ci;   // current index
	int    L2L3; // max2*max3
    int    endtype;

  public:

	int len1;           // number of grid points in z direction
	int len2;           // number of grid points in r direction
	int len3;           // number of grid points in p direction


	// Blank Constructer
    array3D_int() {}

    // Deconstructor
    ~array3D_int() {}

    // Init - defines the array and its dimensions - MUST BE CALLED BEFORE USING
    void Init(int m1, int m2, int m3, int type)
	{
      len1 = m1; len2 = m2; len3 = m3;
      L2L3 = m2*m3;
      a = new int[m1*m2*m3];
	  endtype = type;
	  clear();
	  return;
	}

	// Return value at i1, i2, i3
	int val(int i1, int i2, int i3)
	{
      return a[(i1*L2L3)+(i2*len3)+i3];
	}

	// Set value at i1, i2, i3
	void set(int i1, int i2, int i3, int val)
	{
	  a[(i1*L2L3)+(i2*len3)+i3] = val;
	  return;
	}

	// quick access methods
    void setindx(int i1, int i2, int i3) {	ci = (i1*L2L3)+(i2*len3)+i3; }
	void incindx() {	ci = ci+1; }
    int v()    { return a[ci];      }  // equiv of a[i1][i2][i3]
    int v1p()  { return a[ci+L2L3]; }  // equiv of a[i1-1][i2][i3]
    int v1m()  { return a[ci-L2L3]; }  // equiv of a[i1+1][i2][i3]
    int v2p()  { return a[ci+len3]; }  // equiv of a[i1][i2+1][i3]
    int v2m()  { return a[ci-len3]; }  // equiv of a[i1][i2-1][i3]
    int v3p()  { return a[ci+1];    }  // equiv of a[i1][i2][i3+1]
    int v3m()  { return a[ci-1];    }  // equiv of a[i1][i2][i3-1]
    int v13p()  {return a[ci+1+L2L3];}
    int v23p()  {return a[ci+1+len3];}
    int v12p() { return a[ci+L2L3+len3];}
    // clear - sets all values = 0;
	void clear()
	{
		//std::cout << "type " << endtype << "\n";
		for (int i = 0; i< L2L3*len1; i++)
			{	a[i] = 0;}
	    
        for (int i1 = 0; i1< len1; i1++)
          for (int i3 = 0; i3< len3; i3++)
		  {
			set(i1,0,i3,2);  //set x2 boundary
			set(i1,1,i3,1);  
			//set(i1,2,i3,1);  
            
            set(i1,len2-1,i3,2);
           	set(i1,len2-2,i3,2);
			set(i1,len2-3,i3,1);
         
}
		  
		  
		   for (int i1 = 0; i1< len1; i1++)
          for (int i2 = 0; i2< len2; i2++)
		  {
             
			set(i1,i2,0,2); 
			set(i1,i2,1,1);
		
			  set(i1,i2,len3-1,2);//set x3 boundary
             set(i1,i2,len3-2,2);
             
              if ((i2>0) && (i2<len1-2)){
              set(i1,i2,len3-3,1);  } //	set(i1,i2,2,1); }
            
             
		  }


		for (int i2 = 0; i2< len2; i2++)
         for (int i3 = 0; i3< len3; i3++)
		  {
			 if (endtype == 1)
			 {
                         
                set(0,i2,i3,2);  // if at actual end of space in z set boundary
				 //set(1,i2,i3,2);
                   if ( (i2 > 0) && (i2<(len2-2)) && (i3 > 0) && (i3<(len3-2))){
				
			 	 set(1,i2,i3,1);   }
			     
			 }
             if (endtype == 3)
			 {
                         
                         set(len1-1,i2,i3,2);
				  set(len1-2,i2,i3,2);
                          if ( (i2 > 0) && (i2<(len2-2)) && (i3 >0) && (i3<(len3-2))){
				 
    		    
	   		  set(len1-3,i2,i3,1);  
	   		     
			 }
          }
          }


}
	
};
