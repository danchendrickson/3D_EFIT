#include <iostream>

class array3D
{
  private:
	int    ci;   // current index
	int    ci3;  // current i3 index (used for contenous Boundary)
	int    L2L3; // max2*max3

  public:
    double *a;
	int len1;           // number of grid points in r direction
	int len2;           // number of grid points in z direction
	int len3;           // number of grid points in p direction

	// ===================================================================================
	// Blank Constructer
    array3D() {}

	// ===================================================================================
    // Blank Deconstructor
    ~array3D() {}

	// ===================================================================================
    // Init - defines the array and its dimensions - MUST BE CALLED BEFORE USING
    void Init(int m1, int m2, int m3)
	{
      len1 = m1; len2 = m2; len3 = m3;
      L2L3 = m2*m3;
      a = new double[m1*m2*m3];
	  clear();
	  return;
	}
	
	 void Initt(int m1, int m2, int m3, double def)
	{
      len1 = m1; len2 = m2; len3 = m3;
      L2L3 = m2*m3;
      a = new double[m1*m2*m3];
	  cleart(def);
	  return;
	}

	// ===================================================================================
	// Return value at i1, i2, i3
	double val(int i1, int i2, int i3)
	{
      return a[(i1*L2L3)+(i2*len3)+i3];
	}

	// ===================================================================================
	// Set value at i1, i2, i3
	void set(int i1, int i2, int i3, double val)
	{
	  a[(i1*L2L3)+(i2*len3)+i3] = val;
	  return;
	}

	// ===================================================================================
	// quick access methods
    void   setindx(int i1, int i2, int i3) {	ci = (i1*L2L3)+(i2*len3)+i3; ci3 = i3; }
    void   incindx() 
	{	ci = ci+1; ci3 = ci3+1; 
	    if (ci3==len3) ci3=0;
	}
	void   sv(double x)                    {	a[ci] = x;     }  // set value at ci
    double v()    { return a[ci];        }   // equiv of a[i1][i2][i3]
    double v1p()  { return a[ci+L2L3];   }   // equiv of a[i1-1][i2][i3]
    double v1m()  { return a[ci-L2L3];   }   // equiv of a[i1+1][i2][i3]
    double v12p()  { return a[ci+L2L3+len3];   }   // equiv of a[i1+1][i2+1][i3]
    double v2p()  { return a[ci+len3];   }   // equiv of a[i1][i2+1][i3]
	double v2p2() { return a[ci+2*len3]; }   // equiv of a[i1][i2+2][i3]
    double v2m()  { return a[ci-len3];   }   // equiv of a[i1][i2-1][i3]
    double v3p()  {//   return a[ci+1]; }                       // equiv of a[i1][i2][i3+1]
		if (ci3 == len3-1)
			return a[ci-len3+1];
        else 
		    return a[ci+1];
	}   
    double v3m()  { //return a[ci-1]; }                         // equiv of a[i1][i2][i3-1]
		if (ci3 == 0)
			return a[ci+len3-1];
	    else 
			return a[ci-1];
	}  
	
	 double v13p()  {//   return a[ci+1]; }                       // equiv of a[i1][i2][i3+1]
		if (ci3 == len3-1)
			return a[ci-len3+1+L2L3];
        else 
		    return a[ci+1+L2L3];
	}   
	
	 double v23p()  {//   return a[ci+1]; }                       // equiv of a[i1][i2][i3+1]
		if (ci3 == len3-1)
			return a[ci-len3+1+len3];
        else 
		    return a[ci+1+len3];
	}   


    // ===================================================================================
    // clear - sets all values = 0;
	void clear()
	{
		for (int i = 0; i< L2L3*len1; i++)
			a[i] = 0;
	}
	
		void cleart(double def)
	{
		for (int i = 0; i< L2L3*len1; i++)
			a[i] = def;
	}

	// ===================================================================================
	// returns 2D slice through 3D array at fixed index 2
    double* slice_fix2(int i2)
	{
	   double *x = new double[(len1-2)*len3];
	   x[0] = (len1-2)*len3;

	   int c = 0;
	   for (int i1 = 1; i1<len1-1; i1++)  // does not return ends
         for (int i3 = 0; i3<len3; i3++)
		 {
			 x[c]=val(i1, i2, i3);
			 c++;
		 }

	   return x;
	}
	int slice_fix2_count() { return (len1-2)*len3;  }
	
	
		// ===================================================================================
	// returns 3D volume returning only the even indexes
    double* GetEvenVol(int start)
	{
	   int len = GetEvenVolLen(start);
	   double *x = new double[len];
	   
	   int c = 0;
	 //  for (int i1 = 1+(start%2); i1<len1-1; i1=i1+2)  // does not return ends
	  for (int i1 = 1+(start%2); i1<len1-1; i1=i1+2) 
		   for (int i2 = 0; i2<len2; i2=i2+2)  
              for (int i3 = 0; i3<len3; i3=i3+2)
			  {
			    // x[c]=(val(i1-1, i2, i3)+val(i1, i2, i3))/2;
			      x[c]=val(i1, i2, i3);
			     c++;
			  }
			  
	   //std::cout << len << " " << c-1 << "\n";


	   return x;
	}

	int GetEvenVolLen(int start)
	{
		int len;
		if (start%2 == 0) 
		   len = (len1-1)/2*(len2/2)*(len3/2);
	    else
           len = (len1-2)/2*(len2/2)*(len3/2);
		return len;
	}

	// ===================================================================================
	// returns 1D slice through 3D array at fixed index i2, i3
    double* slice_fix1(int i1, int i2)
	{
	   double *x = new double[(len3)];
	   x[0] = len3;

	   int c = 0;
	   for (int i3 = 0; i3<len3-1; i3++)  // does not return ends
        		 {
			 x[c]=val(i1, i2, i3);
			 c++;
		 }

	   return x;
	}
	int slice_fix1_count() { return len3;  }
	


};
