/*
'setup_cart_space.h'
Sets up the transducers and reflectors for the Cartesian EFIT simulation
Cleaned up and modified version of 'spipe.h' (Bertoncini, Campbell-Leckey, Miller, etc)
Includes code to generate Rayleigh waves via IDTs (Miller)

Eric A. Dieckman (WM)
Last edited: 21 Sept 2011 EAD
*/

#include <iostream>
#include "array3D.h"
#include "array3D_int.h"
#include "transducer.h"
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

class setup_cart_space
{

public:
  setup_cart_space()  { }
  ~setup_cart_space() { }
	
  int num1;        // number of grid points in x direction
  int num2;        // number of grid points in y direction
  int num3;        // number of grid points in z direction
    
  int abc;         // number of abc points on each end 

  double ds;       // spatial step size (meters)
  double dt;       // time step size (seconds)

  double den;      // density (kg/m^3)
  double lm;       // lame constant - lamda
  double mu;       // lame constant - mu
 
  double xbeg;     // x start position (meters)

  int pipetype;    // pipe type  1 = left end , 2 = middle, 3 = right end 

  array3D v1;      // x velocities
  array3D v2;      // y velocities
  array3D v3;      // z velocities
  array3D T11;     // xx normal stress 
  array3D T22;     // yy normal stress 
  array3D T33;     // zz normal stress
  array3D T12;     // yx sheer stress
  array3D T23;     // yz sheer stress
  array3D T13;     // xz sheer stress
    
  array3D muu;     // mu (shear moduli)
  array3D lmd;     // first lame parameter
  array3D d;       // density

  array3D_int B;   // boundary array
  
  int time; 

  transducer *trans;
  int numtrans;
  
private:
  double dtods;
  double lmdtods;
  double l2mdtods;
  double mdtods;
  
  int x,y,z,pp1,pm1;
  double ro,ri,rr,co,ci,cc;
    
  double PIo2;

public:
// ===================================================================================
// Initialize
// ===================================================================================  
  void Init()
  {
    v1.Init(num1,num2,num3);
    v2.Init(num1,num2,num3);
    v3.Init(num1,num2,num3);
    T11.Init(num1,num2,num3);
    T22.Init(num1,num2,num3);
    T33.Init(num1,num2,num3);
    T12.Init(num1,num2,num3);
    T23.Init(num1,num2,num3);
    T13.Init(num1,num2,num3);
    B.Init(num1+2,num2+2,num3+2, pipetype);
    d.Init(num1,num2,num3);
    muu.Init(num1,num2,num3);
    lmd.Init(num1,num2,num3);

    dtods   = (dt)/(ds);
    lmdtods  = (lm*dt)/ds;
    l2mdtods = ((lm+2*mu)*dt)/ds;
    mdtods   = (mu*dt)/ds;

    PIo2 = 3.14159265358979/2;
  
    numtrans=0;
    time = 0;
    abc = 80;
  }


  void setmaterialprops()
  {  
    int x,y,z;
    for(x=0;x<num1;x++) 
      for(y=0;y<num2;y++)
	for(z=0;z<num3;z++)
        {
	  muu.set(x,y,z,mu);
	  d.set(x,y,z,den);
	  lmd.set(x,y,z,lm);
        }    
  }
   

// ===================================================================================
// Update velocities and stresses
// ===================================================================================
  void UpdateVs(int xs, int xend)
  {
    for (x = xs; x <= xend; x++)
    {
      v2.setindx(x,0,0);  v1.setindx(x,0,0);  v3.setindx(x,0,0);
      T22.setindx(x,0,0); T11.setindx(x,0,0); T33.setindx(x,0,0);
      T12.setindx(x,0,0); T23.setindx(x,0,0); T13.setindx(x,0,0);
      d.setindx(x,0,0);muu.setindx(x,0,0);lmd.setindx(x,0,0);

      for (y = 0; y < num2; y++)
      {
	B.setindx(x+1,y+1,1);
	for (z = 0; z < num3; z++)
	{     

	  // --- v2 ---
          if (B.v()==0)  
	  { 
	    v2.sv( v2.v() + 2*dtods/(d.v()+d.v2p())*((T12.v()-T12.v1m())+(T22.v2p()-T22.v())+(T23.v()-T23.v3m()))); 
	  }
	  else if (B.v()==2 | B.v1m()==2 | B.v3m()==2) {} // because v2 requires elements in direction v3m() and v1m()
	  //vr.sv( vr.v() + dtodsp*(trans[B.v()-1000].drivef(time)+(1/rr)*((1/(2*(rr)))*(3*Tpp.v()-Tpp.v2m()))) );
          else if (B.v2p()==2)  
	  { 
	    v2.sv( v2.v() + 2*dtods/(d.v()+d.v())*(-2*(T22.v()))); 
	  }
	  else if (B.v2m()==2)  
	  { 
	    v2.sv( v2.v() + 2*dtods/(d.v()+d.v())*(2*(T22.v2p()))); 
	  }
	  else
	  { 
	    v2.sv( v2.v() + 2*dtods/(d.v()+d.v2p())*((T12.v()-T12.v1m())+(T22.v2p()-T22.v())+(T23.v()-T23.v3m())));
	  }
	                           
	  
	  // --- v1 ---
	  if (B.v()==0)   
	  { 
	    v1.sv( v1.v() + 2*dtods/(d.v()+d.v1p())*((T11.v1p()-T11.v())+(T12.v()-T12.v2m())+(T13.v()-T13.v3m())));  // middle	    
	  }
          else if (B.v()==2 | B.v2m()==2 | B.v3m()==2) {}  //b/c v1 needs elements in the v2m() and v3m() directions
	  else if (B.v1p()==2)   
	  { 
	    v1.sv( v1.v() + 2*dtods/(d.v()+d.v())*(-2*(T11.v()))); //right end (only called if at last node)
	  }
          else if (B.v1m()==2)  
	  { 
	    v1.sv( v1.v() + 2*dtods/(d.v()+d.v())*(2*(T11.v1p()))); //left end (only called if at first node)
	  }
	  else 
	  { 
	    v1.sv( v1.v() + 2*dtods/(d.v()+d.v1p())*((T11.v1p()-T11.v())+(T12.v()-T12.v2m())+(T13.v()-T13.v3m()))); 
	  }


	  // --- v3 ---
	  if (B.v()==0)   	  
	  {  
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v3p())*((T13.v()-T13.v1m())+(T23.v()-T23.v2m())+(T33.v3p()-T33.v())) );
	  }
	  else if (B.v()==2 | B.v2m()==2 | B.v1m()==2) {}
	  else if ( B.v()>=9000)   
	  { 
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(-2*(T33.v()))+ 2*dtods/(den+den)*(-trans[B.v()-9000].drivef(time))); 
	  }
          else if ( B.v()>=1000)   
	  { 
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(-2*(T33.v()))+ 2*dtods/(den+den)*(trans[B.v()-1000].drivef(time))); 
	  }
          else if (B.v3p()==2)    
	  { 
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(-2*(T33.v()))); //top
	  }
          else if (B.v3m()==2)   
	  { 
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(2*(T33.v3p()))); //bottom
	  }
          else                  
	  { 
	    v3.sv( v3.v() + 2*dtods/(d.v()+d.v3p())*((T13.v()-T13.v1m())+(T23.v()-T23.v2m())+(T33.v3p()-T33.v())) );
	  }

	  
	  B.incindx(); v2.incindx();  v1.incindx();  v3.incindx();
          T22.incindx(); T11.incindx(); T33.incindx(); T12.incindx(); T23.incindx(); T13.incindx();
          d.incindx(); lmd.incindx();muu.incindx();
	}
      }
    }
  }


  void UpdateTs(int xs, int xend)
  {
    for (x = xs; x <= xend; x++)
    {
      v2.setindx(x,0,0);  v1.setindx(x,0,0);  v3.setindx(x,0,0);
      T22.setindx(x,0,0); T11.setindx(x,0,0); T33.setindx(x,0,0);
      T12.setindx(x,0,0); T23.setindx(x,0,0); T13.setindx(x,0,0);
      d.setindx(x,0,0);muu.setindx(x,0,0);lmd.setindx(x,0,0);
	
      for (y = 0; y < num2; y++)
      {
	B.setindx(x+1,y+1,1);
	for (z = 0; z < num3; z++)
	{
	  // --- Tii ---
	  if (B.v()==2 | B.v2m()==2 | B.v3m()==2 | B.v1m()==2) {} // 1 is x, 2 is y, 3 is z
          else    
	  {
	    T11.sv(T11.v()+dtods*(lmd.v()+(2*muu.v()))*(v1.v()-v1.v1m())+dtods*(lmd.v())*((v2.v()-v2.v2m())+(v3.v()-v3.v3m())));
            T22.sv(T22.v()+dtods*(lmd.v()+(2*muu.v()))*(v2.v()-v2.v2m())+dtods*(lmd.v())*((v1.v()-v1.v1m())+(v3.v()-v3.v3m())));
            T33.sv(T33.v()+dtods*(lmd.v()+(2*muu.v()))*(v3.v()-v3.v3m())+dtods*(lmd.v())*((v1.v()-v1.v1m())+(v2.v()-v2.v2m())));
	  }   
          
          // --- T13 --- 
          if (B.v()==2 | B.v1p()==2 | B.v1m()==2 | B.v3p()==2 | B.v3m()==2 | B.v2m()==2 | B.v13p()==2) {} // need all terms needed by v1 and v3 as well as any additional terms
          else 
	  {
	    T13.sv(T13.v()+dtods*(4/(((1/muu.v())+(1/muu.v1p())+(1/muu.v3p())+(1/muu.v13p()))))*((v1.v3p()-v1.v())+(v3.v1p()-v3.v())));
	  }

	  // --- T23 ---
	  if (B.v()==2 | B.v2p()==2 | B.v2m()==2 | B.v3p()==2 | B.v3m()==2 | B.v1m()==2 | B.v23p()==2) {} 
          else 
	  {
	    T23.sv(T23.v()+dtods*(4/((1/(muu.v())+(1/muu.v2p())+(1/muu.v3p())+(1/muu.v23p()))))*((v2.v3p()-v2.v())+(v3.v2p()-v3.v())));
	  }
                  
	  // --- T12 ---
	  if (B.v()==2 | B.v2p()==2 | B.v2m()==2 | B.v1p()==2 | B.v1m()==2 | B.v3m()==2 | B.v12p()==2) {} 
          else 
	  {
	    T12.sv(T12.v()+dtods*((4/((1/muu.v())+(1/muu.v1p())+(1/muu.v2p())+(1/muu.v12p()))))*((v1.v2p()-v1.v())+(v2.v1p()-v2.v())));
	  }
 
	  B.incindx(); v2.incindx();  v1.incindx();  v3.incindx();
          T22.incindx(); T11.incindx(); T33.incindx(); T12.incindx(); T23.incindx(); T13.incindx();
          d.incindx(); lmd.incindx();muu.incindx();
	}
      }
    }
  }

// ===================================================================================
// Add/update transducers (both original and IDT methods below)
// ===================================================================================
  void UpdateTransducers(int t)
  { 
    int tr;
    for (int i1 = 1; i1<num1-1; i1++)
      for (int i2 = 0; i2<num2-1; i2++)
	if (B.val(i1+1,i2+1,num3) >= 1000) 
	{
	  tr = B.val(i1+1,i2+1,num3)-1000; //transducer ID
	  trans[tr].record[t] =  trans[tr].record[t] + v3.val(i1,i2,num3-1); //record holds recorded info
	}
  }


  void addTransducer(transducer t)
  { 
        
    // --- Original way, with a circular transducer ---
//    /*
    numtrans = numtrans+1;
    transducer *temp = new transducer[numtrans]; // of class transducer, length numtrans
    
    for (int i = 0; i<numtrans-1; i++) // loops from 0 to actual number of transducers
    {
      temp[i] = trans[i];
    }
    
    temp[numtrans-1] = t;
    trans = temp;
    int nnodes = 0;
    
    for (int i1 = 1; i1<num1-1; i1++)
      for (int i2 = 0; i2<num2-1; i2++)
	if (((i1+xbeg-trans[numtrans-1].posi1)*(i1+xbeg-trans[numtrans-1].posi1)+(i2-trans[numtrans-1].posi2)*(i2-trans[numtrans-1].posi2)) <= (trans[numtrans-1].radius*trans[numtrans-1].radius))
	{  
	  B.set(i1+1,i2+1,t.posi3+1,1000+numtrans-1); //executes if (x-xo)^2+(y-yo)^2<rad^2, ie, if within a circle of transducer radius then set B >1000
	  nnodes++;
	} 
	
	//if (t.posi1-xbeg+1 >0) B.set(t.posi1-xbeg+1,t.posi2+1,t.posi3+1,1000+numtrans-1);
	temp[numtrans-1].numnodes = nnodes;
//    */
	  
	  
    // --- IDTs (generates Rayleigh waves) ---
    /*
    numtrans = numtrans+1;
    transducer *temp = new transducer[numtrans];//of class transducer, length numtrans
    for (int i = 0; i<numtrans-1; i++)//loops from 0 to actual number of transducers
    {
      temp[i] = trans[i];
    }

    temp[numtrans-1] = t;
    trans = temp;
    int nnodes = 0;
    int pos_left;
    int pos_right;
    int neg_left;
    int neg_right;
    int wlength = 16; // 16 grid points per wavelength
    int aperture = 10*wlength; // Ten wavelengths in the grid
    int top_y = trans[numtrans-1].posi2;
    int bot_y = top_y - aperture;
    
    for (int i1 = 1; i1<num1-1; i1++) 
    {
      // cout << "i1 = " << i1 << " and i1+xbeg-1 = " << i1+xbeg-1 << endl;  // Debugging only
      for (int i2 = 0; i2<num2-1; i2++)
	for (int i = -5; i <= 4; i++) // Ten finger pairs, centered at transducer_x
	{
	  pos_left  = trans[numtrans-1].posi1 + wlength*(i);
          pos_right = pos_left + wlength/4;
          neg_left  = pos_right + wlength/4;
          neg_right = neg_left + wlength/4;
	  
	  if ((i1+xbeg >= pos_left) && (i1+xbeg <= pos_right) && (i2 >= bot_y) && (i2 <= top_y)) 
	  {
	    B.set(i1+1,i2+1,t.posi3+1,1000+numtrans-1);//sets B value to 1000 or higher for x and y and z positions of positive transducer
	    nnodes++;
          }

	  else if ((i1+xbeg >= neg_left) && (i1+xbeg <= neg_right) && (i2 >= bot_y) && (i2 <= top_y)) 
	  {
	    B.set(i1+1,i2+1,t.posi3+1,9000+numtrans-1);//sets B value to 9000 or higher for x and y and z positions of negative transducer
	    nnodes++;
	  }
	}
    }

    temp[numtrans-1].numnodes = nnodes;
    */
  }

// ===================================================================================
// Add/update relectors - different types defined below
// ===================================================================================
  void addReflector(double typ, double p1, double p2, int start3, int end3, double rad, double dd, double mu, double lambda)
  {
    // --- VERTICAL 2D CRACK ---
    if (typ == 0)
    {
      // Below is for easy of use only, no real purpose in defining new variables
      int x1 = p2;
      int y1 = start3;
      int x2 = end3;
      int y2 = rad;
      int depth  = dd; // Crack depth in z direction
            
      cout<<"Crack set from p1 = ["<<x1<<" "<<y1<<"] to p2 = ["<<x2<<" "<<y2<<"] at a depth of z = "<<depth<<endl;

      double slope;
      int y_int;
      int x;
      int y;
      int begin;
      int finish;
      int temp;
            
      if (x2 == x1) 
      {
	slope = 9999; // approximates infinity
      }

      else 
      {
	slope = (double)(y2-y1)/(x2-x1);
      }
      
      y_int = y1 - slope*x1;
      cout<<" Slope = "<<slope<<" and y_int = "<<y_int<<endl;
      
      if ((x1-x2)*(x1-x2) > (y1-y2)*(y1-y2)) // Fixes stretching
      {
	temp = 1;
	if (x2 > x1) 
	{
	  begin = x1;
          finish = x2;
	}

	else 
	{
	  begin = x2;
          finish = x1;
	}
	
	// cout<<"Restricting to x begin = "<<begin<<" to finish = "<<finish<<endl;
      }

      else 
      {
	temp = 2;
        if (y2 > y1) 
	{
	  begin = y1;
          finish = y2;
	}
	
	else 
	{
	  begin = y2;
          finish = y1;
	}
	
	// cout<<"Restricting to y begin = "<<begin<<" to finish = "<<finish<<endl;
      }

      for (int i1 = 0; i1 < num1; i1++)
	for (int i2 = 0; i2 < num2; i2++)
	  for (int i3 = 0; i3 < depth; i3++) 
	  {
	    if (temp == 1) 
	    {
	      y = round(slope*i1+xbeg-1 + y_int);
	      if ((i2 == y) && (i1+xbeg-1 >= begin) && (i1+xbeg-1 <= finish)) 
	      {
            	// cout<<"Setting crack at ["<<i1+xbeg-1<<" "<<i2<<" "<<i3<<"]!!"<<endl;
                B.set(i1,i2,i3,2);
              }
	    }

	    if (temp == 2) 
	    {
	      x = round((i2-y_int)/slope);
              if ((i1+xbeg-1 == x) && (i2 >= begin) && (i2 <= finish)) 
	      {
		// cout<<"Setting crack at ["<<i1+xbeg-1<<" "<<i2<<" "<<i3<<"]!!"<<endl;
		B.set(i1,i2,i3,2);
	      }
	    }
	  }  
    }                  
           
           
    // --- SPHERE ---
    else if (typ == 1)
    {
      for (int i1 = 0; i1 < num1; i1++)
	for (int i2 = 0; i2 < num2; i2++)
	  for (int i3 = 0; i3 < num3; i3++) 
	    if (((i1+xbeg-1-p1)*(i1+xbeg-1-p1) + (i2-p2)*(i2-p2) + (i3-start3)*(i3-start3)) < rad*rad)
	      if ((dd == -1))
	      {
		B.set(i1,i2,i3,2);//rigid scatterer
	      }
	      
	      else
	      {
		muu.set(i1,i2,i3,mu);
		d.set(i1,i2,i3,dd);
		lmd.set(i1,i2,i3,lambda);
	      }
    }


    // --- HALF SPACE ---
    else if (typ == 2)
    {
      for (int i1 = 0; i1 < num1; i1++)
	for (int i2 = 0; i2 < num2; i2++)
	  for (int i3 =0; i3 <= num3; i3++)  
	  {
	    if (i3<=round(p1*(i1+xbeg)))
	    {
	      cout<<"p1 is "<<p1<<"\n";
              cout<<"half space is "<<round(p1*(i1+xbeg))<<"\n";         
	      B.set(i1,i2,i3,2);
	    }
	  }
    }


    // --- 3D RECTANGULAR VOID ---
    else if (typ == 3)
    {
      for (int i1 = 0; i1 < num1; i1++)
	for (int i2 = 0; i2 < num2; i2++)
	  for (int i3 = 0; i3 < num3; i3++)  
	    if ((i1+xbeg-1 >= p1) && (i1+xbeg-1 <= p2) && (i2 >= start3) && (i2 <= end3) && (i3 >= rad) && (i3 <= dd))
	    {
	      B.set(i1,i2,i3,2); // Stress-free void
	    }
    } 
    
    
    // --- ROUNDED RECTANGLE (EXPERIMENTAL)
    else if (typ == 101)
    {
      //double p1, double p2, int start3, int end3, double rad, double dd, double mu, double lambda)
      //p1=startx
      //p2=lenx
      //start3=starty
      //end3=leny
      int startx=p1;
      int starty=start3;
      int endx=p1+p2;
      int endy=start3+end3;
      //rad=radius of spheres making up edges
      //dd=thickness
      
      for (int i1 = 0; i1 < num1; i1++)
	for (int i2 = 0; i2 < num2; i2++)
	  for (int i3 = 0; i3 < num3; i3++) 
	  {
	    if ((i1+xbeg-1)>(p1+rad) && (i1+xbeg-1)<(endx-rad) && i2>starty && i2<endy && i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2); // rigid scatterer, or create boundary
            }
            
            if ((i1+xbeg-1)>(p1) && (i1+xbeg-1)<(p1+rad) && i2>(starty+rad) && i2<(endy-rad) && i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
            }
            
            if ((i1+xbeg-1)>(endx-rad) && (i1+xbeg-1)<(endx) && i2>(starty+rad) && i2<(endy-rad) && i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
	    }

	    if (((i1+xbeg-1-(endx-rad))*(i1+xbeg-1-(endx-rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
            }   

	    if (((i1+xbeg-1-(startx+rad))*(i1+xbeg-1-(startx+rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
            }   

	    if (((i1+xbeg-1-(endx-rad))*(i1+xbeg-1-(endx-rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
            }   
            
            if (((i1+xbeg-1-(startx+rad))*(i1+xbeg-1-(startx+rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num3-dd))
	    {
	      B.set(i1,i2,i3,2);
            }
	  }
    }

}


  int round(double a)
  {
    return int(a+0.5);
  }

  
};
