#include <iostream>
#include "array3D.h"
#include "array3D_int.h"
#include "transducer.h"
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

class spipe
{

  public:
	spipe()  { }
    ~spipe() { }
	
	int num2;        // number of grid points in r direction
	int num1;        // number of grid points in z direction
	int num3;        // number of grid points in p direction
    
    int abc;         // number of abc points on each end 

	double ds;       // spatial step size in r and z direction (meters)
	double dp;       // angular step size in p direction (radians)
	double dt;       // time step size (seconds)

	double den;      // density (kg/m^3)
    double lm;       // lame constant - lamda
	double mu;       // lame constant - mu

 
	double zbeg;     // z start position (meters)
	double rbeg;     // r start position (meters) !inner pipe radius!

    int pipetype;    // pipe type  1 = left end , 2 = middle, 3 = right end 

	array3D v2;      //  r  - velocities
	array3D v1;      //  z  - velocities
    array3D v3;      //  p  - velocities
    array3D T22;     //  rr - normal stress 
    array3D T11;     //  zz - normal stress 
    array3D T33;     //  pp - normal stress
    array3D T12;     //  rz - sheer stress
    array3D T23;     //  rp - sheer stress
    array3D T13;     //  zp - sheer stress
    
    array3D muu;    //mu (shear modulus) matrix
    array3D lmd;    // first lame parameter matrix
    array3D d;      // density matrix

    array3D_int B;   //  Boundary Array
  
	int time; 

    transducer *trans;
	int        numtrans;
  
  private:

    double dtods;
    double lmdtods;
    double l2mdtods;
    double mdtods;
  
	int r,z,p,pp1,pm1;
	double ro,ri,rr,co,ci,cc;
    
	double PIo2;

  public:


	void Init()
	{
       v2.Init(num1,num2,num3);
	   v1.Init(num1,num2,num3);
       v3.Init(num1,num2,num3);
       T22.Init(num1,num2,num3);
       T11.Init(num1,num2,num3);
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

  void setmaterialprops(){
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
   

	void UpdateVs(int zs, int zend)
	{       
	   

	   for (z = zs; z <= zend; z++)
       {
        v2.setindx(z,0,0);  v1.setindx(z,0,0);  v3.setindx(z,0,0);
         T22.setindx(z,0,0); T11.setindx(z,0,0); T33.setindx(z,0,0);
		 T12.setindx(z,0,0); T23.setindx(z,0,0); T13.setindx(z,0,0);
		 d.setindx(z,0,0);muu.setindx(z,0,0);lmd.setindx(z,0,0);
		 
		 for (r = 0; r < num2; r++)
		 {

		   B.setindx(z+1,r+1,1);

		   for (p = 0; p < num3; p++)
		   {     
             //v2----------
              
              if (B.v()==0)  { v2.sv( v2.v() + 2*dtods/(d.v()+d.v2p())*((T12.v()-T12.v1m())+(T22.v2p()-T22.v())+(T23.v()-T23.v3m()))); }
              else if (B.v()==2 | B.v1m()==2 | B.v3m()==2) {}// because v2 requires elements in direction v3m() and v1m()
              //vr.sv( vr.v() + dtodsp*(trans[B.v()-1000].drivef(time)+(1/rr)*((1/(2*(rr)))*(3*Tpp.v()-Tpp.v2m()))) );
			  else if (B.v2p()==2)  { v2.sv( v2.v() + 2*dtods/(d.v()+d.v())*(-2*(T22.v()))); }
			  else if (B.v2m()==2)  { v2.sv( v2.v() + 2*dtods/(d.v()+d.v())*(2*(T22.v2p()))); }
              else            { v2.sv( v2.v() + 2*dtods/(d.v()+d.v2p())*((T12.v()-T12.v1m())+(T22.v2p()-T22.v())+(T23.v()-T23.v3m()))); }
	                           
              
			 
              // v1 ----------
              if (B.v()==0)   { v1.sv( v1.v() + 2*dtods/(d.v()+d.v1p())*((T11.v1p()-T11.v())+(T12.v()-T12.v2m())+(T13.v()-T13.v3m()))); }//middle
              else if (B.v()==2 | B.v2m()==2 | B.v3m()==2) {}  //b/c v1 needs elements in the v2m() and v3m() directions
			  else if (B.v1p()==2)   { v1.sv( v1.v() + 2*dtods/(d.v()+d.v())*(-2*(T11.v()))); }//right end (only called if at last node)
              else if (B.v1m()==2)  { v1.sv( v1.v() + 2*dtods/(d.v()+d.v())*(2*(T11.v1p()))); }//left end (only called if at first node)
              else         { v1.sv( v1.v() + 2*dtods/(d.v()+d.v1p())*((T11.v1p()-T11.v())+(T12.v()-T12.v2m())+(T13.v()-T13.v3m()))); }

			  // v3 ----------
              if (B.v()==0)   	  {  v3.sv( v3.v() + 2*dtods/(d.v()+d.v3p())*((T13.v()-T13.v1m())+(T23.v()-T23.v2m())+(T33.v3p()-T33.v())) );}
              else if (B.v()==2 | B.v2m()==2 | B.v1m()==2) {}
              else if ( B.v()>=1000)   { v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(-2*(T33.v()))+ 2*dtods/(den+den)*(trans[B.v()-1000].drivef(time)));}
              else if (B.v3p()==2)    { v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(-2*(T33.v()))); }//top
              else if (B.v3m()==2)   { v3.sv( v3.v() + 2*dtods/(d.v()+d.v())*(2*(T33.v3p()))); }//bottom
              else                  { v3.sv( v3.v() + 2*dtods/(d.v()+d.v3p())*((T13.v()-T13.v1m())+(T23.v()-T23.v2m())+(T33.v3p()-T33.v())) );}


			 
			  B.incindx(); v2.incindx();  v1.incindx();  v3.incindx();
              T22.incindx(); T11.incindx(); T33.incindx(); T12.incindx(); T23.incindx(); T13.incindx();
               d.incindx(); lmd.incindx();muu.incindx();
           }
         }
	   }
	 }

     void UpdateTs(int zs, int zend)
	{       
       
	   for (z = zs; z <= zend; z++)
       {
        v2.setindx(z,0,0);  v1.setindx(z,0,0);  v3.setindx(z,0,0);
         T22.setindx(z,0,0); T11.setindx(z,0,0); T33.setindx(z,0,0);
		 T12.setindx(z,0,0); T23.setindx(z,0,0); T13.setindx(z,0,0);
		 d.setindx(z,0,0);muu.setindx(z,0,0);lmd.setindx(z,0,0);
		 
		 for (r = 0; r < num2; r++)
		 {
	 
		   B.setindx(z+1,r+1,1);
	       
		   for (p = 0; p < num3; p++)
		   {   
    	    
			 // Tii
             if (B.v()==2 | B.v2m()==2 | B.v3m()==2 | B.v1m()==2) {}// 1 is z, 2 is r, 3 is phi
             else    
               {
                T11.sv(T11.v()+dtods*(lmd.v()+(2*muu.v()))*(v1.v()-v1.v1m())+dtods*(lmd.v())*((v2.v()-v2.v2m())+(v3.v()-v3.v3m())));
                T22.sv(T22.v()+dtods*(lmd.v()+(2*muu.v()))*(v2.v()-v2.v2m())+dtods*(lmd.v())*((v1.v()-v1.v1m())+(v3.v()-v3.v3m())));
                T33.sv(T33.v()+dtods*(lmd.v()+(2*muu.v()))*(v3.v()-v3.v3m())+dtods*(lmd.v())*((v1.v()-v1.v1m())+(v2.v()-v2.v2m())));
               }   
                                            
             // T13 
             if (B.v()==2 | B.v1p()==2 | B.v1m()==2 | B.v3p()==2 | B.v3m()==2 | B.v2m()==2 | B.v13p()==2) {} //need all terms needed by v1 and v3 as well as any additional terms
             else {T13.sv(T13.v()+dtods*(4/(((1/muu.v())+(1/muu.v1p())+(1/muu.v3p())+(1/muu.v13p()))))*((v1.v3p()-v1.v())+(v3.v1p()-v3.v())));}
			    

		
			 
			 //  T23
             if (B.v()==2 | B.v2p()==2 | B.v2m()==2 | B.v3p()==2 | B.v3m()==2 | B.v1m()==2 | B.v23p()==2) {} 
               else {T23.sv(T23.v()+dtods*(4/((1/(muu.v())+(1/muu.v2p())+(1/muu.v3p())+(1/muu.v23p()))))*((v2.v3p()-v2.v())+(v3.v2p()-v3.v())));}
                  
             //   T12
             if (B.v()==2 | B.v2p()==2 | B.v2m()==2 | B.v1p()==2 | B.v1m()==2 | B.v3m()==2 | B.v12p()==2) {} 
             else {T12.sv(T12.v()+dtods*((4/((1/muu.v())+(1/muu.v1p())+(1/muu.v2p())+(1/muu.v12p()))))*((v1.v2p()-v1.v())+(v2.v1p()-v2.v())));}
                
                      
           
			  B.incindx(); v2.incindx();  v1.incindx();  v3.incindx();
              T22.incindx(); T11.incindx(); T33.incindx(); T12.incindx(); T23.incindx(); T13.incindx();
               d.incindx(); lmd.incindx();muu.incindx();

           }
         }
	   }
	 }



	void UpdateTransducers(int t)
	{ 
        
		int tr;

        for (int i1 = 1; i1<num1-1; i1++)
           for (int i2 = 0; i2<num2-1; i2++)
                  if (B.val(i1+1,i2+1,num3) >= 1000) 
			 { 
                  
				 tr = B.val(i1+1,i2+1,num3)-1000; //transducer ID
			
				 trans[tr].record[t] =  trans[tr].record[t] + v3.val(i1,i2,num3-1);//record holds recorded info
			
			 }
	}



	void addTransducer(transducer t)
	{ 
        
		numtrans = numtrans+1;
		transducer *temp = new transducer[numtrans];//of class transducer, length numtrans
		for (int i = 0; i<numtrans-1; i++)//loops from 0 to actual number of transducers
			temp[i] = trans[i];
		temp[numtrans-1] = t;
		trans = temp;
	//	cout<<"trans in add tranducer is: "<<trans;
        int nnodes = 0;
		//if (trans[numtrans-1].driven)
        for (int i1 = 1; i1<num1-1; i1++)
           for (int i2 = 0; i2<num2-1; i2++)
                 if (((i1+zbeg-trans[numtrans-1].posi1)*(i1+zbeg-trans[numtrans-1].posi1)+(i2-trans[numtrans-1].posi2)*(i2-trans[numtrans-1].posi2))<=(trans[numtrans-1].radius*trans[numtrans-1].radius)  )
			   {  //if line above says: if (x-xo)^2+(y-yo)^2<rad^2, ie, if within a circle of transducer radius then set B >1000
			//   cout<<"inside addTrans if \n";
				   B.set(i1+1,i2+1,t.posi3+1,1000+numtrans-1);//sets B value to 1000 or higher for x and y and z positions of transducer
				   //cout<<"position marked as trans is "<<t.posi3<<"\n";
				   
				//   cout<<"value B is set to "<<B.val(i1+1,i2+1,i3+1)<<"\n";
				   nnodes++;
			   }
		//if (t.posi1-zbeg+1 >0) B.set(t.posi1-zbeg+1,t.posi2+1,t.posi3+1,1000+numtrans-1);

        temp[numtrans-1].numnodes = nnodes;
        
	}




void addReflector(double typ, double p1, double p2, int start3, int end3, double rad, double dd, double mu, double lambda)
	{
                         
        if (typ == 0)  //crack 2d
        {
           for (int i1 = (dd); i1 <= (mu); i1++)
                for (int i3 =(start3); i3 <= lambda; i3++) 
                {
                    if (i3>=round(p1*i1+end3) && i3<=round(p1*i1+end3+4) ) {
                    //cout<<"making crack \n";
                    //    cout<<"p1 is "<<p1<<"\n";
                    // cout<<"half space is "<<round(p1*(i1+zbeg))<<"\n";        
                    B.set(i1,996,i3,2);
                    }
                }
        }                  
                         
		else if (typ == 1)  //sphere
         {
           	for (int i1 = 0; i1 < num1; i1++)
              for (int i2 = 0; i2 < num2; i2++)
                for (int i3 = 0; i3 < num3; i3++) 
				  
                  
                  if (((i1+zbeg-1-p1)*(i1+zbeg-1-p1) + (i2-p2)*(i2-p2) + (i3-start3)*(i3-start3)) < rad*rad)
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
       else if (typ == 2) //half space
          {
	      	for (int i1 = 0; i1 < num1; i1++)
              for (int i2 = 0; i2 < num2; i2++)
                   for (int i3 =0; i3 <= num3; i3++)  
                   { 
                       if (i3<=round(p1*(i1+zbeg))){
                                                    cout<<"p1 is "<<p1<<"\n";
                   cout<<"half space is "<<round(p1*(i1+zbeg))<<"\n";         
				      B.set(i1,i2,i3,2);}
                  }
                     
	      }
	   else if (typ == 3) // 3D rectanglar void
          {
	      	for (int i1 = 0; i1 < num1; i1++)
              for (int i2 = 0; i2 < num2; i2++)
                for (int i3 = 0; i3 < num3; i3++)  
				  if ((i1+zbeg-1 >= p1) && (i1+zbeg-1 <= p2) && (i2 >= start3) && (i2 <= end3) && (i3 >= rad) && (i3 <= dd))
				  {
                     B.set(i1,i2,i3,2); // Stress-free void
				  }



	      }
	      
	        else if (typ == 101) //Rounded rectangle, experimental flaw
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
                for (int i3 = 0; i3 < num3; i3++) {
				                    
                  if ((i1+zbeg-1)>(p1+rad) && (i1+zbeg-1)<(endx-rad) && i2>starty && i2<endy && i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }
                 if ((i1+zbeg-1)>(p1) && (i1+zbeg-1)<(p1+rad) && i2>(starty+rad) && i2<(endy-rad) && i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }
                      
                   if ((i1+zbeg-1)>(endx-rad) && (i1+zbeg-1)<(endx) && i2>(starty+rad) && i2<(endy-rad) && i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }
                  if (((i1+zbeg-1-(endx-rad))*(i1+zbeg-1-(endx-rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }   
                      
                      if (((i1+zbeg-1-(startx+rad))*(i1+zbeg-1-(startx+rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }   
                      
                      if (((i1+zbeg-1-(endx-rad))*(i1+zbeg-1-(endx-rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }   
                      
                      if (((i1+zbeg-1-(startx+rad))*(i1+zbeg-1-(startx+rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num3-dd))
					{
					  B.set(i1,i2,i3,2);//rigid scatterer, or create boundary
                      }    
                              
                      }
                     
	      }

}

int round(double a){
    return int(a+0.5);}
};
