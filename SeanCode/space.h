/*
 * space.h
 * Sets up the transducers and reflectors for the Cartesian EFIT simulation
 * Redesigned integration of 'setup_cart_space.h' and 'array3D.h' (Dieckman)
 *
 *  Created on: March 27, 2019
 *      Author: Sean M. Raley (UNH)
 */

#ifndef SPACE_H_
#define SPACE_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include "transducer.h"

using namespace std;

// TODO: Are these defs necessary?
//#define min(a,b) (((a)<(b))?(a):(b))
//#define max(a,b) (((a)>(b))?(a):(b))

class space{

public:
	space(double *params);
	~space();

	int num_z, num_y, num_x; // number of grid points in each direction
	int num_zB, num_yB, num_xB; // number of Boundary grid points in each direction

	int abc;    // number of abc points on each end

	double ds;	// spatial step size (meters)
	double dt;	// time step size (seconds)

//	int abc;    // number of abc points on each end TODO: delete?

	double den;	// default density (kg/m^3)
	double lm;	// default Lame constant - lambda
	double mu;	// default Lame constant - mu

	int zbeg;	// z start position (where divvacuum ends)
	int type;	// type -> 1 = left, 2 = middle, 3 = right (for MPI)

	double *vz;	// velocities in z-dir
	double *vy;	// velocities in y-dir
	double *vx;	// velocities in x-dir
	double *T11;// zz normal stress
	double *T22;// yy normal stress
	double *T33;// xx normal stress
	double *T12;// zy normal stress
	double *T23;// yx normal stress
	double *T13;// zy normal stress

	double *d;	// density
	double *lmd;// Lame parameter - lambda
	double *muu;// Lame parameter - mu

	int *B;		// Boundary array

	int drivetime;

	transducer *trans{nullptr};
	int numtrans;

	int Lyx;	// array length of entire x-y plane
	int Lzyx;	// array length of entire x-y-z volume
	int LyxB;	// array length of entire x-y plane (Boundary array)
	int LzyxB;	// array length of entire x-y-z volume (Boundary array)

	double dtods;	// time step over spatial step

private:
	double lmdtods;
	double l2mdtods;
	double mdtods;

	int iz, iy, ix;	// counters for spatial loops
	int pp1,pm1;

	double PIo2;

	//========================================================
	// low-level initialization function
	// TODO: (low priority) may want to overload to initialize without a value for clear()
	template <class SomeType>
	SomeType *init(SomeType def){
		SomeType *temparray = new SomeType[Lzyx];
		clear(temparray,def);
		return temparray;
	}

	// special version for B array
	int *initB(int def){
		int *temparray = new int[LzyxB];
		clearB(temparray, def);
		return temparray;
	}

public:
	// ===================================================================================
	// Update velocities and stresses
	// ===================================================================================
	void UpdateVs(int zs, int zend){
		int tau_time;
		int tau_count = 0;

		for (iz = zs; iz <= zend; iz++){
			setindx(iz,0,0);
			for (iy = 0; iy < num_y; iy++){
				setindxB(iz+1,iy+1,1);
				for (ix = 0; ix < num_x; ix++){

		  // --- vy ---
					if (vB(B)==0){
						sv(vy, v(vy) + 2*dtods/(v(d)+vyp(d))*((v(T12)-vzm(T12))+(vyp(T22)-v(T22))+(v(T23)-vxm(T23))));
					}
					else if (vB(B)==2 || vzmB(B)==2 || vxmB(B)==2) {} // because vy requires elements in direction vxm() and vzm()
		  //vr.sv( vr.v() + dtodsp*(trans[vB(B)-1000].drivef(drivetime)+(1/rr)*((1/(2*(rr)))*(3*Tpp.v()-Tpp.vym()))) );
					else if (vypB(B)==2){
						sv(vy, v(vy) + 2*dtods/(v(d)+v(d))*(-2*(v(T22))));
					}
					else if (vymB(B)==2){
						sv(vy, v(vy) + 2*dtods/(v(d)+v(d))*(2*(vyp(T22))));
					}
					else{
						sv(vy, v(vy) + 2*dtods/(v(d)+vyp(d))*((v(T12)-vzm(T12))+(vyp(T22)-v(T22))+(v(T23)-vxm(T23))));
					}


		  // --- vz ---
					if (vB(B)==0){
						sv(vz, v(vz) + 2*dtods/(v(d)+vzp(d))*((vzp(T11)-v(T11))+(v(T12)-vym(T12))+(v(T13)-vxm(T13))));  // middle
					}
					else if (vB(B)==2 || vymB(B)==2 || vxmB(B)==2) {}  //b/c vz needs elements in the vym() and vxm() directions
					else if (vzpB(B)==2){
						sv(vz, v(vz) + 2*dtods/(v(d)+v(d))*(-2*(v(T11)))); //right end (only called if at last node)
					}
					else if (vzmB(B)==2){
						sv(vz, v(vz) + 2*dtods/(v(d)+v(d))*(2*(vzp(T11)))); //left end (only called if at first node)
					}
					else{
						sv(vz, v(vz) + 2*dtods/(v(d)+vzp(d))*((vzp(T11)-v(T11))+(v(T12)-vym(T12))+(v(T13)-vxm(T13))));
					}


		  // --- vx ---
					if (vB(B)==0){
						sv(vx, v(vx) + 2*dtods/(v(d)+vxp(d))*((v(T13)-vzm(T13))+(v(T23)-vym(T23))+(vxp(T33)-v(T33))) );
					}
					else if (vB(B)==2 || vymB(B)==2 || vzmB(B)==2) {}
					else if ( vB(B)>=9000){
						sv(vx, v(vx) + 2*dtods/(v(d)+v(d))*(-2*(v(T33)))+ 2*dtods/(den+den)*(-trans[vB(B)-9000].drivef(drivetime)));
					}
//					else if ( vB(B)>=1000){ // original definition
//						sv(vx, v(vx) + 2*dtods/(v(d)+v(d))*(-2*(v(T33)))+ 2*dtods/(den+den)*(trans[vB(B)-1000].drivef(drivetime)));
//					}
					else if ( vB(B)>=1000){
//						*currentloc = giveindx();
						giveindx();
						tau_count = currentloc[0]*num_y+currentloc[1];
//			    		if (((currentloc[0]+zbeg-trans[vB(B)-1000].posi1)*(currentloc[0]+zbeg-trans[vB(B)-1000].posi1)+(currentloc[1]-trans[vB(B)-1000].posi2)*(currentloc[1]-trans[vB(B)-1000].posi2)) <= (trans[vB(B)-1000].radius*trans[vB(B)-1000].radius)){
							if(trans[vB(B)-1000].tau[tau_count]<drivetime){
								tau_time = drivetime - trans[vB(B)-1000].tau[tau_count];
							}
							else{
								tau_time=0;
							}
//							tau_count++;
//			    		}
							// debug lines
//							cout<<"@tau_count="<<tau_count<<", ci="<<ci<<", (y,z)=("<<currentloc[1]<<","<<currentloc[0]<<"), (iy,iz)=("<<iy<<","<<iz<<"), drivetime="<<drivetime<<", tau_time="<<tau_time<<endl;
//							cout<<"@(y,z)=("<<currentloc[1]<<","<<currentloc[0]<<"), tau_time="<<tau_time<<endl;
						sv(vx, v(vx) + 2*dtods/(v(d)+v(d))*(-2*(v(T33)))+ 2*dtods/(den+den)*(trans[vB(B)-1000].drivef(tau_time)));
					}
					else if (vxpB(B)==2){
						sv(vx,  v(vx) + 2*dtods/(v(d)+v(d))*(-2*(v(T33)))); //top
					}
					else if (vxmB(B)==2){
						sv(vx, v(vx) + 2*dtods/(v(d)+v(d))*(2*(vxp(T33)))); //bottom
					}
					else{
						sv(vx, v(vx) + 2*dtods/(v(d)+vxp(d))*((v(T13)-vzm(T13))+(v(T23)-vym(T23))+(vxp(T33)-v(T33))) );
					}

					incindx();
					incindxB();
				}
			}
		}
	}

	void UpdateTs(int zs, int zend){
		for (iz = zs; iz <= zend; iz++){
			setindx(iz,0,0);
			for (iy = 0; iy < num_y; iy++){
				setindxB(iz+1,iy+1,1);
				for (ix = 0; ix < num_x; ix++){

		  // --- Tii ---
					if (vB(B)==2 || vymB(B)==2 || vxmB(B)==2 || vzmB(B)==2) {} // 1 is iz, 2 is iy, 3 is ix
					else{
						sv(T11, v(T11)+dtods*(v(lmd)+(2*v(muu)))*(v(vz)-vzm(vz))+dtods*(v(lmd))*((v(vy)-vym(vy))+(v(vx)-vxm(vx))));
						sv(T22, v(T22)+dtods*(v(lmd)+(2*v(muu)))*(v(vy)-vym(vy))+dtods*(v(lmd))*((v(vz)-vzm(vz))+(v(vx)-vxm(vx))));
						sv(T33, v(T33)+dtods*(v(lmd)+(2*v(muu)))*(v(vx)-vxm(vx))+dtods*(v(lmd))*((v(vz)-vzm(vz))+(v(vy)-vym(vy))));
					}

		  // --- T13 ---
					if (vB(B)==2 || vzpB(B)==2 || vzmB(B)==2 || vxpB(B)==2 || vxmB(B)==2 || vymB(B)==2 || vzxpB(B)==2) {} // need all terms needed by vz and vx as well as any additional terms
					else{
						sv(T13, v(T13)+dtods*(4/(((1/v(muu))+(1/vzp(muu))+(1/vxp(muu))+(1/vzxp(muu)))))*((vxp(vz)-v(vz))+(vzp(vx)-v(vx))));
					}

		  // --- T23 ---
					if (vB(B)==2 || vypB(B)==2 || vymB(B)==2 || vxpB(B)==2 || vxmB(B)==2 || vzmB(B)==2 || vyxpB(B)==2) {}
					else{
						sv(T23, v(T23)+dtods*(4/((1/(v(muu))+(1/vyp(muu))+(1/vxp(muu))+(1/vyxp(muu)))))*((vxp(vy)-v(vy))+(vyp(vx)-v(vx))));
					}

		  // --- T12 ---
					if (vB(B)==2 || vypB(B)==2 || vymB(B)==2 || vzpB(B)==2 || vzmB(B)==2 || vxmB(B)==2 || vzypB(B)==2) {}
					else{
						sv(T12, v(T12)+dtods*((4/((1/v(muu))+(1/vzp(muu))+(1/vyp(muu))+(1/vzyp(muu)))))*((vyp(vz)-v(vz))+(vzp(vy)-v(vy))));
					}

					incindx();
					incindxB();
				}
			}
		}
	}

	// ===================================================================================
	// Add/update transducers (both original and IDT methods below)
	// ===================================================================================
	void UpdateTransducers(int t){
		int tr;
	    for (int i1 = 1; i1<num_z-1; i1++){
	    	for (int i2 = 0; i2<num_y-1; i2++){
				if (valB(B,i1+1,i2+1,num_x) >= 1000){
					tr = valB(B,i1+1,i2+1,num_x)-1000; //transducer ID
					trans[tr].record[t] =  trans[tr].record[t] + val(vx,i1,i2,num_x-1); //record holds recorded info
				}
	    	}
	    }
	}


	void addTransducer(transducer t){

	    // --- Original way, with a circular transducer ---
	//    /*
		numtrans = numtrans+1;
	    transducer *temp = new transducer[numtrans]; // of class transducer, length numtrans

	    for (int i = 0; i<numtrans-1; i++){ // loops from 0 to actual number of transducers
	    	temp[i] = trans[i];
	    }

	    temp[numtrans-1] = t;
	    trans = temp;
	    int nelems = 0;
		//debug line
//		cout<<"t.drive[1] = "<<t.drive[1]<<endl;
		//debug line
//		cout<<"temp[numtrans-1].drive[1] = "<<temp[numtrans-1].drive[1]<<endl;
		//debug line
//		cout<<"trans[numtrans-1].drive[1] = "<<trans[numtrans-1].drive[1]<<endl;
//		for (int i1 = 1; i1<num_z-1; i1++){ // working version
//	    	for (int i2 = 0; i2<num_y; i2++){
//	    		if (((i1+zbeg-trans[numtrans-1].posi1)*(i1+zbeg-trans[numtrans-1].posi1)+(i2-trans[numtrans-1].posi2)*(i2-trans[numtrans-1].posi2)) <= (trans[numtrans-1].radius*trans[numtrans-1].radius)){
//	    			setB(B,i1+1,i2+1,t.posi3+1,1000+numtrans-1); //executes if (z-zo)^2+(y-yo)^2<=rad^2, ie, if within a circle of transducer radius then set B >1000
//	    			nelems++;
//	    		}
//	    	}
//		}
	    int offset_marker, xoffset;
		for (int i1 = 1; i1<num_z-1; i1++){ // "x-drop" version for xdcr surfaces not planar on max x
	    	for (int i2 = 0; i2<num_y; i2++){
	    		if (((i1+zbeg-trans[numtrans-1].posi1)*(i1+zbeg-trans[numtrans-1].posi1)+(i2-trans[numtrans-1].posi2)*(i2-trans[numtrans-1].posi2)) <= (trans[numtrans-1].radius*trans[numtrans-1].radius)){
	    			offset_marker=0;
	    			xoffset=0;
	    			while(offset_marker==0){
						if(valB(B,i1+1,i2+1,t.posi3+1-xoffset)==0){
							setB(B,i1+1,i2+1,t.posi3+1-xoffset,1000+numtrans-1); //executes if (z-zo)^2+(y-yo)^2<=rad^2, ie, if within a circle of transducer radius then set B >1000
							nelems++;
							offset_marker++;
						}
						else if(t.posi3-xoffset<0){
							cout<<"WARNING! Transducer x-drop failed at (y,z) = ("<<i2<<","<<(i1+zbeg)<<")"<<endl;
							offset_marker++;
						}
						else{
							xoffset++;
						}
	    			}
	    		}
	    	}
		}

		if(nelems != trans[numtrans-1].tau_values){
			cout<<"WARNING! For worker starting at z="<<zbeg<<"Number of transducer elements defined by addTransducer() does NOT equal number of transducer elements in delay array tau. Transducer errors likely."<<endl;
		}

		//if (t.posi1-zbeg+1 >0) B.set(t.posi1-zbeg+1,t.posi2+1,t.posi3+1,1000+numtrans-1); //TODO: delete?
		temp[numtrans-1].numelems = nelems; //TODO: delete?
//		delete[] temp; // TODO: why does adding this line cause ALL simulation data to be 0?
	//    */

// NOTE: THIS SECTION NOT CONVERTED FROM OLD EFIT_cart
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

	    for (int i1 = 1; i1<num_z-1; i1++)
	    {
	      // cout << "i1 = " << i1 << " and i1+zbeg-1 = " << i1+zbeg-1 << endl;  // Debugging only
	      for (int i2 = 0; i2<num_y-1; i2++)
		for (int i = -5; i <= 4; i++) // Ten finger pairs, centered at transducer_x
		{
		  pos_left  = trans[numtrans-1].posi1 + wlength*(i);
	          pos_right = pos_left + wlength/4;
	          neg_left  = pos_right + wlength/4;
	          neg_right = neg_left + wlength/4;

		  if ((i1+zbeg >= pos_left) && (i1+zbeg <= pos_right) && (i2 >= bot_y) && (i2 <= top_y))
		  {
		    B.set(i1+1,i2+1,t.posi3+1,1000+numtrans-1);//sets B value to 1000 or higher for x and y and z positions of positive transducer
		    nnodes++;
	          }

		  else if ((i1+zbeg >= neg_left) && (i1+zbeg <= neg_right) && (i2 >= bot_y) && (i2 <= top_y))
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
	void addReflector(double typ, double p1, double p2, int start3, int end3, double rad, double dd, double mu, double lambda){
	    // --- VERTICAL 2D CRACK ---
		if (typ == 0){
	      // Below is for easy of use only, no real purpose in defining new variables
			int z1 = static_cast<int>(p2);
			int y1 = start3;
			int z2 = end3;
			int y2 = static_cast<int>(rad);
			int depth  = static_cast<int>(dd); // Crack depth in x direction

			cout<<"Crack set from p1 = ["<<z1<<" "<<y1<<"] to p2 = ["<<z2<<" "<<y2<<"] at a depth of x = "<<depth<<endl;

			double slope;
			int y_int;
			int z;
			int y;
			int begin;
			int finish;
			int temp;

			if (z2 == z1){
				slope = 9999; // approximates infinity
			}
			else{
				slope = static_cast<double>((y2-y1)/(z2-z1));
			}

			y_int = y1 - slope*z1; // TODO: add static_cast<int> ?
			cout<<" Slope = "<<slope<<" and y_int = "<<y_int<<endl;

			if ((z1-z2)*(z1-z2) > (y1-y2)*(y1-y2)){ // Fixes stretching
				temp = 1;
				if (z2 > z1){
					begin = z1;
					finish = z2;
				}
				else{
					begin = z2;
					finish = z1;
				}

				// cout<<"Restricting to z begin = "<<begin<<" to finish = "<<finish<<endl;
			}
			else{
				temp = 2;
				if (y2 > y1){
					begin = y1;
					finish = y2;
				}
				else{
					begin = y2;
					finish = y1;
				}

				// cout<<"Restricting to y begin = "<<begin<<" to finish = "<<finish<<endl;
			}

			for (int i1 = 0; i1 < num_z; i1++){
				for (int i2 = 0; i2 < num_y; i2++){
					for (int i3 = 0; i3 < depth; i3++){
						if (temp == 1){
							y = round(slope*i1+zbeg-1 + y_int);
							if ((i2 == y) && (i1+zbeg-1 >= begin) && (i1+zbeg-1 <= finish)){
								// cout<<"Setting crack at ["<<i1+zbeg-1<<" "<<i2<<" "<<i3<<"]!!"<<endl;
								setB(B,i1+1,i2+1,i3+1,2);
							}
						}

						if (temp == 2){
							z = round((i2-y_int)/slope);
							if ((i1+zbeg-1 == z) && (i2 >= begin) && (i2 <= finish)){
								// cout<<"Setting crack at ["<<i1+zbeg-1<<" "<<i2<<" "<<i3<<"]!!"<<endl;
								setB(B,i1+1,i2+1,i3+1,2);
							}
						}
					}
				}
			}
		}

	    // --- SPHERE ---
	    else if (typ == 1){
	    	for (int i1 = 0; i1 < num_z; i1++){
	    		for (int i2 = 0; i2 < num_y; i2++){
	    			for (int i3 = 0; i3 < num_x; i3++){
	    				if (((i1+zbeg-1-p1)*(i1+zbeg-1-p1) + (i2-p2)*(i2-p2) + (i3-start3)*(i3-start3)) < rad*rad){
	    					if ((dd == -1)){
	    						setB(B,i1+1,i2+1,i3+1,2);//rigid scatterer
	    					}
	    					else{
	    						set(muu,i1,i2,i3,mu);
	    						set(d,i1,i2,i3,dd);
	    						set(lmd,i1,i2,i3,lambda);
	    					}
	    				}
	    			}
	    		}
	    	}
		}

	    // --- HALF SPACE ---
	    else if (typ == 2){
	    	for (int i1 = 0; i1 < num_z; i1++){
	    		for (int i2 = 0; i2 < num_y; i2++){
	    			for (int i3 =0; i3 <= num_x; i3++){
	    				if (i3<=round(p1*(i1+zbeg))){
	    					cout<<"p1 is "<<p1<<"\n";
	    					cout<<"half space is "<<round(p1*(i1+zbeg))<<"\n";
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    			}
	    		}
	    	}
		}

	    // --- 3D RECTANGULAR VOID ---
	    else if (typ == 3){
	    	for (int i1 = 0; i1 < num_z; i1++){
	    		for (int i2 = 0; i2 < num_y; i2++){
	    			for (int i3 = 0; i3 < num_x; i3++){
	    				if ((i1+zbeg-1 >= p1) && (i1+zbeg-1 <= p2) && (i2 >= start3) && (i2 <= end3) && (i3 >= rad) && (i3 <= dd)){
	    					setB(B,i1+1,i2+1,i3+1,2); // Stress-free void
	    				}
	    			}
	    		}
	    	}
		}

		// --- arb 3d scatterer from STL file ---
		else if (typ == 5)
		{
			char inputFilename[] = "arbscatt.file";
			ifstream inFile;
			inFile.open("arbscatt.file", ios::in); //arbscatt file follows (x,y,z) orientation

			if (!inFile){
       			cerr << "Can't open input file " << inputFilename << endl;
       			exit(1);
			}
			//debug line
//			cout<<"Arbscat debug 1"<<endl;

			int tss = rad; // total size of space
			int *scatterspace = new int[tss];
			for (int i = 0; i<tss; i++){
    			inFile >> scatterspace[i]; // read in entire space including scatterer
    		}

    		// Nested loops through width, height, depth to compute array indices for linear array
    		// here, num1 is z, num2 is y, and num3 is x -> (y,x,z) == (num2, num3, num1)
	        //for (i2 = 0; i2 < num2; i2++) // width (num2/i2)
            //  for (i3 = 0; i3 < num3; i3++) // height (num3/i3)
            //   	for (i1 = 0; i1< num1; i1++) { //depth (num1/i1)
			//		if (scatterspace[(long)i2*(long)num3*(long)num1 + (long)i3*(long)num1 + (long)i1] == 1)

			int tempB;
	        for (iz = 0; iz < num_z; iz++){
	        	//debug line
//	        	if(type==1){cout<<"Arbscat debug 2 @iz="<< iz <<endl;}
	        	for (iy = 0; iy < num_y; iy++){
	        		for (ix = 0; ix< num_x; ix++){
	        			tempB = valB(B,iz+1,iy+1,ix+1);
	        			if ((scatterspace[(iz+zbeg)*num_y*num_x + iy*num_x + ix] != 1) && (tempB != 1) && (tempB != 2)){
	        				if ((lambda == -1) && (mu == -1) && (dd == -1)){
	        					setB(B,iz+1,iy+1,ix+1,2);
	        				}
	        				else{
	        					set(lmd,iz,iy,ix,lambda);
	        					set(muu,iz,iy,ix,mu);
	        					set(d,iz,iy,ix,dd);
	        				}
	        			}
	        			//i++;
	        		}
	        	}
	        }
			inFile.close();
			delete[] scatterspace;
		}

		// --- ROUNDED RECTANGLE (EXPERIMENTAL)
	    else if (typ == 101){
	    	//double p1, double p2, int start3, int end3, double rad, double dd, double mu, double lambda)
	    	//p1=startx
	    	//p2=lenx
	    	//start3=starty
	    	//end3=leny
	    	int startx=static_cast<int>(p1);
	    	int starty=start3;
	    	int endx=static_cast<int>(p1+p2);
	    	int endy=start3+end3;
	    	//rad=radius of spheres making up edges
	    	//dd=thickness

	    	for (int i1 = 0; i1 < num_z; i1++){
	    		for (int i2 = 0; i2 < num_y; i2++){
	    			for (int i3 = 0; i3 < num_x; i3++){
	    				if ((i1+zbeg-1)>(p1+rad) && (i1+zbeg-1)<(endx-rad) && i2>starty && i2<endy && i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2); // rigid scatterer, or create boundary
	    				}
	    				if ((i1+zbeg-1)>(p1) && (i1+zbeg-1)<(p1+rad) && i2>(starty+rad) && i2<(endy-rad) && i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    				if ((i1+zbeg-1)>(endx-rad) && (i1+zbeg-1)<(endx) && i2>(starty+rad) && i2<(endy-rad) && i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    				if (((i1+zbeg-1-(endx-rad))*(i1+zbeg-1-(endx-rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    				if (((i1+zbeg-1-(startx+rad))*(i1+zbeg-1-(startx+rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    				if (((i1+zbeg-1-(endx-rad))*(i1+zbeg-1-(endx-rad)) + (i2-(starty+rad))*(i2-(starty+rad)) ) < rad*rad  &&  i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    				if (((i1+zbeg-1-(startx+rad))*(i1+zbeg-1-(startx+rad)) + (i2-(endy-rad))*(i2-(endy-rad)) ) < rad*rad  &&  i3>(num_x-dd)){
	    					setB(B,i1+1,i2+1,i3+1,2);
	    				}
	    			}
	    		}
	    	}
		}

	}


	int round(double a){
		return int(a+0.5);
	}

	//========================================================
	// Array Manipulation

private:
	int ci; // current index
	int ciB; // current Boundary index (used for Boundary array only)
	int currentloc[3]; // for storing current zyx location

public:
	// Returns value at i_z, i_y, i_x
	template <class SomeType>
	SomeType val(SomeType *a, int i_z, int i_y, int i_x){
		return a[(i_z*Lyx)+(i_y*num_x)+i_x];
	}

	// Sets value at i_z, i_y, i_x
	template <class SomeType>
	void set(SomeType *a, int i_z, int i_y, int i_x, SomeType val){
		a[(i_z*Lyx)+(i_y*num_x)+i_x] = val;
	}

	// Returns Boundary Array value at i_z, i_y, i_x
	template <class SomeType>
	SomeType valB(SomeType *a, int i_z, int i_y, int i_x){
		return a[(i_z*LyxB)+(i_y*num_xB)+i_x];
	}

	// Sets Boundary Array value at i_z, i_y, i_x
	template <class SomeType>
	void setB(SomeType *a, int i_z, int i_y, int i_x, SomeType val){
		a[(i_z*LyxB)+(i_y*num_xB)+i_x] = val;
	}
	//========================================================
	// Quick Access Methods

	// sets index ci
	void setindx(int iz, int iy, int ix){
		ci = (iz*Lyx)+(iy*num_x)+ix;
	}
	// increments index counter by 1
	void incindx(){
		ci = ci+1;
	}
	void sv(double *a, double x){a[ci]=x;				}	// sets value at ci

	template <class SomeType>
	SomeType v(SomeType *a)		{return a[ci];			}	// equiv of a[iz][iy][ix]
	template <class SomeType>
	SomeType vzp(SomeType *a)	{return a[ci+Lyx];		}	// equiv of a[iz+1][iy][ix]
	template <class SomeType>
	SomeType vzm(SomeType *a)	{return a[ci-Lyx];		}	// equiv of a[iz-1][iy][ix]
	template <class SomeType>
	SomeType vyp(SomeType *a)	{return a[ci+num_x];	}	// equiv of a[iz][iy+1][ix]
	template <class SomeType>
	SomeType vyp2(SomeType *a)	{return a[ci+2*num_x];	}	// equiv of a[iz][iy+2][ix]
	template <class SomeType>
	SomeType vym(SomeType *a)	{return a[ci-num_x];	}	// equiv of a[iz][iy-1][ix]
	template <class SomeType>
	SomeType vxp(SomeType *a)	{return a[ci+1]; 		}	// equiv of a[iz][iy][ix+1]
	template <class SomeType>
	SomeType vxm(SomeType *a)	{return a[ci-1];		}	// equiv of a[iz][iy][ix-1]
	template <class SomeType>
	SomeType vzxp(SomeType *a)	{return a[ci+1+Lyx];	}	// equiv of a[iz+1][iy][ix+1]
	template <class SomeType>
	SomeType vyxp(SomeType *a)	{return a[ci+1+num_x];	}	// equiv of a[iz][iy+1][ix+1]
	template <class SomeType>
	SomeType vzyp(SomeType *a)	{return a[ci+Lyx+num_x];}	// equiv of a[iz+1][iy+1][ix]

	// sets array for z,y,x at current index ci
//	int* giveindx(){
	void giveindx(){
//		int *indxarray = new int[3];
		currentloc[0] = ((ci / num_x) / num_y) % num_z; // z position
		currentloc[1] = (ci / num_x) % num_y; // y position
		currentloc[2] = ci % num_x; // x position
//		return indxarray;
	}

	//========================================================
	// Quick Access Methods (Boundary array)

	// sets index ciB
	void setindxB(int iz, int iy, int ix){
		ciB= (iz*LyxB)+(iy*num_xB)+ix;
	}

	// increments Boundary index counter by 1
	void incindxB(){
		ciB = ciB+1;
	}
	void svB(double *a, double x){a[ciB]=x;				}	// sets value at ciB

	template <class SomeType>
	SomeType vB(SomeType *a)	{return a[ciB];			}	// equiv of a[iz][iy][ix]
	template <class SomeType>
	SomeType vzpB(SomeType *a)	{return a[ciB+LyxB];		}	// equiv of a[iz+1][iy][ix]
	template <class SomeType>
	SomeType vzmB(SomeType *a)	{return a[ciB-LyxB];		}	// equiv of a[iz-1][iy][ix]
	template <class SomeType>
	SomeType vypB(SomeType *a)	{return a[ciB+num_xB];	}	// equiv of a[iz][iy+1][ix]
	template <class SomeType>
	SomeType vyp2B(SomeType *a)	{return a[ciB+2*num_xB];	}	// equiv of a[iz][iy+2][ix]
	template <class SomeType>
	SomeType vymB(SomeType *a)	{return a[ciB-num_xB];	}	// equiv of a[iz][iy-1][ix]
	template <class SomeType>
	SomeType vxpB(SomeType *a)	{return a[ciB+1]; 		}	// equiv of a[iz][iy][ix+1]
	template <class SomeType>
	SomeType vxmB(SomeType *a)	{return a[ciB-1];		}	// equiv of a[iz][iy][ix-1]
	template <class SomeType>
	SomeType vzxpB(SomeType *a)	{return a[ciB+1+LyxB];	}	// equiv of a[iz+1][iy][ix+1]
	template <class SomeType>
	SomeType vyxpB(SomeType *a)	{return a[ciB+1+num_xB];	}	// equiv of a[iz][iy+1][ix+1]
	template <class SomeType>
	SomeType vzypB(SomeType *a)	{return a[ciB+LyxB+num_xB];}	// equiv of a[iz+1][iy+1][ix]

	//========================================================
	// sets all values = def
	template <class SomeType>
	void clear(SomeType *a,SomeType def){
		for(int i=0; i<Lzyx; i++){
			a[i]=def;
		}
	}

	// special version for Boundary array
	void clearB(int *a, int def){
		for(int i=0; i<LzyxB; i++){
			a[i]=def;
		}
	}

	// sets all outer boundary surfaces for Boundary array
	void setBoundaries(){
	    for (int i1 = 0; i1< num_zB; i1++){
	    	for (int i3 = 0; i3< num_xB; i3++){
	    		setB(B,i1,0,i3,2);  //set x2 boundary
	    		setB(B,i1,1,i3,1);
	    		//setB(B,i1,2,i3,1);
	    		setB(B,i1,num_yB-1,i3,2);
	    		setB(B,i1,num_yB-2,i3,2);
	    		setB(B,i1,num_yB-3,i3,1);
	    	}
	    }

	    for (int i1 = 0; i1< num_zB; i1++){
	    	for (int i2 = 0; i2< num_yB; i2++){
	    		setB(B,i1,i2,0,2);
	    		setB(B,i1,i2,1,1);
	    		setB(B,i1,i2,num_xB-1,2);//set x3 boundary
	    		setB(B,i1,i2,num_xB-2,2);

	    		if ((i2>0) && (i2<num_zB-2)){
	    			setB(B,i1,i2,num_xB-3,1);
	    			// setB(B,i1,i2,2,1);
	    		}
	    	}
	    }

	    for (int i2 = 0; i2< num_yB; i2++){
	    	for (int i3 = 0; i3< num_xB; i3++){
	    		if (type == 1){
	    			setB(B,0,i2,i3,2);  // if at actual end of space in z set boundary
	    			//setB(B,1,i2,i3,2);
	    			if ( (i2 > 0) && (i2<(num_yB-2)) && (i3 > 0) && (i3<(num_xB-2))){
	    				setB(B,1,i2,i3,1);
	    			}
	    		}

	    		if (type == 3){
	    			setB(B,num_zB-1,i2,i3,2);
	    			setB(B,num_zB-2,i2,i3,2);
	    			if ( (i2 > 0) && (i2<(num_yB-2)) && (i3 >0) && (i3<(num_xB-2))){
	    				setB(B,num_zB-3,i2,i3,1);
	    			}
	    		}
	    	}
	    }

	}

	//========================================================
	// returns 2D slice through 3D array at fixed index y
	double* slice_fixy(double *a, int y_slice){
		double *slice = new double[(num_z-2)*num_x];
		//slice[0] = (num_z-2)*num_x; // ACTION ITEM: figure out what this line does
		int count = 0;
		for(iz=1; iz<num_z-1; iz++) // does not return ends
			for(ix=0; ix<num_x; ix++){
				slice[count] = val(a,iz,iy,ix);
				count++;
			}
		return slice;
	}

	//========================================================
	// returns 3D volume returning only the even indexes
	double* GetEvenVol(double *a, int len){ // removed "start" input
//		int len = GetEvenVolLen();
		double *EvenArray = new double[len];
		int count = 0;
		for(iz=1+(num_z%2); iz<num_z-1; iz+=2) // does not return ends
			for(iy=0; iy<num_y-1; iy+=2)
				for(ix=0; ix<num_x-1; ix+=2){
//					EvenArray[count] = (val(a, iz-1, iy, ix)+val(a, iz, iy, ix))/2;
					EvenArray[count] = val(a, iz, iy, ix);
					//debug line
//					if(type==1) cout<<count<<" ";

//					EvenArray[count] = val(a, iz-2, iy, ix); // to check drive plane
					count++;
				}
		// debug line
//		cout<<"End of GetEvenVol"<<endl;
		return EvenArray;
	}

	// overload to retrieve magnitude of full velocity vector
	double* GetEvenVol(int len){ // removed "start" input
//		int len = GetEvenVolLen();
		double *EvenArray = new double[len];
		int count = 0;
		double xtemp, ytemp, ztemp;
		for(iz=1+(num_z%2); iz<num_z-1; iz+=2) // does not return ends
			for(iy=0; iy<num_y-1; iy+=2)
				for(ix=0; ix<num_x-1; ix+=2){
//					EvenArray[count] = (val(a, iz-1, iy, ix)+val(a, iz, iy, ix))/2;
					xtemp = val(vx, iz, iy, ix);
					ytemp = val(vy, iz, iy, ix);
					ztemp = val(vz, iz, iy, ix);
					EvenArray[count] = sqrt(xtemp*xtemp+ytemp*ytemp+ztemp*ztemp);
					//debug line
//					if((ix*iy*iz)%1000000==0) cout<<"(x,y,z) = "<<xtemp<<ytemp<<ztemp <<endl;

					count++;
				}
		return EvenArray;
	}

	int GetEvenVolLen(){ // removed "start" input
		int len;
		if(num_z%2==0)
			len = (num_z-1)/2*(num_y/2)*(num_x/2); // len if num_z is even
		else
			len = (num_z-2)/2*(num_y/2)*(num_x/2); // len if num_z is odd
// TODO: (num_z-2) and (num_z-3), or (num_z-1) and (num_z-2)? Probably shouldn't matter, actually.

		// for debugging incorrect length received from node 3
//		cout<<"type "<<type<<" len returned from GetEvenVolLen is: "<<len<<endl;
//		cout<<"type "<<type<<" num_z from GetEvenVolLen is: "<<num_z<<endl;
//		cout<<"type "<<type<<" num_y from GetEvenVolLen is: "<<num_y<<endl;
//		cout<<"type "<<type<<" num_x from GetEvenVolLen is: "<<num_x<<endl;

		return len;
	}
};

space::space(double *params){
	num_z = params[0]+2; 	// number of nodes in z-direction
	num_y = params[1]; 		// number of nodes in y-direction
	num_x = params[2]; 		// number of nodes in x-direction
	ds = params[3]; 		// spatial step size (m)
	dt = params[4]; 		// time step size (s)

	den = params[5]; 		// density
	lm = params[6]; 		// Lame constant - lambda
	mu = params[7]; 		// Lame constant - mu
	zbeg = params[8]; 		// simspace z-starting position for each node

	Lyx = num_y*num_x;
	Lzyx = num_z*num_y*num_x;

	vz = init(0.0);
	vy = init(0.0);
	vx = init(0.0);
	T11 = init(0.0);
	T22 = init(0.0);
	T33 = init(0.0);
	T12 = init(0.0);
	T23 = init(0.0);
	T13 = init(0.0);

	d = init(den);
	lmd = init(lm);
	muu = init(mu);

	num_zB = num_z+2; 	// number of Boundary nodes in z-direction
	num_yB = num_y+2; 	// number of Boundary nodes in y-direction
	num_xB = num_x+2; 	// number of Boundary nodes in x-direction
	LyxB = num_yB*num_xB;
	LzyxB = num_zB*num_yB*num_xB;

	B = initB(0);
	setBoundaries();

	dtods = dt/ds;
    lmdtods  = (lm*dt)/ds;
    l2mdtods = ((lm+2*mu)*dt)/ds;
    mdtods   = (mu*dt)/ds;

    PIo2 = 3.14159265358979/2;

    numtrans=0;
	drivetime = 0;
    abc = 80;
}

space::~space() {
	delete[] vz; // velocities in z-dir
	delete[] vy; // velocities in y-dir
	delete[] vx; // velocities in x-dir
	delete[] T11;
	delete[] T22;
	delete[] T33;
	delete[] T12;
	delete[] T23;
	delete[] T13;

	delete[] d;  // density
	delete[] muu;
	delete[] lmd;

	delete[] B;  // Boundary array

//	delete[] trans;	// Transducer array TODO: Why does this cause code to hang up at the end?
}

#endif /* SPACE_H_ */
