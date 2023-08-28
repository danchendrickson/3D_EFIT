/*
'EFIT_cart.cpp'
Run EFIT simulation on a cartesian grid
Cleaned up and modified version of 'EFIT.cpp' (Bertoncini, Campbell-Leckey, Miller, etc)
NOT COMPATIBLE WITH PREVIOUS INPUT FILES (uses x,y,z instead of y,x,z)

DEPENDS: mpi.h, setup_cart_space.h, array3d.h, array3d_int.h, transducer.h, in.file, trans.file

Eric A. Dieckman (WM)
Last edited: 21 Sept 2011 EAD
*/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include "setup_cart_space.h"

using namespace std;
//const int syncevery = 100;

void master();
void slave();
int* DistributeSimulationParameters();
void DistributeTransducers(int *xposs);
void dumpTopPlate(int t);
void dumpTopPlateThreeD(int t);
void SyncNodes();

int rank, numworkers;

int max1, maxt;
int outputevery;
int numtransducers;

int main(int argc, char *argv[]) // initialize MPI
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numworkers);  // get number of nodes
  numworkers--;

  if (rank == 0)
    master();
  else
    slave();

  MPI_Finalize();
  return 0;
}

// ===================================================================================
// Master node! -- Distributes simulation space and receives data for output
// ===================================================================================
void master()
{
  time_t start,end;
  time (&start);

  int *xstartpos = new int[numworkers];  

//  MPI_Status status;
  cout << "Master node is online! \n";

  xstartpos = DistributeSimulationParameters(); // Initialize each node
  DistributeTransducers(xstartpos);
  
  for (int t=0; t<maxt; t++)
  {
    if (t%outputevery == 0) 
    {
      dumpTopPlateThreeD(t); // for 3D simulation
      //dumpTopPlate(t); // for 2D simulation
      cout << "Collecting Slices at time: " << t << "\n";
    }
  }

  time (&end);
  printf ("Total Run Time: %.2lf seconds\n", difftime (end,start));
  return;
}


// ===================================================================================
// Slave node! -- Does the grunt work
// ===================================================================================
void slave()
{
  // --- Receive simulation parameters from master and initialize ---
  MPI_Status status;
  MPI_Status status1;
  MPI_Status status2;
  MPI_Status status3;
  MPI_Status status4; 
  MPI_Status status5;
  MPI_Status status6; 
  MPI_Request request1; 
  MPI_Request request2;
  MPI_Request request3;
  MPI_Request request4;
  MPI_Request request5;
  MPI_Request request6;

  double simparams[11];
  MPI_Recv(&simparams, 11, MPI_DOUBLE, 0, 201, MPI_COMM_WORLD, &status);
	    
  setup_cart_space simspace;	
  simspace.num1 = simparams[0]+2; // number of nodes in x direction 
  simspace.num2 = simparams[1];   // number of nodes in y direction 
  simspace.num3 = simparams[2];   // number of nodes in z direction 
  simspace.ds   = simparams[3];   // spatial step size in x and y (meters)
  
  simspace.dt   = simparams[4];   // time step size (seconds)
  simspace.den  = simparams[5];   // density
  simspace.lm   = simparams[6];	  // Lame constant - lambda
  simspace.mu   = simparams[7];   // Lame constant - mu
  
  simspace.xbeg = simparams[8];   // simspace starting x position for each node (set to zero)
  maxt      	= simparams[9];   // number of time steps
  outputevery 	= simparams[10];  // output every time steps

  int max2,max3, m2m3;
  max2 = simparams[1];
  max3 = simparams[2];
  m2m3 = max2*max3;

  if (rank == 1) simspace.pipetype = 1;
  else if (rank == numworkers) simspace.pipetype = 3;
  else simspace.pipetype = 2;

  simspace.Init();
	
  simspace.setmaterialprops();
	

  // --- Receive reflector parameters ---
  int nr;   double *rpars = new double[9];
  MPI_Recv(&nr, 1, MPI_INT, 0, 203, MPI_COMM_WORLD, &status);
  for (int i = 0; i < nr; i++)
  {
    MPI_Recv(&rpars[0], 9, MPI_DOUBLE, 0, 204, MPI_COMM_WORLD, &status);
    simspace.addReflector(rpars[0],rpars[1],rpars[2],rpars[3],rpars[4],rpars[5],rpars[6],rpars[7],rpars[8]);
  }

//	 stringstream strmst; strmst <<rank ;
//    string fname = "reflector_props" +strmst.str()+ ".ascii";
//  	ofstream outFile(fname.c_str(), ios::out);
//  	int hn1,hn2,hn3;
//  	for(hn1=0;hn1<max1; hn1++)
//  	for(hn2=0;hn2<max2; hn2++)
//  	for(hn3=0;hn3<max3; hn3++){
//  	 outFile << simspace.d.val(hn1,hn2,hn3) << " ";
//    }
  
//	outFile.close();
	
	

  // --- Receive transducer parameters ---
  double tparams[6];
  bool done = false;

  while (done == false)
  {
    MPI_Recv(&tparams, 6, MPI_DOUBLE, 0, 211, MPI_COMM_WORLD, &status);
    //cout << "  " << tparams[5] << "  " << rank << "  \n";
    int drvlen;
    drvlen=tparams[4];
    if (tparams[0] == -1) done = true;
    else
    {
      transducer t(tparams[0],tparams[1],tparams[2],tparams[3],tparams[5],maxt);
      //transducer t(tparams[0],5,tparams[2],tparams[3],(int)tparams[5], maxt);
      if (drvlen > 0)
      {
	double *drive = new double[drvlen];
	MPI_Recv(&drive[0], tparams[4], MPI_DOUBLE, 0, 212, MPI_COMM_WORLD, &status);
	t.setDriveFunction(tparams[4],drive);
      }
	simspace.addTransducer(t);
    }
  }


  // --- Perform simulation ---
  // double *toplate2; int len2; // Needed for 2D case
  for (int t = 0; t<maxt; t++)
  {
    if (rank == 1 && t%10==0) cout << "timestep: " << t << " " << simspace.num1 << " " << simspace.num2 << " " << simspace.num3 << "\n";
    simspace.time=t;
    simspace.UpdateTransducers(t);

    //Send output to master node (2D)
    /*
    if (t%outputevery == 0) 
    {
      len2 = simspace.v3.slice_fix2_count();
      toplate2 = new double[len2];
      toplate2 = simspace.v3.slice_fix2(simspace.num3-5);
      
      MPI_Send(&len2, 1, MPI_INT, 0, 401, MPI_COMM_WORLD);
      MPI_Send(&toplate2[0], len2, MPI_DOUBLE, 0, 402, MPI_COMM_WORLD);
      delete toplate2;
    }
    */
		
    // Send output to master node (3D)
    if (t%outputevery == 0) 
    {
      time_t start,end;
      time (&start);
      //tosend = ar.pp;
      int len = simspace.v3.GetEvenVolLen(simspace.xbeg);  //start at beginning in x for each node
      double* x=simspace.v3.GetEvenVol(simspace.xbeg);
	  
      MPI_Send(&len, 1, MPI_INT, 0, 1101, MPI_COMM_WORLD);
      MPI_Send(&x[0], len, MPI_DOUBLE, 0, 1102, MPI_COMM_WORLD);
	
      time (&end);
      printf (" Total time for sending output to master node: %.2lf seconds\n",
      difftime (end,start) );
      delete [] x;  // dynamically allocated by GetEvenVol()
    }


    // --- Update V's ---
    simspace.UpdateVs(1,1);                        	// Update left boundary
    simspace.UpdateVs(simspace.num1-2,simspace.num1-2);    	// Update right boundary 
		
    if (rank>1)					// send left
    {
      MPI_Isend(&simspace.v2.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 301, MPI_COMM_WORLD, &request1);
      MPI_Isend(&simspace.v3.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 302, MPI_COMM_WORLD, &request2);
    }
		  
    if (rank<numworkers)			// send  right
    {
      MPI_Isend(&simspace.v1.a[(simspace.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 303, MPI_COMM_WORLD, &request3);
    }

    simspace.UpdateVs(2,simspace.num1-3);              // update inner nodes

    if (rank<numworkers)                       // receive  from right
    {
      MPI_Recv(&simspace.v2.a[(simspace.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 301, MPI_COMM_WORLD, &status1);
      MPI_Recv(&simspace.v3.a[(simspace.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 302, MPI_COMM_WORLD, &status2);	
    }

    if (rank>1)                                // receive from left
    {	
      MPI_Recv(&simspace.v1.a[0], m2m3, MPI_DOUBLE, (rank-1), 303, MPI_COMM_WORLD, &status3);
    }
	
    
    if (rank>1) 				// Wait for sends to complete
    {
      MPI_Wait(&request1,&status1);
      MPI_Wait(&request2,&status2);
    }
	
    if (rank<numworkers) 
    {
      MPI_Wait(&request3, &status3);
    }

    //simspace.doABCs(max1);


    // --- Update T's ---
    simspace.UpdateTs(1,1);                        // Update left boundary
    simspace.UpdateTs(simspace.num1-2,simspace.num1-2);    // Update right boundary 

    if (rank>1)                                // send Trz, Tzp left
    {
      MPI_Isend(&simspace.T11.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 311, MPI_COMM_WORLD, &request4);			
    }

    if (rank<numworkers)                       // send Tzz, Trp right
    {
      MPI_Isend(&simspace.T12.a[(simspace.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 313, MPI_COMM_WORLD, &request5);
      MPI_Isend(&simspace.T13.a[(simspace.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 314, MPI_COMM_WORLD, &request6);
    }
    
    simspace.UpdateTs(2,simspace.num1-3);              // update inner nodes    
		
    if (rank<numworkers)                       // reveive Trz, Tzp from right
    {
      MPI_Recv(&simspace.T11.a[(simspace.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 311, MPI_COMM_WORLD, &status4);
    }

    if (rank>1)                                // receive Tzz, Trp from left
    {
      MPI_Recv(&simspace.T12.a[0], m2m3, MPI_DOUBLE, (rank-1), 313, MPI_COMM_WORLD, &status5);
      MPI_Recv(&simspace.T13.a[0], m2m3, MPI_DOUBLE, (rank-1), 314, MPI_COMM_WORLD, &status6);
    }
		
		
    if (rank>1) 				// Wait for sends to complete
    {
      MPI_Wait(&request4, &status4);
    }

    if (rank<numworkers) 
    {
      MPI_Wait(&request5, &status5);
      MPI_Wait(&request6, &status6);
    }

  }
  return;
}


// ===================================================================================
// Reads in parameter file (in.file), distributes parameters to all workers, 
// and divides up the simulation space  
// ===================================================================================
int* DistributeSimulationParameters()
{
  char inputFilename[] = "in.file";
  ifstream inFile;
  inFile.open("in.file", ios::in);

  if (!inFile) 
  {
    cerr << "Can't open input file " << inputFilename << endl;
    exit(1);
  }
  
  double *simparams = new double[11];

  inFile >> simparams[0];   //simspace.num1;   // number of nodes in x direction 
  inFile >> simparams[1];   //simspace.num2;   // number of nodes in y direction 
  inFile >> simparams[2];   //simspace.num3;   // number of nodes in z direction 
  inFile >> simparams[3];   //simspace.ds;     // spatial step size (meters)
 
  inFile >> simparams[4];   //simspace.dt;     // time step size (seconds)
  inFile >> simparams[5];   //simspace.den;    // density
  inFile >> simparams[6];   //simspace.lm;     // Lame constant - lambda
  inFile >> simparams[7];   //simspace.mu;     // Lame constant - mu
  inFile >> maxt;   			      // Number of time steps
  inFile >> outputevery;  // outputevery     // Output every so many time steps
  
  simparams[9] = maxt;
  simparams[10] = outputevery;  
  max1 = simparams[0];
  
  // --- Send initial data to each node ---
  int div, divaccum = 0; 
  int* xpos = new int[numworkers];
  for (int n = 1; n <= numworkers; n++)
  {
    div = (max1/(numworkers));
    if ((n-1) < (max1%(numworkers))) div++; // Divide space along x (num1) direction
    simparams[0] = div;
    cout<<"Divided space (div) is "<<div<<"\n";
    simparams[8] = divaccum; // tells the worker where its starting x location is
    cout<<"Worker's starting locations (divaccum) is "<<divaccum<<"\n";
    MPI_Send(&simparams[0], 11, MPI_DOUBLE, n, 201, MPI_COMM_WORLD); 
    xpos[n-1] = simparams[8];
    divaccum = divaccum+div;
  }


  // --- Reads and distributes reflectors to all workers ---
  int numref; inFile >> numref;
  double *rpars = new double[9];
  cout <<"  Number of reflectors:  " << numref << endl;
  for (int n = 1; n <= numworkers; n++)	
  {
    MPI_Send(&numref, 1, MPI_INT, n, 203, MPI_COMM_WORLD); 
  }
  
  for (int i = 0; i < numref; i++)
  {
    inFile >> rpars[0];  // reflector type
    inFile >> rpars[1];  // reflector position in x1
    inFile >> rpars[2];  // reflector position in x2
    inFile >> rpars[3];  // reflector position in x3 - (start for cylinder)
    inFile >> rpars[4];  // reflector position in x3 - (end for cylinder)
    inFile >> rpars[5];  // refector radius
    inFile >> rpars[6];  // refector density
    inFile >> rpars[7];  // refector mu
    inFile >> rpars[8];  // reflector lambda
    for (int n = 1; n <= numworkers; n++)
    {  
      MPI_Send(&rpars[0], 9, MPI_DOUBLE, n, 204, MPI_COMM_WORLD); 
    }
  }

  inFile.close();
  return xpos;
}

// ===================================================================================
// Reads in transducer file (trans.file) and distributes to the correct workers
// ===================================================================================
void DistributeTransducers(int *xposs)
{
  double *drive; 
  double tparams[6];
  int numtrans, worker;
	
  char inputFilename[] = "trans.file";
  ifstream inFile;
  inFile.open("trans.file", ios::in);
	
  if (!inFile) 
  {
    cerr << "Can't open input file " << inputFilename << endl;
    exit(1);
  }

  inFile >> numtrans;
  cout << "  number of transducers: " << numtrans << endl;
  numtransducers = numtrans;

  for (int tr = 0; tr<numtrans; tr++)
  {
    inFile >> tparams[0];  // transducer x location
    inFile >> tparams[1];  // transducer y location
    inFile >> tparams[2];  // transducer z location 
    inFile >> tparams[3];  // transducer radius
    inFile >> tparams[4];  // len of drive function
    tparams[5] = tr;
    int drvlen;
    drvlen=tparams[4];
    
    if (tparams[4]>0)
    {
      drive = new double[drvlen];
      for (int i = 0; i<drvlen; i++)
      {
	inFile >> drive[i];
      }
    }

    // --- Figure out which workers get the transducer ---
    worker = 0;
    for (int tosend = 1; tosend<numworkers; tosend++)
      if (tparams[0] >= xposs[tosend-1] && tparams[0] < xposs[tosend]) worker = tosend;
      if (tparams[0] >= xposs[numworkers-1] && tparams[0] < max1) worker = numworkers;
      else if (worker == 0) cout << "error: transducer postion not found: xpos - " <<  tparams[0] << ", " << xposs[numworkers-1] << ", " << max1 << endl;

    // --- Send transducer info to worker ---
    if (worker > 1) 
    {  
      if ((tparams[0] - tparams[3]) <= xposs[worker-1])
      {
	MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker-1, 211, MPI_COMM_WORLD);
	if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-1, 212, MPI_COMM_WORLD);
      }
	
      if ((tparams[0] - tparams[3]) <= xposs[worker-2])
      {
	MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker-2, 211, MPI_COMM_WORLD);
	if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-2, 212, MPI_COMM_WORLD);
      }
    }

    MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker, 211, MPI_COMM_WORLD);
    if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker, 212, MPI_COMM_WORLD);

    if (worker < numworkers)
    {
      if ((tparams[0] + tparams[3]) >= xposs[worker])
      {	
	MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker+1, 211, MPI_COMM_WORLD);
	if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+1, 212, MPI_COMM_WORLD);
      }

      if ((tparams[0] + tparams[3]) >= xposs[worker+1])
      {	
	MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker+2, 211, MPI_COMM_WORLD);
	if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+2, 212, MPI_COMM_WORLD);
      }
    }
  }
    
  // --- Let all workers know we are done distributing transducers ---
  tparams[0] = -1;tparams[1] = -1;tparams[2] = -1;tparams[3] = -1;tparams[4] = -1;tparams[5] = -1;
  for (int n = 1; n <= numworkers; n++)
  {
    MPI_Send(&tparams[0], 5, MPI_DOUBLE, n, 211, MPI_COMM_WORLD);
  }

  inFile.close();
  delete drive;
  return;
}


// ===================================================================================
// Sync nodes OUTDATED????
// ===================================================================================
void SyncNodes()
{
  int s; MPI_Status  status;
  for (int n = 1; n <= numworkers; n++)
	MPI_Recv(&s, 1, MPI_INT, n, 721, MPI_COMM_WORLD, &status);
  for (int n = 1; n <= numworkers; n++)
	MPI_Send(&n, 1, MPI_INT, n, 722, MPI_COMM_WORLD);
  cout << "  nodes synced \n"; 
}


// ===================================================================================
// Dump data to file!
// ===================================================================================
void dumpTopPlate(int t) // For 2D data; illustrates outputting ASCII files
{
  MPI_Status  status; 
  double *topplate2;
  int len2;

  stringstream strm; strm << t;
  string fname = "toplate_at_t" +strm.str()+ ".ascii";
  ofstream outFile(fname.c_str(), ios::out);

  for (int n = 1; n <= numworkers; n++)
  {
    MPI_Recv(&len2, 1, MPI_INT, n, 401, MPI_COMM_WORLD, &status);
    topplate2 = new double[len2];
    MPI_Recv(&topplate2[0], len2, MPI_DOUBLE, n, 402, MPI_COMM_WORLD, &status);
	
    //cout << " << " << h[0] << "\n";
    for (int i = 1; i <= len2; i++)
    {
      outFile << topplate2[i] << " ";
    }
    delete [] topplate2;
  }
  
  outFile.close();
  return;    
}


void dumpTopPlateThreeD(int t) // For 3D data; illustrates outputting binary data
{
  MPI_Status  status;
  double *topplate;
  int len;

  stringstream strm; strm << t;
  string fname = "Threetoplate_at_t" +strm.str()+ ".bin";
  ofstream fout(fname.c_str(), ios::binary);

  for (int n = 1; n <= numworkers; n++)
  { 
    MPI_Recv(&len, 1, MPI_INT, n, 1101, MPI_COMM_WORLD, &status);
    //	if (n==1) topplate = new double[len];
    topplate = new double[len];
    MPI_Recv(&topplate[0], len, MPI_DOUBLE, n, 1102, MPI_COMM_WORLD, &status);
    
    for (int i = 0; i < len; i++)
    {
      fout.write((char *)(&topplate[i]), sizeof(topplate[i]));
    }
    
    delete [] topplate;		   
  }
  
  fout.close();
  return;    
}
