#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include "spipe.h"

using namespace std;
const int syncevery = 100;

void master();
void slave();
int* DistributeSimulationParameters();    // sends out simulation params to workers
void DistributeTransducers(int *zposs);   // distributes transducers to the appropriate workers
void dumpTopPlate(int t);
void dumpTopPlateThreeD(int t);
void SyncNodes();

int rank, numworkers;

int maxt, max1, m2m3;  // max number of time steps
int outputevery;       // output every
int numtransducers;    //   
int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numworkers);  /* get number of nodes */
  numworkers--;

  if (rank == 0)
	  master();
  else
	  slave();

  MPI_Finalize();
  return 0;
}

// ===================================================================================
// master node! -- distribures simulation space and receives data for output
// ===================================================================================
void master()
{
   time_t start,end;
   time (&start);

   int rank, div, n, i, max1;
   int *zstartpos = new int[numworkers];  

   MPI_Status status;
   cout << "master node is online! \n";

   //------------------------------------------------------------------ 
   /* Send out initialization messages to each node */

   zstartpos = DistributeSimulationParameters(); 
   DistributeTransducers(zstartpos);
   
   for (int t=0; t<maxt; t++)
   {
	   if (t%outputevery == 0 ) 
	   {
         dumpTopPlateThreeD(t);
        cout << "Collecting Slices at time: " << t << "\n";
	   }
	   
	   
	  
   }
   
//    for (int t=0; t<maxt; t++)
//   {
//   
//    if (t%outputevery == 0 ) 
//	   {
//         dumpTopPlate(t);
//         cout << "Collecting 2D Slices at time: " << t << "\n";
//	   }
	   //if (t%syncevery == 0)
		//   SyncNodes();
//}
  

   time (&end);
   printf ("Total Run Time: %.2lf seconds\n", difftime (end,start) );
   return;
}

// ===================================================================================
// slave node! -- Does the grunt work
// ===================================================================================
void slave()
{
    // ---------------------------------------------------------------------	
	// receive simulation parameters from master and initialize pipe section
	MPI_Status  status;     MPI_Request request[2];
	MPI_Status  status1; 	MPI_Status status2;
    MPI_Status status3; 	MPI_Status status4; 
	MPI_Status status5; 	MPI_Status status6; 
    MPI_Status status9;     MPI_Status status10; 
	MPI_Request request1; 
    MPI_Request request2;
	MPI_Request request3;
	MPI_Request request4;
    MPI_Request request5;
    MPI_Request request6;
    MPI_Request request9;
    MPI_Request request10;
	double simparams[15];
    MPI_Recv(&simparams, 15, MPI_DOUBLE, 0, 201, MPI_COMM_WORLD, &status);
	    
	spipe pipe;	
	pipe.num2 = simparams[0];    // number of nodes in r direction 
    pipe.num1 = simparams[1]+2;  // number of nodes in z direction 
    pipe.num3 = simparams[2];    // number of nodes in p direction 
    pipe.ds   = simparams[3];    // spatial step size in r and z (meters)
    pipe.dp   = simparams[4];    // spatial step size in phi (radians)
    pipe.dt   = simparams[5];    // time step size (seconds)
    pipe.den  = simparams[6];    // density
    pipe.lm   = simparams[7];    // Lame constant - lambda
    pipe.mu   = simparams[8];    // Lame constant - mu
    pipe.rbeg = simparams[10];   // pipe inner radius (in ds units)
    pipe.zbeg = simparams[11];   // pipe starting z position for each node (set to zero)
    maxt      = simparams[12];   // number of time steps
	outputevery = simparams[13]; // output every time steps
	max1      = simparams[14];   // total number of z across entire simulation
	m2m3 = simparams[0]*simparams[2];
int max2,max3;
max2=simparams[0];
max3=simparams[2];

	if (rank == 1) pipe.pipetype = 1;
	else if (rank == numworkers) pipe.pipetype = 3;
	else pipe.pipetype = 2;



	//cout << pipe.numr << " " << pipe.numz << " " << pipe.nump << " " << pipe.curvem[1] <<" 99 \n";
	pipe.Init();
	
	pipe.setmaterialprops();
//cout<<"after set mat props \n";
	
//-- Receive Reflectors
          int nr;   double *rpars = new double[9];
	MPI_Recv(&nr, maxt, MPI_INT, 0, 203, MPI_COMM_WORLD, &status);
	for (int i = 0; i < nr; i++)
	{
		MPI_Recv(&rpars[0], 9, MPI_DOUBLE, 0, 204, MPI_COMM_WORLD, &status);
		
    	    pipe.addReflector(rpars[0],rpars[1],rpars[2],rpars[3],rpars[4],rpars[5],rpars[6],rpars[7],rpars[8]);
        
	}
	
	
	//
//	 stringstream strmst; strmst <<rank ;
//    string fname = "reflector_props" +strmst.str()+ ".ascii";
//  	ofstream outFile(fname.c_str(), ios::out);
//  	int hn1,hn2,hn3;
//  	for(hn1=0;hn1<max1; hn1++)
//  	for(hn2=0;hn2<max2; hn2++)
//  	for(hn3=0;hn3<max3; hn3++){
//  	 outFile << pipe.d.val(hn1,hn2,hn3) << " ";
//    }
  
//	outFile.close();
  	
  	

	// ---------------------------------------------------------------------
	// receive transducer parameters from master and add them to the pipe
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
		   pipe.addTransducer(t);
		}
	 }

  

    // ---------------------------------------------------------------------
	// perform simulation
	double *toplate; int len;
		double *toplate2; int len2;
	for (int t = 0; t<maxt; t++)
	{
		



        if (rank == 1 && t%10==0) cout << "timestep: " << t << " " << pipe.num1 << " " << pipe.num2 << " " << pipe.num3 <<"\n";
	    	    
		pipe.time=t;
        pipe.UpdateTransducers(t);

		// ------  Send Output to Master -------------
// BELOW is for 2D output , uncomment for that purpose
//--------------------------------
//	if (t%outputevery == 0) 
//		{
//            len2 = pipe.v3.slice_fix2_count();
//			toplate2 = new double[len2];
//			toplate2 = pipe.v3.slice_fix2(pipe.num3-5);
//			MPI_Send(&len2, 1, MPI_INT, 0, 401, MPI_COMM_WORLD);
//			MPI_Send(&toplate2[0], len2, MPI_DOUBLE, 0, 402, MPI_COMM_WORLD);
//			delete toplate2;
	//	}
		
		
		//BELOW is for 3D output
		//--------------------------
		 if (t%outputevery == 0) 
            {

            time_t start,end;

            time (&start);
            //tosend = ar.pp;
            int len = pipe.v3.GetEvenVolLen(pipe.zbeg);  //start at beginning in x for each node
            double* x=pipe.v3.GetEvenVol(pipe.zbeg);
        
            MPI_Send(&len, 1, MPI_INT, 0, 1101, MPI_COMM_WORLD);
            MPI_Send(&x[0], len, MPI_DOUBLE, 0, 1102, MPI_COMM_WORLD);
	
            time (&end);
            printf (" Total time for sending output to master node: %.2lf seconds\n",
                    difftime (end,start) );
            delete [] x;  // dynamically allocated by GetEvenVol()
            }


		// ------  Update V's -------------
		pipe.UpdateVs(1,1);                        // Update left boundary
		pipe.UpdateVs(pipe.num1-2,pipe.num1-2);    // Update right boundary 
		


		if (rank>1)  {                              // send  left
		   MPI_Isend(&pipe.v2.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 301, MPI_COMM_WORLD,
                     &request1);
		   MPI_Isend(&pipe.v3.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 302, MPI_COMM_WORLD,
                     &request2);
         }
		   
		if (rank<numworkers)                       // send  right
           {  		
           MPI_Isend(&pipe.v1.a[(pipe.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 303,
                     MPI_COMM_WORLD, &request3);
           }

		pipe.UpdateVs(2,pipe.num1-3);              // update inner nodes

		if (rank<numworkers)                       // receive  from right
           {
           MPI_Recv(&pipe.v2.a[(pipe.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 301,
                    MPI_COMM_WORLD, &status1);
           MPI_Recv(&pipe.v3.a[(pipe.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 302,
                    MPI_COMM_WORLD, &status2);	
           }

        if (rank>1)                                // receive from left
           {	
           MPI_Recv(&pipe.v1.a[0], m2m3, MPI_DOUBLE, (rank-1), 303, MPI_COMM_WORLD, &status3);
           }
	    
	    // Wait for sends to complete
        if (rank>1) {
           MPI_Wait(&request1,&status1);
           MPI_Wait(&request2,&status2);
           }
        if (rank<numworkers) {
           MPI_Wait(&request3, &status3);
           }


        //pipe.doABCs(max1);
       
        // ------  Update T's -------------
		pipe.UpdateTs(1,1);                        // Update left boundary
		pipe.UpdateTs(pipe.num1-2,pipe.num1-2);    // Update right boundary 


		if (rank>1)                                // send Trz, Tzp left
           {
           MPI_Isend(&pipe.T11.a[m2m3], m2m3, MPI_DOUBLE, (rank-1), 311, MPI_COMM_WORLD,
                     &request4);			
           }
		if (rank<numworkers)                       // send Tzz, Trp right
           {
           MPI_Isend(&pipe.T12.a[(pipe.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 313,
                     MPI_COMM_WORLD, &request5);
           MPI_Isend(&pipe.T13.a[(pipe.num1-2)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 314
                     , MPI_COMM_WORLD, &request6);
           }
		
		pipe.UpdateTs(2,pipe.num1-3);              // update inner nodes    
		
		if (rank<numworkers)                       // reveive Trz, Tzp from right
           {
           MPI_Recv(&pipe.T11.a[(pipe.num1-1)*m2m3], m2m3, MPI_DOUBLE, (rank+1), 311,
                    MPI_COMM_WORLD, &status4);
           }
		if (rank>1)                                // receive Tzz, Trp from left
           {
           MPI_Recv(&pipe.T12.a[0], m2m3, MPI_DOUBLE, (rank-1), 313, MPI_COMM_WORLD, &status5);
           MPI_Recv(&pipe.T13.a[0], m2m3, MPI_DOUBLE, (rank-1), 314, MPI_COMM_WORLD, &status6);
           }
		
		
        if (rank>1) {
           MPI_Wait(&request4,&status4);
           }
        if (rank<numworkers) {
           MPI_Wait(&request5,&status5);
           MPI_Wait(&request6, &status6);
           }

	}


	return;
}


// ===================================================================================
// Reads in parameter file and distributes parameters to all workers.  This is also
// where the simulation space is divided up.  
// ===================================================================================
int* DistributeSimulationParameters()
{
	char inputFilename[] = "in.file";
	ifstream inFile;
	inFile.open("in.file", ios::in);
        
	if (!inFile) {
       cerr << "Can't open input file " << inputFilename << endl;
       exit(1);
	}
    
	double *simparams = new double[15];

	inFile >> simparams[0];   //pipe.num2;   // number of nodes in r direction 
    inFile >> simparams[1];   //pipe.num1;   // number of nodes in z direction 
    inFile >> simparams[2];   //pipe.num3;   // number of nodes in p direction 
    inFile >> simparams[3];   //pipe.ds;     // spatial step size in r and z (meters)
    inFile >> simparams[4];   //pipe.dp;     // spatial step size in phi (radians)
    inFile >> simparams[5];   //pipe.dt;     // time step size (seconds)
    inFile >> simparams[6];   //pipe.den;    // density
    inFile >> simparams[7];   //pipe.lm;     // Lame constant - lambda
    inFile >> simparams[8];   //pipe.mu;     // Lame constant - mu
    inFile >> simparams[10];  //pipe.rbeg;   // pipe inner radius (in ds units)
    inFile >> maxt;           // number of time steps
    inFile >> outputevery;    // number of nodes in x3 direction 
    simparams[12] = maxt; simparams[13] = outputevery;
	simparams[14] = simparams[1];
	m2m3 = simparams[0]*simparams[2];

    max1 = simparams[1];


    // send initial data to each node
	int div, divaccum = 0; 
	int* zpos = new int[numworkers];
    for (int n = 1; n <= numworkers; n++)
    {
	  div  = (max1/(numworkers)); if ((n-1)< (max1%(numworkers))) div++;   /* divide space along x1 direstion */
	  simparams[1]  = div;
	  cout<<"div is "<<div<<"\n";
      simparams[11] = divaccum;    // tells the worker where its starting z location is
      cout<<"divaccum is "<<divaccum<<"\n";
	  MPI_Send(&simparams[0], 15, MPI_DOUBLE, n, 201, MPI_COMM_WORLD); 

	 
	  zpos[n-1] = simparams[11]; divaccum = divaccum+div;
    }


// read in reflectors and distribute to all workers
	// ------------------------------------------------------------------
	int numref; inFile >> numref;
	double *rpars = new double[9];
	cout <<"  Number of reflectors:  " << numref << endl;
	
	for (int n = 1; n <= numworkers; n++)
       MPI_Send(&numref, 1, MPI_INT, n, 203, MPI_COMM_WORLD); 

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
   		inFile >> rpars[8];    // reflector lambda
	    for (int n = 1; n <= numworkers; n++)
		    MPI_Send(&rpars[0], 9, MPI_DOUBLE, n, 204, MPI_COMM_WORLD); 
      
	}

	inFile.close();

	return zpos;
}
// ===================================================================================
// Reads in transducer file and distributes transducers to the correct workers. 
// ===================================================================================
void DistributeTransducers(int *zposs)
{
	double *drive; 
	double tparams[6];
	int drivelen, numtrans, worker;
	
	char inputFilename[] = "trans.file";
	ifstream inFile;
	inFile.open("trans.file", ios::in);
        
	if (!inFile) {
       cerr << "Can't open input file " << inputFilename << endl;
       exit(1);
	}

    inFile >> numtrans;
	cout << "  number of transducers: " << numtrans << endl;
	numtransducers = numtrans;

	for (int tr = 0; tr<numtrans; tr++)
	{
      inFile >> tparams[0];  // tposz;    // transducer z location
	  inFile >> tparams[1];  // tposr;    // transducer r location
      inFile >> tparams[2];  // tposp;    // transducer p location 
      inFile >> tparams[3];  // trad;     // transducer radius
	  inFile >> tparams[4];  // drivelen; // len of drive function
	  tparams[5] = tr;
	  
      int drvlen;
      drvlen=tparams[4];
      
      if (tparams[4]>0)
	  {
	    drive = new double[drvlen];
        for (int i = 0; i<drvlen; i++)
		  inFile >> drive[i];
	  }

	  // find which worker gets the transducer
      worker = 0;
      for (int tosend = 1; tosend<numworkers; tosend++)
	    if (tparams[0] >= zposs[tosend-1] && tparams[0] < zposs[tosend]) worker = tosend;
	  if (tparams[0] >= zposs[numworkers-1] && tparams[0] < max1) worker = numworkers;
      else if (worker == 0) cout << "error: transducer postion not found: zpos - " <<  tparams[0] << ", " << zposs[numworkers-1] << ", " << max1 << endl;
      
	  // send the transducer info to worker
	  if (worker > 1){
          if ((tparams[0] - tparams[3]) <= zposs[worker])
	  {
	  		  MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker-1, 211, MPI_COMM_WORLD);
	          if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-1, 212, MPI_COMM_WORLD);
		  }
		  
		        if ((tparams[0] - tparams[3]) <= zposs[worker-1])
	  {
	  		  MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker-2, 211, MPI_COMM_WORLD);
	          if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-2, 212, MPI_COMM_WORLD);
		  }
       }
		  
		  

	  MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker, 211, MPI_COMM_WORLD);
	  if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker, 212, MPI_COMM_WORLD);



 	  if (worker < numworkers){
          if ((tparams[0] + tparams[3]) >= zposs[worker+1])
		  {
			  MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker+1, 211, MPI_COMM_WORLD);
              if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+1, 212, MPI_COMM_WORLD);
		  }
		  
		  if ((tparams[0] + tparams[3]) >= zposs[worker+2])
		  {
			  MPI_Send(&tparams[0], 6, MPI_DOUBLE, worker+2, 211, MPI_COMM_WORLD);
              if (tparams[4]>0) MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+2, 212, MPI_COMM_WORLD);
		  }
		  
       }

	}
    
	// send all workers a message letting them know we are done distributing transducers
	tparams[0] = -1;tparams[1] = -1;tparams[2] = -1;tparams[3] = -1;tparams[4] = -1;tparams[5] = -1;
	for (int n = 1; n <= numworkers; n++)
      MPI_Send(&tparams[0], 5, MPI_DOUBLE, n, 211, MPI_COMM_WORLD);

	inFile.close();
    delete drive;

 	return;
}

void SyncNodes()
{
   int s; MPI_Status  status;
   for (int n = 1; n <= numworkers; n++)
	 MPI_Recv(&s, 1, MPI_INT, n, 721, MPI_COMM_WORLD, &status);
   for (int n = 1; n <= numworkers; n++)
	 MPI_Send(&n, 1, MPI_INT, n, 722, MPI_COMM_WORLD);
   cout << "  nodes synced \n"; 
}



//  dump topplate
void dumpTopPlate(int t)
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
	       outFile << topplate2[i] << " ";
	delete [] topplate2;

	}

	outFile.close();

	return;    
}

void dumpTopPlateThreeD(int t)
{
	MPI_Status  status;
    double *topplate;
    int len;

    stringstream strm; strm << t;
    string fname = "Threetoplate_at_t" +strm.str()+ ".ascii";
  	ofstream outFile(fname.c_str(), ios::out);

    for (int n = 1; n <= numworkers; n++)
       { 
		MPI_Recv(&len, 1, MPI_INT, n, 1101, MPI_COMM_WORLD, &status);

        //	if (n==1) topplate = new double[len];
        topplate = new double[len];

		MPI_Recv(&topplate[0], len, MPI_DOUBLE, n, 1102, MPI_COMM_WORLD, &status);
		
        for (int i = 0; i < len; i++)
		  outFile << topplate[i] << " ";
       delete [] topplate;		   
       }

 
	outFile.close();
//cout<<"after outfile close before 2nd delete\n";
	
//cout<<"immediately before end of dumpTopPlate\n";
	return;    

}
