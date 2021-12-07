/* efit.cpp
* EFIT simulation - development copy
*
* Sean M. Raley (UNH)
* based on code by Eric A. Dieckman (WM)
* 27 March 2019
* Last edited: 27 March 2019 SMR
*
* log
*/

//#include "pch.h" // for Microsoft VS only

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include "space.h"

using namespace std;

void master();
void slave();
int* DistributeSimulationParameters();
void DistributeTransducers(int *xposs);
void dump3Dbin(int t);

int mpirank, numworkers;
int maxt, outputevery, maxz, numtransducers;
int* EvenVolDims = new int[3];

/* TODO: delete?
//int whohasaline = 0;
int recordalineat = 50; // z-point to save pressure data (x and y points set to middle of sim space) - should be less than div
*/
int main(int argc, char *argv[]) // initialize MPI
{
	MPI_Init(NULL, NULL); // TODO: change? to EFIT_cart syntax: MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size(MPI_COMM_WORLD, &numworkers);  /* get number of nodes */
	numworkers--;

	if (mpirank == 0)
		master();
	else
		slave();

	MPI_Finalize();
	return 0;
}

// ==========================================================================
// Master node! -- Distributes simulation space and receives data for output
// ==========================================================================
void master() {
	time_t start, end;
	time(&start);

	int *zstartpos = new int[numworkers];

//	MPI_Status  status; // TODO: delete?
	cout << "Master node is online! \n";

	zstartpos = DistributeSimulationParameters(); // Initialize each node
	DistributeTransducers(zstartpos);

	// writes text file for binary->bov bash command
	string fname = "bin2bovCmd.ascii";
	ofstream outFile(fname.c_str(), ios::out);
	outFile<<"bash makebovs ";
	outFile<<maxt<<" ";
	outFile<<outputevery<<" ";
	for(int e=2;e>=0;e--)
		outFile<<EvenVolDims[e]<<" ";
	outFile.close();

	/* TODO: delete?
	ofstream pfile("pressuredata.ascii", ios::out);
	double al = 0;
	*/

	// TODO: probably replace with lines 71-79 from EFIT_cart
	for (int t = 0; t < maxt; t++) {
		//MPI_Recv(&al, 1, MPI_DOUBLE, whohasaline, 858, MPI_COMM_WORLD, &status);
		//outFile << al << " ";
		//cout << "Time in Master node is: " << t << endl;

		/* TODO: delete?
		if (outputevery == 1) { // Save pressure data at specific time point
			MPI_Recv(&al, 1, MPI_DOUBLE, 1, 858, MPI_COMM_WORLD, &status);
			// pulls from 1st node - recordalineat must be < div 
			pfile << al << " ";
			cout << "Saved single pressure measurement at time: " << t << "\n";
		}
		*/

		if (t%outputevery == 0 && outputevery != 1) {
			//dump3Dascii(t);
			//cout << "Saved pressure data as ASCII at time: " << t << "\n";
			dump3Dbin(t);
			cout << "Saved pressure data as binary at time: " << t << "\n";
			//dump3Dvtk(t);
			//cout << "Saved pressure data as vtk at time: " << t << "\n";
		}
	}

//	pfile.close(); TODO: delete?

	time(&end);
	printf("Total Run Time: %.2lf seconds\n", difftime(end, start));
	return;
}

// ========================================================================
// Slave node! -- Does the grunt work
// ========================================================================
void slave() {
	// --- Receive sim parameters from master and initialize ---
	MPI_Status  status;
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

//	int tosend0 = 0; // to use when sending a 0 in MPI is desired

	double simparams[11];
	MPI_Recv(&simparams, 11, MPI_DOUBLE, 0, 201, MPI_COMM_WORLD, &status);

	space simspace(simparams);
	maxt    	  = static_cast<int>(simparams[9]); // total number of time steps
	outputevery   = static_cast<int>(simparams[10]); // output every this many time steps
	double cL = sqrt((simparams[6]+2*simparams[7])/simparams[5]); // default speed of sound, calculated

	//	int simspace.Lyx = simspace.num_y*simspace.num_x; // max size of y*x dim
	//
	if (mpirank == 1){ // node is on left
		simspace.type = 1;
		cout<<"slave is type "<<simspace.type<<endl;
	}
	else if (mpirank == numworkers){ // node is on right
		simspace.type  = 3;
		cout<<"slave is type "<<simspace.type<<endl;
	}
	else{
		simspace.type  = 2; // node is in middle
		cout<<"slave is type "<<simspace.type<<endl;
	}

/* TODO: delete? should be taken care of by transducer.h, assuming input strategy isn't changed
		// --- Receive drive function ---
	if (mpirank == 1) {
		double *drive = new double[maxt];
		MPI_Recv(&drive[0], maxt, MPI_DOUBLE, 0, 202, MPI_COMM_WORLD, &status);
		simspace.df = drive;
	}
*/
	// -- Receive reflector parameters ---
	int nr;
	double *rpars = new double[9];
	MPI_Recv(&nr, 1, MPI_INT, 0, 203, MPI_COMM_WORLD, &status);

	for (int i = 0; i < nr; i++) {
		MPI_Recv(&rpars[0], 9, MPI_DOUBLE, 0, 204, MPI_COMM_WORLD, &status);
		simspace.addReflector(rpars[0], rpars[1], rpars[2], static_cast<int>(rpars[3]), static_cast<int>(rpars[4]), rpars[5], rpars[6], rpars[7],rpars[8]);
	}
	delete[] rpars;
	cout << "num_z in slave(mpirank="<<mpirank<<") is: " << simspace.num_z << endl;

	// --- Receive transducer parameters ---
	double tparams[8];
	bool done = false;

	while (done == false){
		MPI_Recv(&tparams, 8, MPI_DOUBLE, 0, 211, MPI_COMM_WORLD, &status);
		int drvlen=static_cast<int>(tparams[4]);
		int min_tau_local, min_tau_global;

		if (tparams[0] == -1) done = true;
		else{
			transducer t(tparams[0],tparams[1],tparams[2],tparams[3],static_cast<int>(tparams[5]),maxt,tparams[6],tparams[7],simspace.num_y,simspace.num_z,simspace.zbeg,simspace.dtods);

			if (drvlen > 0){
				double *drive = new double[drvlen];
				MPI_Recv(&drive[0], drvlen, MPI_DOUBLE, 0, 212, MPI_COMM_WORLD, &status);
				min_tau_local = t.setDelays(cL);
				MPI_Send(&min_tau_local, 1, MPI_INT, 0, 213, MPI_COMM_WORLD);
				// debug line
//				cout << "Slave "<<mpirank<<" Debug Code 3.6: min_tau_local="<<min_tau_local << endl;
				MPI_Recv(&min_tau_global, 1, MPI_INT, 0, 214, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// debug line
//				cout << "Slave "<<mpirank<<" Debug Code 3.7: min_tau_global="<<min_tau_global << endl;

				t.setDriveFunction(drvlen,drive,min_tau_global);

	//			delete[] drive; //TODO: deleting this array is causing the drive function to be 0!!!!!
			}
			simspace.addTransducer(t);
		}
	}

	// --- Run simulation ---
	for (int t = 0; t < maxt; t++) {
		if (mpirank == 1 && t%10==0)
			cout << " timestep: " << t << "     " << simspace.num_z << ", " << simspace.num_y << ", " << simspace.num_x << endl;
		simspace.drivetime = t;
		//debug line
//		if(mpirank==2) cout<<"drivetime = "<<simspace.drivetime<<"; drivef = "<<simspace.trans[0].drivef(simspace.drivetime)<<endl;
	    simspace.UpdateTransducers(t);

/* TODO: delete?
		if (outputevery == 1) {
			double al = simspace.val(simspace.pp, recordalineat, simspace.num_y / 2, simspace.num_x / 2); // AT WHICH POINT TO SAVE PRESSURE DATA (Z,Y,X)?
			MPI_Send(&al, 1, MPI_DOUBLE, 0, 858, MPI_COMM_WORLD);
		}
*/
		if (t%outputevery == 0) { // sends output to master node
			int len = simspace.GetEvenVolLen(); // get even num_z length for each node
//			double* x = simspace.GetEvenVol(simspace.vx, len);
			double* x = simspace.GetEvenVol(len);

			MPI_Send(&len, 1, MPI_INT, 0, 1101, MPI_COMM_WORLD);
			MPI_Send(&x[0], len, MPI_DOUBLE, 0, 1102, MPI_COMM_WORLD);

			delete[] x;
		}

	    // --- Update V's ---
	    simspace.UpdateVs(1,1);					// Update left boundary
	    simspace.UpdateVs(simspace.num_z-2,simspace.num_z-2);    	// Update right boundary

	    if (mpirank>1){							// send left
	    	MPI_Isend(&simspace.vy[simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 301, MPI_COMM_WORLD, &request1);
	    	MPI_Isend(&simspace.vx[simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 302, MPI_COMM_WORLD, &request2);
	    }

	    if (mpirank<numworkers){				// send  right
	    	MPI_Isend(&simspace.vz[(simspace.num_z-2)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 303, MPI_COMM_WORLD, &request3);
	    }

	    simspace.UpdateVs(2,simspace.num_z-3);	// update inner nodes

	    if (mpirank<numworkers){				// receive  from right
			MPI_Recv(&simspace.vy[(simspace.num_z-1)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 301, MPI_COMM_WORLD, &status1);
			MPI_Recv(&simspace.vx[(simspace.num_z-1)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 302, MPI_COMM_WORLD, &status2);
	    }

	    if (mpirank>1){							// receive from left
	    	MPI_Recv(&simspace.vz[0], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 303, MPI_COMM_WORLD, &status3);
	    }


	    if (mpirank>1){ 							// Wait for sends to complete
	    	MPI_Wait(&request1,&status1);
	    	MPI_Wait(&request2,&status2);
	    }

	    if (mpirank<numworkers){
	      MPI_Wait(&request3, &status3);
	    }

	    // --- Update T's ---
	    simspace.UpdateTs(1,1);                        // Update left boundary
	    simspace.UpdateTs(simspace.num_z-2,simspace.num_z-2);    // Update right boundary

	    if (mpirank>1){                                // send Trz, Tzp left
	    	MPI_Isend(&simspace.T11[simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 311, MPI_COMM_WORLD, &request4);
	    }

	    if (mpirank<numworkers){                       // send Tzz, Trp right
	    	MPI_Isend(&simspace.T12[(simspace.num_z-2)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 313, MPI_COMM_WORLD, &request5);
	    	MPI_Isend(&simspace.T13[(simspace.num_z-2)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 314, MPI_COMM_WORLD, &request6);
	    }

	    simspace.UpdateTs(2,simspace.num_z-3);              // update inner nodes

	    if (mpirank<numworkers){                       // reveive Trz, Tzp from right
	    	MPI_Recv(&simspace.T11[(simspace.num_z-1)*simspace.Lyx], simspace.Lyx, MPI_DOUBLE, (mpirank+1), 311, MPI_COMM_WORLD, &status4);
	    }

	    if (mpirank>1){                                // receive Tzz, Trp from left
	    	MPI_Recv(&simspace.T12[0], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 313, MPI_COMM_WORLD, &status5);
	    	MPI_Recv(&simspace.T13[0], simspace.Lyx, MPI_DOUBLE, (mpirank-1), 314, MPI_COMM_WORLD, &status6);
	    }


	    if (mpirank>1){ 				// Wait for sends to complete
	    	MPI_Wait(&request4, &status4);
	    }

	    if (mpirank<numworkers){
	    	MPI_Wait(&request5, &status5);
	    	MPI_Wait(&request6, &status6);
	    }
	}
	return;
}

// ========================================================================
// Reads in parameter file (in.file), distributes to all workers, 
// and divides up the simulation space  
// ========================================================================
int* DistributeSimulationParameters() {
	char inputFilename[] = "in.file";
	ifstream inFile;
	inFile.open("in.file", ios::in);

	if (!inFile) {
		cerr << "Can't open input file " << inputFilename << endl;
		exit(1);
	}

	double *simparams = new double[11];
	inFile >> simparams[0];     //maxz
	inFile >> simparams[1];     //maxy
	inFile >> simparams[2];     //maxx
	inFile >> simparams[3];     //ds
	inFile >> simparams[4];     //dt
	inFile >> simparams[5];     //default den
	inFile >> simparams[6];     //lm - default Lame constant - lambda
	inFile >> simparams[7];     //mu - default Lame constant - mu

	inFile >> simparams[9];     //maxt
	inFile >> simparams[10];    //outevery

	maxt = static_cast<int>(simparams[9]);
	outputevery = static_cast<int>(simparams[10]);
	maxz = static_cast<int>(simparams[0]);

	// Send initial data to each node
	int div, divaccum = 0;
	int* xpos = new int[numworkers];
	EvenVolDims[0]=0;
	EvenVolDims[1]=static_cast<int>(simparams[1]/2);
	EvenVolDims[2]=static_cast<int>(simparams[2]/2);

	for (int n = 1; n <= numworkers; n++) {
		div = (maxz/numworkers);
		if ((n-1) <= (maxz%numworkers))
			div++;   // divide space along z direction
		simparams[0] = static_cast<double>(div);
		cout << "Divided space (div) is " << div << "\n";
		simparams[8] = static_cast<double>(divaccum); // tells the worker where its starting z location is
		cout << "Worker's starting locations (divaccum) is " << divaccum << "\n";
		MPI_Send(&simparams[0], 11, MPI_DOUBLE, n, 201, MPI_COMM_WORLD);
	    xpos[n-1] = static_cast<int>(simparams[8]);
		divaccum = divaccum + div;

		if(div%2==0)
			EvenVolDims[0]+=(div/2); // len if num_z is even
		else
			EvenVolDims[0]+=((div-1)/2); // len if num_z is even

/* TODO: delete?
		//if ((whohasaline==0)&&(divaccum>=recordalineat))
		//	whohasaline=n;
*/

	}

	cout << "Total simulation timesteps (maxt) = " << maxt << endl;
	delete[] simparams;

/* TODO: delete?
	//cout << "whohasaline = " << whohasaline << "\n";

	// --- Read in drive function and send to worker number 1 ---
	double *drive = new double[maxt];
	for (int i = 0; i < maxt; i++) {
		inFile >> drive[i];
	}
	MPI_Send(&drive[0], maxt, MPI_DOUBLE, 1, 202, MPI_COMM_WORLD);
*/

	// --- Read in reflectors and distribute to all workers ---
	int numref;
	inFile >> numref;
	double *rpars = new double[9];
	cout << "  Number of reflectors:  " << numref << endl;
	for (int n = 1; n <= numworkers; n++) {
		MPI_Send(&numref, 1, MPI_INT, n, 203, MPI_COMM_WORLD);
	}

	for (int i = 0; i < numref; i++) {
		inFile >> rpars[0];  // reflector type
		inFile >> rpars[1];  // reflector position in x1
		inFile >> rpars[2];  // reflector position in x2
		inFile >> rpars[3];  // reflector position in x3 - (start for cylinder)
		inFile >> rpars[4];  // reflector position in x3 - (end for cylinder)
		inFile >> rpars[5];  // refector radius
		inFile >> rpars[6];  // refector density
		inFile >> rpars[7];  // refector mu
		inFile >> rpars[8];  // refector lambda

		for (int n = 1; n <= numworkers; n++) {
			MPI_Send(&rpars[0], 9, MPI_DOUBLE, n, 204, MPI_COMM_WORLD);
		}
	}
	delete[] rpars;

	inFile.close();
	return xpos;
}

// ===================================================================================
// Reads in transducer file (trans.file) and distributes to the correct workers
// ===================================================================================
void DistributeTransducers(int *zposs){
	double tparams[8];
	int numtrans, worker;
	int min_tau=0;
	int min_tau_recv=0;

	char inputFilename[] = "trans.file";
	ifstream inFile;
	inFile.open("trans.file", ios::in);

	if (!inFile){
		cerr << "Can't open input file " << inputFilename << endl;
		exit(1);
	}

	inFile >> numtrans;
	cout << "  number of transducers: " << numtrans << endl;
	numtransducers = numtrans;

	for (int tr = 0; tr<numtrans; tr++){
		inFile >> tparams[0];  // transducer z location
		inFile >> tparams[1];  // transducer y location
		inFile >> tparams[2];  // transducer x location
		inFile >> tparams[3];  // transducer radius
		inFile >> tparams[4];  // len of drive function
		int drvlen = static_cast<int>(tparams[4]);

		tparams[5] = static_cast<double>(tr);

		inFile >> tparams[6];  // transducer theta
		inFile >> tparams[7];  // transducer phi


		double *drive{nullptr};
		if (tparams[4]>0){
			drive = new double[drvlen];
			for (int i = 0; i<drvlen; i++){
				inFile >> drive[i];
			}
		}

		// --- Figure out which workers get the transducer ---
		worker = 0;
		for (int tosend = 1; tosend<numworkers; tosend++)
		  if (tparams[0] >= zposs[tosend-1] && tparams[0] < zposs[tosend]) worker = tosend;
		  if (tparams[0] >= zposs[numworkers-1] && tparams[0] < maxz) worker = numworkers;
		  else if (worker == 0) cout << "error: transducer postion not found: zpos - " <<  tparams[0] << ", " << zposs[numworkers-1] << ", " << maxz << endl;

		// --- Send transducer info to worker ---
		if (worker > 1){
			if ((tparams[0] - tparams[3]) <= zposs[worker-1]){
				MPI_Send(&tparams[0], 8, MPI_DOUBLE, worker-1, 211, MPI_COMM_WORLD);
				if (tparams[4]>0){
					MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-1, 212, MPI_COMM_WORLD);
					MPI_Recv(&min_tau_recv, 1, MPI_INT, worker-1, 213, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if(min_tau_recv<min_tau) min_tau = min_tau_recv;
					MPI_Send(&min_tau, 1, MPI_INT, worker-1, 214, MPI_COMM_WORLD);
				}
			}
			if ((tparams[0] - tparams[3]) <= zposs[worker-2]){
				MPI_Send(&tparams[0], 8, MPI_DOUBLE, worker-2, 211, MPI_COMM_WORLD);
				if (tparams[4]>0){
					MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker-2, 212, MPI_COMM_WORLD);
					MPI_Recv(&min_tau_recv, 1, MPI_INT, worker-2, 213, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if(min_tau_recv<min_tau) min_tau = min_tau_recv;
					MPI_Send(&min_tau, 1, MPI_INT, worker-2, 214, MPI_COMM_WORLD);
				}
			}
		}

		MPI_Send(&tparams[0], 8, MPI_DOUBLE, worker, 211, MPI_COMM_WORLD);
		if (tparams[4]>0){
			MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker, 212, MPI_COMM_WORLD);
			MPI_Recv(&min_tau_recv, 1, MPI_INT, worker, 213, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(min_tau_recv<min_tau) min_tau = min_tau_recv;
			MPI_Send(&min_tau, 1, MPI_INT, worker, 214, MPI_COMM_WORLD);
		}
		if (worker < numworkers){
			if ((tparams[0] + tparams[3]) >= zposs[worker]){ // TODO: should this be zposs[worker+1] ?
				MPI_Send(&tparams[0], 8, MPI_DOUBLE, worker+1, 211, MPI_COMM_WORLD);
				if (tparams[4]>0){
					MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+1, 212, MPI_COMM_WORLD);
					MPI_Recv(&min_tau_recv, 1, MPI_INT, worker+1, 213, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if(min_tau_recv<min_tau) min_tau = min_tau_recv;
					MPI_Send(&min_tau, 1, MPI_INT, worker+1, 214, MPI_COMM_WORLD);
				}
			}
			if ((tparams[0] + tparams[3]) >= zposs[worker+1]){ // TODO: should this be zposs[worker+2] ?
				MPI_Send(&tparams[0], 8, MPI_DOUBLE, worker+2, 211, MPI_COMM_WORLD);
				if (tparams[4]>0){
					MPI_Send(&drive[0], tparams[4], MPI_DOUBLE, worker+2, 212, MPI_COMM_WORLD);
					MPI_Recv(&min_tau_recv, 1, MPI_INT, worker+2, 213, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if(min_tau_recv<min_tau) min_tau = min_tau_recv;
					MPI_Send(&min_tau, 1, MPI_INT, worker+2, 214, MPI_COMM_WORLD);
				}
			}
		}
		delete[] drive;
	}

	// --- Let all workers know we are done distributing transducers ---
	tparams[0] = -1;tparams[1] = -1;tparams[2] = -1;tparams[3] = -1;tparams[4] = -1;tparams[5] = -1;
	for (int n = 1; n <= numworkers; n++){
		MPI_Send(&tparams[0], 8, MPI_DOUBLE, n, 211, MPI_COMM_WORLD);
	}

	inFile.close();
	return;
}

/* TODO: delete?
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
*/

// ========================================================================
// Dump data to file!
// ========================================================================
void dump3Dbin(int t) // pressure data as binary
{
	MPI_Status  status;
	double *data3d{nullptr};
	int len;

	stringstream strm;
	strm << t;
	string fname = "data3d_at_t_" + strm.str() + ".bin";
	ofstream outFile(fname.c_str(), ios::binary);

	for (int n = 1; n <= numworkers; n++) {
		MPI_Recv(&len, 1, MPI_INT, n, 1101, MPI_COMM_WORLD, &status);
		data3d = new double[len];
		MPI_Recv(&data3d[0], len, MPI_DOUBLE, n, 1102, MPI_COMM_WORLD, &status);

		// debug line
//		cout<<"For worker "<<n<<", binary data string length is "<<len<<endl; // for debugging noncubic output matrix
		for (int i = 0; i < len; i++) {
			outFile.write((char *)(&data3d[i]), sizeof(data3d[i]));
		}

		delete[] data3d;
	}


	outFile.close();
	return;
}
