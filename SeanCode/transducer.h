/*
'transducer.h'
Custom transducer class for the Cartesian EFIT simulation
Cleaned up and modified version of 'transducer.h' (Bertoncini, Campbell-Leckey, Miller, etc)

Eric A. Dieckman (WM)
Last edited: 10 Apr 2019 SMR
*/

// TODO: Rule of 5: Add a copy constructor, copy assignment operator, move constructor, and move assignment
using namespace std;
class transducer{
//private:
public:
	double *drive{nullptr};	// array that holds drive function
	int dflen=0;			// length of drivefunc

//public:
	double posi1=0;			// transducer center (z-direction)
	double posi2=0;			// transducer center (y-direction)
	double posi3=0;			// transducer center (x-direction)
	double radius=0;		// transducer radius - meters
	bool driven=false;		// driven: true=active (pitch or pitch/catch), false=passive (catch)
	int transID=0;
	int numelems=0;			// number of elements in simulation space
	double *record{nullptr};// array that holds recorded value

	double theta=0; 		// beam angle from normal / x-direction - radians
	double phi=0; 			// beam angle from (+)z-direction - radians
	int num_y = 0;			// y elements in this worker
	int num_z = 0;			// z elements in this worker
	int *tau{nullptr}; 		// yz array of element time delays for drive function
	int zbeg = 0;			// starting location for this worker in z-direction
//	int tau_size = 0;		// number of elements in the tau array (aka how many transducer elements are on this node)
	int tau_values=0;		// number of elements which have values in the tau array (equal to number of transducer elements on this node)
	double dtods=0;

	transducer() { // blank constructor
//		driven = false;
	}

	transducer(double x1, double x2, double x3, double rad, int tID, int maxt, double angle1, double angle2, int dim_y, int dim_z, int z_start, double timebyspace){
		posi1  = x1;
		posi2  = x2;
		posi3  = x3;
		radius = rad;
		transID = tID;
		driven = false;
//		dflen=0;
		record = new double[maxt];
		for (int i = 0; i< maxt; i++) record[i] = 0;

		theta = angle1;
		phi = angle2;
		num_y = dim_y;
		num_z = dim_z;
		zbeg = z_start;
		dtods=timebyspace;
//		for(int k=0; k<num_z; k++){
//			for(int j=0; j<num_y; j++){
//				if(((k+zbeg-posi1)*(k+zbeg-posi1)+(j-posi2)*(j-posi2)) <= (radius*radius)){
//					tau_size++; //executes if (z-zo)^2+(y-yo)^2<=rad^2, ie, if within a circle of transducer radius
//				}
//			}
//		}
//		tau = new int[tau_size];

		tau = new int[num_z*num_y];

//		numnodes = 0;
	}

	~transducer() { // blank deconstructor
//		delete[] drive;
//		delete[] record;
//		delete[] tau;
	}

// ===================================================================================
// Initialize (define array and dimensions - call before using!)
// ===================================================================================
	void setDriveFunction(int len, double df[], int min_tau_absolute){
		// debug line
//		int tempz, tempy;
//		for(int tau_count=0; tau_count<tau_size; tau_count++){
		for(int tau_count=0; tau_count<num_z*num_y; tau_count++){
			tau[tau_count]=tau[tau_count]-min_tau_absolute;
			// debug lines
//			tempz=tau_count%num_z;
//			tempy=(tau_count/num_z)%num_y;
//			cout<<"@(tau_count)=("<<tau_count<<"), (y,z)=("<<tempy<<","<<tempz<<"), tau="<<tau[tau_count]<<endl;
		}
		drive = new double[len];
		drive = df;
		dflen = len;
		driven = true;
		return;
	}

	double drivef(int t){
		if (t<dflen){
			return drive[t];
		}
		else{
			return 0;
		}
	}

	double setDelays(int speed){
		int cL = speed;		// default speed of sound, calculated in slave()
		int min_tau_local=0;	// minimum value of tau in this node
		int temp_tau;
//		int tau_count=0; 	// tracks location in tau array
		//debug line
//		cout<<"num_z="<<num_z<<", num_y="<<num_y<<", zbeg="<<zbeg<<", posi1="<<posi1<<", posi2="<<posi2<<", radius="<<radius<<", cL="<<cL<<endl;

		for(int k=0; k<num_z; k++){
			for(int j=0; j<num_y; j++){
				if(((k+zbeg-posi1)*(k+zbeg-posi1)+(j-posi2)*(j-posi2)) <= (radius*radius)){
					temp_tau=static_cast<int>( floor(((k+zbeg-posi1)*cos(phi) + (j-posi2)*sin(phi)) * sin(theta)/cL/dtods) );
//					tau[tau_count]=temp_tau;
					tau[k*num_y+j]=temp_tau;
					if(temp_tau<min_tau_local){
						min_tau_local=temp_tau;
						// debug line
//						cout<<"@(y,z)=("<<j<<","<<k<<"), tau="<<temp_tau<<endl;

					}
//					tau_count++;
					tau_values++;
				}
				else{
					tau[k*num_y+j]=0; // if not within the circle of the transducer
				}
			}
		}
		return min_tau_local;
	}
};
