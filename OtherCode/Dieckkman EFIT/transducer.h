/*
'transducer.h'
Custom transducer class for the Cartesian EFIT simulation
Cleaned up and modified version of 'transducer.h' (Bertoncini, Campbell-Leckey, Miller, etc)

Eric A. Dieckman (WM)
Last edited: 20 Sept 2011 EAD
*/

class transducer
{

private:
  double *drive;   // array that holds drive function
  int dflen;       // length of drivefunc

public:
  double posi1;            // transducer center (x-direction)
  double posi2;            // transducer center (y-direction)
  double posi3;            // transducer center (z-direction)
  double radius;           // transducer radius - meters
  bool driven;             // driven: true=active (pitch or pitch/catch), false=passive (catch)
  int transID;
  int numnodes;            // number of nodes in simulation space
  double *record;  // array that holds recorded value

  transducer() {driven = false;} // blank constructer

  transducer(double x1, double x2, double x3, double rad, int tID, int maxt)
  {
    posi1  = x1;
    posi2  = x2;
    posi3  = x3;
    radius = rad;
    transID = tID;
    driven = false;
    dflen=0;
    record = new double[maxt];
    for (int i = 0; i< maxt; i++) record[i] = 0;
  }

  ~transducer() {} // blank deconstructor

// ===================================================================================
// Initialize (define array and dimensions - call before using!)
// ===================================================================================
  void setDriveFunction(int len, double df[])
  {
    //drive = new double[len];
    drive = df;
    dflen = len;
    driven = true;
    return;
  }
  
  double drivef(int t) 
  {
    if (t<dflen)
      return drive[t];
    else
      return 0;
  }
};
