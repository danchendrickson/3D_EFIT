class transducer
{
  private:
	double *drive;   // array that holds drive function

    int dflen;       // length of drivefunc

  public:

	double posi1;            // transducer center (r-direction) - meters
	double posi2;            // transducer center (z-direction) - meters
	double posi3;            // transducer center (p-direction) - angle

    double radius;           // transducer radius - meters
	bool driven;             // driven = true  - active (pitch or pitch/catch)
	                         //        = false - passive (catch)
	int transID;
    int numnodes;            // number of nodes in simulation space

	double *record;  // array that holds recorded value

	// Blank Constructer
    transducer() {driven = false;}

    //
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

    // Blank Deconstructor
    ~transducer() {}

    // Init - defines the array and its dimensions - MUST BE CALLED BEFORE USING
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
