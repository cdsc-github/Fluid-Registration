//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//


#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <fstream>

#include <math.h>

#include "BC_3D.h"
#include "FileIO_3D.h"
#include "INIT_img3D.h"
#include "J_Reg3D.h"
#include "Convolution3D.h"
#include "common_routines.h"

#define PI 4.0*atan2(1.0,1.0)

double realtime(void);

//#####################################################################

int main(int argc, char* argv[])
{
	long factor = 256;
	bool regrid = 0;

	FileIO_3D IO;

	Grid3D grid;

	long dataset;
	cout << "Dataset:" << endl;
	cout << "2 - ADNI BASELINE  " << endl;
	cout << "3 - ADNI FOLLOWUP  " << endl;
	
	cin >> dataset;

	long volN = 0;

	long m, n, p;
	long bxi = grid.bxi, byi = grid.byi, bzi = grid.bzi;
	long bxo = grid.bxo, byo = grid.byo, bzo = grid.bzo;
	if( dataset == 2 )				// ADNI BASELINE
	{	// m = 222; n = 222; p = 222;
		m = 128; n = 160; p = 128;	}
	else if( dataset == 3 )				// ADNI FOLLOWUP
	{	// m = 220; n = 220; p = 220;	
		m = 256; n = 256; p = 256;
	}
	else
	{	cout << "Unexpected dataset!!!" << endl;	}
	
	long w = 1;

	grid.m = m;   grid.n = n;   grid.p = p;  grid.w = w;
	grid.bxi = bxi; grid.byi = byi; grid.bzi = bzi;
	grid.bxo = bxo; grid.byo = byo; grid.bzo = bzo;

	grid.dx = 1.0;  grid.dy = 1.0;  grid.dz = 1.0;

	grid.xMin = 0.0;  grid.xMax = (m-1)*grid.dx;
	grid.yMin = 0.0;  grid.yMax = (n-1)*grid.dy;
	grid.zMin = 0.0;  grid.zMax = (p-1)*grid.dz;

	DoubleArray3D  S(m,n,p);
	DoubleArray3D  T(m,n,p);


	// S = Study = Reference = Target
	// T = Template = Source

	if( dataset == 2 )
	{
		cout << "ADNI BASELINE" << endl;
		cout << "Volume # (1-10) = ";	cin >> volN;
		
		initialize( S, grid, 62, volN );	cout << "T1 to T2" << endl;
		initialize( T, grid, 61, volN );	cout << "T1 to T2" << endl;
	}
	else if( dataset == 3 )
	{
		cout << "ADNI FOLLOWUP" << endl;
		cout << "Volume # (1-10) = ";	cin >> volN;
		initialize( S, grid, 171, volN );
		initialize( T, grid, 172, volN );
	}

	findMaxMin( S );
	findMaxMin( T );
	IO.write_bin_usi( S, "S",  factor, bxo, byo, bzo );
	IO.write_bin_usi( T, "T0", factor, bxo, byo, bzo );

	DoubleArray3D u1(m,n,p);
	DoubleArray3D u2(m,n,p);
	DoubleArray3D u3(m,n,p);

	J_Reg3D REG;		

	cout << "Fidelity: " << endl;
	cout << "1 = L2" << endl;
	cout << "2 = MI" << endl;
	cin >> REG.fidelity;

	cout << "Smoothing:             " << endl;
	cout << "1 = Gaussian			" << endl;
	cout << "2 = Heat equation      " << endl;
	cout << "3 = Stationary NS      " << endl;
	cin >> REG.smoothing;

	if( REG.smoothing == 1 )
	{
		cout << "Standard deviation std (e.g. 4.5): " << endl;
		cin >> REG.StD;
		REG.K.initialize(m,n,p);
		REG.K = Gaussian3D( grid, REG.StD );
		REG.K = fftshift3D(REG.K);	
	}


	//
	//*******************************************************
	//  Specify run and output parameters
	//*******************************************************
	//

	long  TimeSteps;
	long  outputCount;

	cout << "Number of TimeSteps:     " << endl;
	cin >> TimeSteps;

	cout << "Output Every nth Step:   " << endl;
	cin >> outputCount;

	DoubleArray3D interpT = T;
	interpT = REG.linearInterpolation( T, u1, u2, u3, grid );

	REG.outputParameters( grid, TimeSteps, outputCount );

	DoubleArray3D v1(m,n,p);
	DoubleArray3D v2(m,n,p);
	DoubleArray3D v3(m,n,p);
	
	REG.J.initialize(m,n,p);	// Jacobian of the transformation

	DoubleArray1D energyL2(TimeSteps+1);
	DoubleArray1D energyMI(TimeSteps+1);

	IO.write_bin_usi( REG.J, "J0", factor, bxo, byo, bzo );

	int NumberRegrids = 0;  if(NumberRegrids != 0 ){ "No regridding implemented"; exit(1); }

	double timeStart, timeTaken;
	timeStart = realtime();

	long kk = 0;
	for( kk = 1; kk <= TimeSteps; kk++ )
	{
		if(      REG.fidelity == 1 )
		{	REG.evaluate_f_L2( interpT, S, u1, u2, u3, grid, v1, v2, v3 );	}
		else if( REG.fidelity == 2 )
		{	REG.evaluate_f_MI( interpT, S, u1, u2, u3, grid, v1, v2, v3 );	}

		REG.evaluate_v( v1, v2, v3, grid );
		REG.update_U( u1, u2, u3, v1, v2, v3, grid );

		interpT = REG.linearInterpolation( T, u1, u2, u3, grid );

		cout << "   " << REG.Jmin << " " << REG.Jmax << endl;

		REG.calculateEnergy_L2( grid, energyL2(kk-1) );
		REG.calculateEnergy_MI( grid, energyMI(kk-1) );

		if( (kk%50) == 0 || kk == 10 )
		{	
			cout << "Step " << kk << endl;

			IO.write_ascii( energyL2, "energyL2", kk );
			IO.write_ascii( energyMI, "energyMI", kk );

			//
			// BINARY OUTPUT:
			//
			outputU( u1, u2, u3, kk, bxo, byo, bzo );

			IO.write_bin_usi( interpT, "T", kk, factor, bxo, byo, bzo  );
			IO.write_bin_usi( REG.J  , "J", kk, factor, bxo, byo, bzo  );

			REG.outputJinv( kk, factor, grid );
		}

	}

	cout << endl << NumberRegrids << " regrids." << endl << endl;

	timeTaken = realtime() - timeStart;
	cout << "Filtering time  : " << timeTaken << " milliseconds.  (timeb.h)" << endl;
	cout << "Filtering time  : " << timeTaken/1000 << " seconds.      (timeb.h)" << endl;

	return 0;
}
