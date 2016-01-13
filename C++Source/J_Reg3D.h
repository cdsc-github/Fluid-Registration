//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//
#ifndef __J_Reg3D__
#define __J_Reg3D__
	
#include "DoubleArray3D.h"
#include "Grid3D.h"


class J_Reg3D
{

public:

	char reg_type;

	long fidelity;
	double triangle_u;
	double lambda;
	long smoothing;

	long heat_eqn_TimeSteps;

	DoubleArray3D J;	// Jacobian values |Dh|

	DoubleArray3D K;	// Gaussian kernel
	double StD;			// standard deviation
	double Parzen_StD;	// std for dGaussian/dy for Parzen windowing kernel

	double Jmin;
	double Jmax;

	double SSD;
	double MI;
	double SKL;
	double KL;

	double totalTime;

	int BCtype;

	J_Reg3D();

	void evaluate_f_L2( const DoubleArray3D& T, const DoubleArray3D& S,
						const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
						const Grid3D& grid,
						DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3 );

	void evaluate_f_MI( const DoubleArray3D& T, const DoubleArray3D& S,
						const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
						const Grid3D& grid,
						DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3 );

	void evaluate_v( DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3,
					 const Grid3D& grid );

	void update_U( DoubleArray3D& u1, DoubleArray3D& u2, DoubleArray3D& u3,
				   const DoubleArray3D& v1, const DoubleArray3D& v2, const DoubleArray3D& v3,
				   const Grid3D& grid );

	DoubleArray3D linearInterpolation( const DoubleArray3D& T, 
									   const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
									   const Grid3D& grid );

	DoubleArray3D getJacobian( const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3, 
							   const Grid3D& grid );

	void calculateEnergy_L2( const Grid3D& grid, double& energy_L2 );

	void calculateEnergy_MI( const Grid3D& grid, double& energy_MI );

	void outputJinv( long kk, long factor, const Grid3D& grid );

	void outputParameters( const Grid3D& grid, const long TimeSteps, const long outputCount );
};

#endif
