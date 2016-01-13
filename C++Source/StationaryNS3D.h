//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//
#ifndef __StationaryNS3D__
#define __StationaryNS3D__

#include "DoubleArray3D.h"
#include "PDE_System3D.h"


class StationaryNS3D : public PDE_System3D
{
	
public:

	double mu;
	double lambda;
	
	StationaryNS3D();

	void operatorL( const DoubleArray3D& v1, 
					const DoubleArray3D& v2, 
					const DoubleArray3D& v3,
					DoubleArray3D& Lop1, 
					DoubleArray3D& Lop2, 
					DoubleArray3D& Lop3, 
					const double& dx );

	void solveExactly( DoubleArray3D& v1,
					   DoubleArray3D& v2,
					   DoubleArray3D& v3,
					   const DoubleArray3D& f1,
					   const DoubleArray3D& f2,
					   const DoubleArray3D& f3,
					   const double& dx );
	
	void applyRelaxationGS( DoubleArray3D& v1, 
							DoubleArray3D& v2,
							DoubleArray3D& v3,
							const DoubleArray3D& f1,
							const DoubleArray3D& f2, 
							const DoubleArray3D& f3, 
							const double& dx );

	void findResidual( const DoubleArray3D& v1, 
					   const DoubleArray3D& v2,
					   const DoubleArray3D& v3,
					   const DoubleArray3D& f1,
					   const DoubleArray3D& f2,
					   const DoubleArray3D& f3,
					   DoubleArray3D& R1,
					   DoubleArray3D& R2,
					   DoubleArray3D& R3,
					   const double& dx );

};

#endif
