//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//
#ifndef __PDE_System3D__
#define __PDE_System3D__


class PDE_System3D
{
	
public:

	long problemNumber;		// problem identifier

	virtual void operatorL( const DoubleArray3D& v1, 
							const DoubleArray3D& v2, 
							const DoubleArray3D& v3,
							DoubleArray3D& Lop1, 
							DoubleArray3D& Lop2, 
							DoubleArray3D& Lop3, 
							const double& dx ) = 0;

	virtual void solveExactly( DoubleArray3D& v1,
							   DoubleArray3D& v2,
							   DoubleArray3D& v3,
							   const DoubleArray3D& f1,
							   const DoubleArray3D& f2,
							   const DoubleArray3D& f3,
							   const double& dx )        = 0;
	
	virtual void applyRelaxationGS( DoubleArray3D& v1, 
									DoubleArray3D& v2,
									DoubleArray3D& v3,
									const DoubleArray3D& f1,
									const DoubleArray3D& f2, 
									const DoubleArray3D& f3, 
									const double& dx )   = 0;

	virtual void findResidual( const DoubleArray3D& v1, 
							   const DoubleArray3D& v2,
							   const DoubleArray3D& v3,
							   const DoubleArray3D& f1,
							   const DoubleArray3D& f2,
							   const DoubleArray3D& f3,
							   DoubleArray3D& R1,
							   DoubleArray3D& R2,
							   DoubleArray3D& R3,
							   const double& dx )        = 0;
};

#endif
