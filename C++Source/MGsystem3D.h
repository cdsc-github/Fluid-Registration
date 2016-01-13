//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//
#ifndef __MGsystem3D__
#define __MGsystem3D__

#include "DoubleArray3D.h"
#include "Grid3D.h"
#include "PDE_System3D.h"


class MGsystem3D
{

private:

	long nr;	// Number of  pre-smoothing relaxations,
				// i.e. while restricting the residual to coarse grid.

	long ns;	// Number of smoothing relaxations.  
				// Appropriate for two-grid, when numerous relaxations done on a coarse grid.

	long np;	// Number of post-smoothing relaxations,
				// i.e. while prolongating the error to fine grid.

	long numVcyclesFMG;	// Number of V-cycles per level for FMG.

	PDE_System3D* pPDE;	// if use virtual, need to use pointer

	void restriction( const DoubleArray3D& R, DoubleArray3D& R2, const long& order );

	void prolongation( const DoubleArray3D& V2, DoubleArray3D& V, const long& order );

	void correct( const DoubleArray3D& V, DoubleArray3D& u );

	
public:

	void setPDE( PDE_System3D* pde );

	long U_RestrictionOrder;
	long R_RestrictionOrder;

	long solver;	// identifier for GS, Two-Grid, or Multigrid.

	MGsystem3D();

	void setNumberGSRelaxations( void );
	void setNumberGSRelaxations( long preSmooth, long Smooth, long postSmooth );

	void initialize( DoubleArray3D& A, const Grid3D& grid, const long& choice );

	void MultiGridCycle( DoubleArray3D& u1, 
						 DoubleArray3D& u2,
						 DoubleArray3D& u3,
						 const DoubleArray3D& f1,
						 const DoubleArray3D& f2,
						 const DoubleArray3D& f3,
						 const double& dx );

	void MultiGridCycle_FAS( DoubleArray3D& u1, 
						 DoubleArray3D& u2,
						 DoubleArray3D& u3,
						 const DoubleArray3D& f1,
						 const DoubleArray3D& f2,
						 const DoubleArray3D& f3,
						 const double& dx );

	void outputParameters( const double& accuracy, const Grid3D& grid, 
						   const long& TimeSteps, const long& outputCount, const double& timeTaken );

};

#endif
