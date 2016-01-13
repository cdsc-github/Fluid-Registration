//
//#####################################################################
// Igor Yanovsky (C) UCLA, JPL                                
// version:  10/09/2008
//#####################################################################
//

#ifndef __common_routines__
#define __common_routines__

#include <vector>
#include "DoubleArray2D.h"
#include "DoubleArray3D.h"

#include "Grid2D.h"
#include "BC_2D.h"

long round_to_int( double a );

void display( const DoubleArray2D& A, const long& i, const long& j );
void display( const DoubleArray3D& A, const long& i, const long& j, const long& k );

void findMaxMin( const DoubleArray2D& A );
void findMaxMin( const DoubleArray3D& A );

DoubleArray2D rotate180( const DoubleArray2D& A );	// rotates structure A by 180 degrees around its center

DoubleArray2D BilinInterp( const DoubleArray2D& T, 
						   const DoubleArray2D& X, const DoubleArray2D& Y,
						   const Grid2D& grid );

DoubleArray2D BilinInterp( const DoubleArray2D& T, 
						   const DoubleArray2D& X, const DoubleArray2D& Y,
						   const double& dx );

DoubleArray2D BilinInterp2( const DoubleArray2D& T, 
						    const DoubleArray2D& u, const DoubleArray2D& v,
						    const Grid2D& grid );

double calculateRMSE( const DoubleArray2D& A, const DoubleArray2D& B );

//#####################################################################

#endif
