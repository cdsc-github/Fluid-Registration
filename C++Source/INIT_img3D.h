#ifndef __INIT_img3D__
#define __INIT_img3D__

#include "DoubleArray3D.h"
#include "Grid3D.h"


void outputU( const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3, 
			  const long& kk, 
			  const long bx, const long by, const long bz );

void initialize( DoubleArray3D& u, const Grid3D& grid, const long IMGnumber, const long volN );

#endif
