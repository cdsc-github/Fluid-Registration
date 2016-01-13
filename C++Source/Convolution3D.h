//
//#####################################################################
//					  Convolution3D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
//
//#####################################################################
//
#ifndef __convolution3D__
#define __convolution3D__

#include "DoubleArray3D.h"
#include "Grid3D.h"


DoubleArray3D fftshift3D( const DoubleArray3D& A );

DoubleArray3D Gaussian3D( const Grid3D& grid, const double std);

void convolution3D( DoubleArray3D& F, const DoubleArray3D& K );

void convGaussian3D( DoubleArray3D& F, const Grid3D& grid, const double& std );

#endif
