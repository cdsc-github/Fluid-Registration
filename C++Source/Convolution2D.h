//
//#####################################################################
//					  Convolution2D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
//
//#####################################################################
//

#ifndef __convolution2D__
#define __convolution2D__

#include "DoubleArray2D.h"
#include "Grid2D.h"


DoubleArray2D fftshift2D( const DoubleArray2D& A );

DoubleArray2D OutOfFocus2D( const Grid2D& grid, const double& radius);
DoubleArray2D OutOfFocus2D_directional( const Grid2D& grid, const double& rx, const double& ry );
DoubleArray2D RectBlur( const Grid2D& grid, const double& width, const double& height );
DoubleArray2D ComplexShapeBlur( const Grid2D& grid, const double& w1, const double& h1, const double& w2, const double& h2 );
DoubleArray2D Gaussian2D( const Grid2D& grid, const double& std);
DoubleArray2D Gaussian2D( const long m, const double std);
DoubleArray2D GeoSTAR( const Grid2D& grid );
DoubleArray2D MER_MI_kernel( const Grid2D& grid, const char kernel_choice );

DoubleArray2D dGdy( const long I, const double std);

void convolution2D( DoubleArray2D& F, const DoubleArray2D& K );

void deconvFidelity2D( DoubleArray2D& result, const DoubleArray2D& F, const DoubleArray2D& K, const DoubleArray2D& Krot, const DoubleArray2D& U );

void convGaussian2D( DoubleArray2D& F, const Grid2D& grid, const double& std );

#endif
