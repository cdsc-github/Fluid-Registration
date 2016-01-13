#ifndef __BC_3D__
#define __BC_3D__

#include "DoubleArray3D.h"

//
//#####################################################################
//						BC_3D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: Dec. 15, 2006
//
//#####################################################################
//

void applyBC( DoubleArray3D& A, const long width, const long type );

void DirichletBC( DoubleArray3D& A, const long width );

void NeumannBC( DoubleArray3D& A, const long width );

void NeumannBC_linear( DoubleArray3D& A, const long width );

void PeriodicBC( DoubleArray3D& A, const long width );

void Periodic_X_Z_Extrapolation_Y( DoubleArray3D& A, const long width );


#endif
