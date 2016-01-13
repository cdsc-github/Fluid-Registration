#include <iostream>
using namespace std;

#include <math.h>
#include "DoubleArray2D.h"
#include "DoubleArray3D.h"

//#####################################################################

double findL2( const DoubleArray2D& A, const long& width )
{
	double L2norm = 0;

	long m = A.getIndex1Size();
	long n = A.getIndex2Size();

	long i; long j;

	for( j = width; j < n-width; j++ )
	{	for( i = width; i < m-width; i++ )
		{	L2norm = L2norm + A(i,j)*A(i,j);	}	}

	L2norm = sqrt( L2norm /(n-2.0*width)/(m-2.0*width) );
	
	return L2norm;
}

//#####################################################################

double findL2( const DoubleArray3D& A, const long& width )
{
	double L2norm = 0;

	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i,j,k;

	for( k = width; k < p-width; k++ )
	{	for( j = width; j < n-width; j++ )
		{	for( i = width; i < m-width; i++ )
			{	L2norm = L2norm + A(i,j,k)*A(i,j,k);	}	}	}

	L2norm = sqrt( L2norm /(n-2.0*width)/(m-2.0*width)/(p-2.0*width) );
	
	return L2norm;
}
