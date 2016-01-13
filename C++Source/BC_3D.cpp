#include <iostream>
using namespace std;

#include "BC_3D.h"

//
//#####################################################################
//						BC_3D.cpp
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: Dec. 15, 2006
//
//#####################################################################
//

void applyBC( DoubleArray3D& A, const long width, const long type )
{
	switch( type )
	{
	case 0:		// None
		break;

	case 1:		// Dirichlet
		DirichletBC( A,width );
		break;

	case 5:		// Constant extrapolation
		NeumannBC( A, width );
		break;

	case 6:		// Linear extrapolation
		NeumannBC_linear( A, width );
		break;

	case 10:	// Periodic
		PeriodicBC( A, width );
		break;

	case 11:	// Periodic in X and Z, linear extrapolation in Y
		Periodic_X_Z_Extrapolation_Y( A, width );
		break;
	
	default: cout << "Incorrect Specification of Boundary Conditions " << endl;
		exit(1);
		break;
	}
}

//#####################################################################

/*
void InvCurvatureFlow3D::NeumannBC( DoubleArray3D& A )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i;  long j;  long k;

	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			A( 0, j,k) = A( 1, j,k);
			A(m-1,j,k) = A(m-2,j,k);
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			A(i, 0, k) = A(i, 1, k);
			A(i,n-1,k) = A(i,n-2,k);
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			A(i,j, 0 ) = A(i,j, 1 );
			A(i,j,p-1) = A(i,j,p-2);
		}
	}

	A( 0 , 0 , 0 ) = A( 1 , 1 , 1 );
	A(m-1, 0 , 0 ) = A(m-2, 1 , 1 );
	A( 0 ,n-1, 0 ) = A( 1 ,n-2, 1 );
	A( 0 , 0 ,p-1) = A( 1 , 1 ,p-2);
	A(m-1,n-1, 0 ) = A(m-2,n-2, 1 );
	A(m-1, 0 ,p-1) = A(m-2, 1 ,p-2);
	A( 0 ,n-1,p-1) = A( 1 ,n-2,p-2);
	A(m-1,n-1,p-1) = A(m-2,n-2,p-2);
}
*/

//#####################################################################

void DirichletBC( DoubleArray3D& A, const long width )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i, j, k;
	
	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			for( i = 0; i < width; i++ )
			{	A(i,j,k) = 0;	}
			for( i = m-width; i < m; i++ )
			{	A(i,j,k) = 0;	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			for( j = 0; j < width; j++ )
			{	A(i,j,k) = 0;	}

			for( j = n-width; j < n; j++ )
			{	A(i,j,k) = 0;	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for( k = 0; k < width; k++ )
			{	A(i,j,k) = 0;	}

			for( k = p-width; k < p; k++ )
			{	A(i,j,k) = 0;	}
		}
	}
}

//#####################################################################

void NeumannBC( DoubleArray3D& A, const long width )
{
	//
	// constant extension:
	//
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i;  long j;  long k;

	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			for( i = 0; i < width; i++ )
			{	A(i,j,k) = A( width,j,k);		}
			for( i = m-width; i < m; i++ )
			{	A(i,j,k) = A((m-1)-width,j,k);	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			for( j = 0; j < width; j++ )
			{	A(i, j, k) = A(i, width, k);	}

			for( j = n-width; j < n; j++ )
			{	A(i,j,k) = A(i,(n-1)-width,k);	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for( k = 0; k < width; k++ )
			{	A(i, j, k) = A(i,j, width );	}

			for( k = p-width; k < p; k++ )
			{	A(i,j,k) = A(i,j,(p-1)-width);	}
		}
	}

	for( i=0; i < width; i++ )
	{
		A(  i  ,  i  ,  i  ) = A( width     ,   width    ,   width    );
		A(m-1-i,  i  ,  i  ) = A((m-1)-width,   width    ,   width    );
		A(  i  ,n-1-i,  i  ) = A( width     , (n-1)-width,   width    );
		A(  i,    i  ,p-1-i) = A( width     ,   width    , (p-1)-width);
		A(m-1-i,n-1-i,  i  ) = A((m-1)-width, (n-1)-width,   width    );
		A(m-1-i,  i  ,p-1-i) = A((m-1)-width,   width    , (p-1)-width);
		A(  i  ,n-1-i,p-1-i) = A( width     , (n-1)-width, (p-1)-width);
		A(m-1-i,n-1-i,p-1-i) = A((m-1)-width, (n-1)-width, (p-1)-width);
	}
}

//#####################################################################

void NeumannBC_linear( DoubleArray3D& A, const long width )
{
	//
	// linear extrapolation:
	//
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i;  long j;  long k;

	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			for( i = 0; i < width; i++ )
			{	A(i,j,k) = 2*A(width,j,k) - A(2*width-i,j,k);			}

			for( i = m-width; i < m; i++ )
			{	A(i,j,k) = 2*A((m-1)-width,j,k) - A( 2*((m-1)-width)-i, j, k );	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			for( j = 0; j < width; j++ )
			{	A(i,j,k) = 2*A(i,width,k) - A(i,2*width-j,k);			}

			for( j = n-width; j < n; j++ )
			{	A(i,j,k) = 2*A(i,(n-1)-width,k) - A( i, 2*((n-1)-width)-j, k );	}
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for( k = 0; k < width; k++ )
			{	A(i,j,k) = 2*A(i,j, width ) - A(i,j,2*width-k);	}

			for( k = p-width; k < p; k++ )
			{	A(i,j,k) = 2*A(i,j,(p-1)-width) - A(i,j,2*((p-1)-width)-k );	}
		}
	}
}

//#####################################################################

void PeriodicBC( DoubleArray3D& A, const long width )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i, j, k;

	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			for( i = 0; i < width; i++ )
			{
				A(i,j,k) = A( m-1-2*width+i , j, k );
			}
			for( i = m-width-1; i < m; i++ )
			{
				A(i,j,k) = A( i+2*width-m+1 , j, k );
			}
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			for( j = 0; j < width; j++ )
			{
				A(i,j,k) = A( i , n-1-2*width+j, k );
			}
			for( j = n-width-1; j < n; j++ )
			{
				A(i,j,k) = A( i , j+2*width-n+1, k );
			}
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for( k = 0; k < width; k++ )
			{	
				A(i,j,k) = A( i, j, p-1-2*width+k );
			}

			for( k = p-width-1; k < p; k++ )
			{
				A(i,j,k) = A( i, j, k+2*width-p+1 );
			}
		}
	}
}

//#####################################################################

void Periodic_X_Z_Extrapolation_Y( DoubleArray3D& A, const long width )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i, j, k;

	// Periodic in X direction:
	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			for( i = 0; i < width; i++ )
			{
				A(i,j,k) = A( m-1-2*width+i , j, k );
			}
			for( i = m-width-1; i < m; i++ )
			{
				A(i,j,k) = A( i+2*width-m+1 , j, k );
			}
		}
	}

	// Linear extrapolation in Y direction:
	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			for( j = 0; j < width; j++ )
			{	A(i,j,k) = 2*A(i,width,k) - A(i,2*width-j,k);			}

			for( j = n-width; j < n; j++ )
			{	A(i,j,k) = 2*A(i,(n-1)-width,k) - A( i, 2*((n-1)-width)-j, k );	}
		}
	}

	// Periodic in Z direction:
	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for( k = 0; k < width; k++ )
			{	
				A(i,j,k) = A( i, j, p-1-2*width+k );
			}

			for( k = p-width-1; k < p; k++ )
			{
				A(i,j,k) = A( i, j, k+2*width-p+1 );
			}
		}
	}
}

//#####################################################################
