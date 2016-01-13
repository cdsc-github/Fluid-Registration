#include <algorithm>
#include <functional> 
#include <iostream>
using namespace std;

#include <math.h>

#include "common_routines.h"
//#include "BC_2D.h"

//
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: October 6, 2008
//
//#####################################################################
//

long round_to_int( double a )
{
	int int_a = (int)floor(a + 0.5);
	return int_a;
}

//#####################################################################

void display( const DoubleArray2D& A, const long& i, const long& j )
{
	cout << "A( " << i << "," << j << ") =" << A(i,j) << endl;
}

//#####################################################################

void display( const DoubleArray3D& A, const long& i, const long& j, const long& k )
{
	cout << "A( " << i << "," << j << "," << k << ") =" << A(i,j,k) << endl;
}

//#####################################################################

void findMaxMin( const DoubleArray2D& A )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();

	double AMAX = fabs(A(1,1));
	double AMIN = fabs(A(1,1));
	double Amax = A(1,1);
	double Amin = A(1,1);

	long i, j;
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			AMAX = (AMAX > fabs(A(i,j))) ? AMAX : fabs(A(i,j));
			AMIN = (AMIN < fabs(A(i,j))) ? AMIN : fabs(A(i,j));

			Amax = (Amax > A(i,j)) ? Amax : A(i,j);
			Amin = (Amin < A(i,j)) ? Amin : A(i,j);
		}
	}
	cout << "MAX = " << AMAX << ", MIN = " << AMIN << ", " << "max = " << Amax << ", min = " << Amin << endl;
}

//#####################################################################

void findMaxMin( const DoubleArray3D& A )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	double AMAX = fabs(A(1,1,1));
	double AMIN = fabs(A(1,1,1));
	double Amax = A(1,1,1);
	double Amin = A(1,1,1);

	long i, j, k;
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			for ( k = 0; k < p; k++ )
			{
				AMAX = (AMAX > fabs(A(i,j,k))) ? AMAX : fabs(A(i,j,k));
				AMIN = (AMIN < fabs(A(i,j,k))) ? AMIN : fabs(A(i,j,k));

				Amax = (Amax > A(i,j,k)) ? Amax : A(i,j,k);
				Amin = (Amin < A(i,j,k)) ? Amin : A(i,j,k);
			}
		}
	}
	cout << "MAX = " << AMAX << ", MIN = " << AMIN << ", " << "max = " << Amax << ", min = " << Amin << endl;
}

//#####################################################################

DoubleArray2D rotate180( const DoubleArray2D& A )		// rotates structure A by 180 degrees around its center
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();

	DoubleArray2D A_rotated(m,n);

	long i, j;
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			A_rotated((m-1)-i,(n-1)-j) = A(i,j);
		}
	}
	return A_rotated;
}


//#####################################################################

DoubleArray2D BilinInterp( const DoubleArray2D& T, 
						   const DoubleArray2D& X, const DoubleArray2D& Y,
						   const Grid2D& grid )
{
	long   m    = grid.m;
	long   n    = grid.n;
	long   w    = grid.w;

	double dx   = grid.dx;
	double dy   = grid.dy;

	double xMin = grid.xMin;
	double xMax = grid.xMax;
	double yMin = grid.yMin;
	double yMax = grid.yMax;

	DoubleArray2D Tx      = T;		// interpolation in x-direction
	DoubleArray2D interpT = T;

	long i, j;

    long newX, newY;
    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newX = long(floor(X(i,j)));
            if( newX < 0 )
            {	Tx(i,j) = T(0,j);	}
            else if( newX+1 >= m )
            {	Tx(i,j) = T(m-1,j);	}
            else
            {	Tx(i,j)      = T(newX,j) + ((T(newX+1,j)-T(newX,j))/dx) * (X(i,j)-newX);	}	// T(x-u1)
        }
    }

    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newY = long(floor(Y(i,j)));
            if( newY < 0 )
            {	interpT(i,j) = T(i,0);	}
            else if( newY+1 >= n )
            {	interpT(i,j) = T(i,n-1); }
            else
            {	interpT(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);	}		// T(y-u2)
        }
    }

    applyBC( interpT, w, 5 );

    return interpT;
}

//#####################################################################

DoubleArray2D BilinInterp( const DoubleArray2D& T, 
        const DoubleArray2D& X1, const DoubleArray2D& Y1,
        const double& dx )
{
    long  m = T.getIndex1Size();
    long  n = T.getIndex2Size();

    double dy = dx;

    DoubleArray2D Tx      = T;		// interpolation in x-direction
    DoubleArray2D interpT = T;


    DoubleArray2D X = X1;	X = X/dx;
    DoubleArray2D Y = Y1;	Y = Y/dx;

    long i, j;

    long newX, newY;
    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newX = long(floor(X(i,j)));
            if( newX < 0 )
            {	Tx(i,j) = T(0,j);	}
            else if( newX+1 >= m )
            {	Tx(i,j) = T(m-1,j);	}
            else
            {	
                //Tx(i,j)      = T(newX,j) + ((T(newX+1,j)-T(newX,j))/dx) * (X(i,j)-newX);	
                Tx(i,j)      = T(newX,j) + (T(newX+1,j)-T(newX,j)) * (X(i,j)-newX);
            }	// T(x-u1)
        }
    }

    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newY = long(floor(Y(i,j)));
            if( newY < 0 )
            {	interpT(i,j) = T(i,0);	}
            else if( newY+1 >= n )
            {	interpT(i,j) = T(i,n-1); }
            else
            {	
                //interpT(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);	
                interpT(i,j) = Tx(i,newY) + (Tx(i,newY+1)-Tx(i,newY)) * (Y(i,j)-newY);
            }		// T(y-u2)
        }
    }

    applyBC( interpT, 1, 5 );

    return interpT;
}

//#####################################################################

DoubleArray2D BilinInterp2( const DoubleArray2D& T, 
        const DoubleArray2D& u, const DoubleArray2D& v,
        const Grid2D& grid )
{
    long   m    = grid.m;
    long   n    = grid.n;
    long   w    = grid.w;

    double dx   = grid.dx;
    double dy   = grid.dy;

    double xMin = grid.xMin;
    double xMax = grid.xMax;
    double yMin = grid.yMin;
    double yMax = grid.yMax;

    DoubleArray2D Tx      = T;		// interpolation in x-direction
    DoubleArray2D interpT = T;

    DoubleArray2D X(m,n), Y(m,n);
    long i, j;
    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            X(i,j) = i + u(i,j);		// X = x-u1
            Y(i,j) = j + v(i,j);		// Y = y-u2
        }
    }


    long newX, newY;
    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newX = long(floor(X(i,j)));
            if( newX < 0 )
            {	Tx(i,j) = T(0,j);	}
            else if( newX+1 >= m )
            {	Tx(i,j) = T(m-1,j);	}
            else
            {	Tx(i,j)      = T(newX,j) + ((T(newX+1,j)-T(newX,j))/dx) * (X(i,j)-newX);	}	// T(x-u1)
        }
    }

    for( j = 1; j < n-1; j++ )
    {
        for( i = 1; i < m-1; i++ )
        {
            newY = long(floor(Y(i,j)));
            if( newY < 0 )
            {	interpT(i,j) = T(i,0);	}
            else if( newY+1 >= n )
            {	interpT(i,j) = T(i,n-1); }
            else
            {	interpT(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);	}		// T(y-u2)
        }
    }

    applyBC( interpT, w, 5 );

    return interpT;
}

//#####################################################################

double calculateRMSE( const DoubleArray2D& A, const DoubleArray2D& B )
{
    long m = A.getIndex1Size();
    long n = A.getIndex2Size();

    double RMSE = 0.0;

    long i, j;
    for(j = 0; j < n; j++)
    {
        for(i = 0; i < m; i++)
        {
            RMSE = RMSE + pow( A(i,j) - B(i,j), 2.0 );
        }
    }
    RMSE = sqrt(RMSE) /m/n;

    cout << "rmse = " << RMSE << endl;

    return RMSE;
}

//#####################################################################

