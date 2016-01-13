//
//#####################################################################
// Igor Yanovsky, Luminita Vese (C) UCLA, JPL
//#####################################################################
//
#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <math.h>

#include "StationaryNS3D.h"

#define PI 4.0*atan2(1.0,1.0)

//#####################################################################

StationaryNS3D::StationaryNS3D()
{
	mu = 1.0;
	lambda = 0.0;

	problemNumber = 4;

	cout << endl << "\t \t ***** 3D Stationary Navier-Stokes Equation ***** " << endl << endl;
	cout << "\t \t \t mu = " << mu << ",  lambda = " << lambda << endl << endl;
}

//#####################################################################

void StationaryNS3D::operatorL( const DoubleArray3D& v1, 
								const DoubleArray3D& v2, 
								const DoubleArray3D& v3,
								DoubleArray3D& Lop1, 
								DoubleArray3D& Lop2, 
								DoubleArray3D& Lop3, 
								const double& dx )
{
	long m = v1.getIndex1Size();
	long n = v1.getIndex2Size();
	long p = v1.getIndex3Size();

	double dx_sqrd = dx*dx;

	long i,j,k;
	
	for( k=1; k < p-1; k++ )
	{
		for( j=1; j < n-1; j++ )	// Do not update the boundary:
		{	
			for( i=1; i < m-1; i++ )
			{
				Lop1(i,j,k) =          mu *( v1(i+1,j,k) + v1(i-1,j,k) + v1(i,j+1,k) + v1(i,j-1,k) + v1(i,j,k+1) + v1(i,j,k-1) - 6.0*v1(i,j,k) ) / dx_sqrd
							+ (mu+lambda) *( v1(i+1,j,k) - 2.0*v1(i,j,k) + v1(i-1,j,k) ) / dx_sqrd
							+ (mu+lambda) *( v2(i+1,j+1,k) - v2(i-1,j+1,k) - v2(i+1,j-1,k) + v2(i-1,j-1,k) ) / (4.0*dx_sqrd)
							+ (mu+lambda) *( v3(i+1,j,k+1) - v3(i-1,j,k+1) - v3(i+1,j,k-1) + v3(i-1,j,k-1) ) / (4.0*dx_sqrd);

				Lop2(i,j,k) =          mu *( v2(i+1,j,k) + v2(i-1,j,k) + v2(i,j+1,k) + v2(i,j-1,k) + v2(i,j,k+1) + v2(i,j,k-1) - 6.0*v2(i,j,k) ) / dx_sqrd
							+ (mu+lambda) *( v2(i,j+1,k) - 2.0*v2(i,j,k) + v2(i,j-1,k) ) / dx_sqrd
							+ (mu+lambda) *( v1(i+1,j+1,k) - v1(i-1,j+1,k) - v1(i+1,j-1,k) + v1(i-1,j-1,k) ) / (4.0*dx_sqrd)
							+ (mu+lambda) *( v3(i,j+1,k+1) - v3(i,j+1,k-1) - v3(i,j-1,k+1) + v3(i,j-1,k-1) ) / (4.0*dx_sqrd);

				Lop3(i,j,k) =          mu *( v3(i+1,j,k) + v3(i-1,j,k) + v3(i,j+1,k) + v3(i,j-1,k) + v3(i,j,k+1) + v3(i,j,k-1) - 6.0*v3(i,j,k) ) / dx_sqrd
							+ (mu+lambda) *( v3(i,j,k+1) - 2.0*v3(i,j,k) + v3(i,j,k-1) ) / dx_sqrd
							+ (mu+lambda) *( v1(i+1,j,k+1) - v1(i-1,j,k+1) - v1(i+1,j,k-1) + v1(i-1,j,k-1) ) / (4.0*dx_sqrd)
							+ (mu+lambda) *( v2(i,j+1,k+1) - v2(i,j-1,k+1) - v2(i,j+1,k-1) + v2(i,j-1,k-1) ) / (4.0*dx_sqrd);
			}
		}
	}
}

//#####################################################################

void StationaryNS3D::solveExactly( DoubleArray3D& v1,
								   DoubleArray3D& v2,
								   DoubleArray3D& v3,
								   const DoubleArray3D& f1,
								   const DoubleArray3D& f2,
								   const DoubleArray3D& f3,
								   const double& dx )
{
	double a, b, c, d;
	a = mu / (8.0*mu + 2.0*lambda);
	b = (mu+lambda) / (8.0*mu + 2.0*lambda);
	c = (mu+lambda) / 4.0 / (8.0*mu + 2.0*lambda);
	d = 1.0 / (8.0*mu + 2.0*lambda);

	double dx_sqrd = dx*dx;

	long i = 1;
	long j = 1;
	long k = 1;

	v1(i,j,k) = a * ( v1(i+1,j,k) + v1(i-1,j,k) + v1(i,j+1,k) + v1(i,j-1,k) + v1(i,j,k+1) + v1(i,j,k-1) ) 
			  + b * ( v1(i+1,j,k) + v1(i-1,j,k) )
			  + c * ( v2(i+1,j+1,k) - v2(i-1,j+1,k) - v2(i+1,j-1,k) + v2(i-1,j-1,k) )
			  + c * ( v3(i+1,j,k+1) - v3(i-1,j,k+1) - v3(i+1,j,k-1) + v3(i-1,j,k-1) )
			  - d * ( f1(i,j,k)*dx_sqrd );

	v2(i,j,k) = a * ( v2(i+1,j,k) + v2(i-1,j,k) + v2(i,j+1,k) + v2(i,j-1,k) + v2(i,j,k+1) + v2(i,j,k-1) ) 
			  + b * ( v2(i,j+1,k) + v2(i,j-1,k) )
			  + c * ( v1(i+1,j+1,k) - v1(i-1,j+1,k) - v1(i+1,j-1,k) + v1(i-1,j-1,k) )
			  + c * ( v3(i,j+1,k+1) - v3(i,j+1,k-1) - v3(i,j-1,k+1) + v3(i,j-1,k-1) )
			  - d * ( f2(i,j,k)*dx_sqrd );

	v3(i,j,k) = a * ( v3(i+1,j,k) + v3(i-1,j,k) + v3(i,j+1,k) + v3(i,j-1,k) + v3(i,j,k+1) + v3(i,j,k-1) ) 
			  + b * ( v3(i,j,k+1) + v3(i,j,k-1) )
			  + c * ( v1(i+1,j,k+1) - v1(i-1,j,k+1) - v1(i+1,j,k-1) + v1(i-1,j,k-1) )
			  + c * ( v2(i,j+1,k+1) - v2(i,j-1,k+1) - v2(i,j+1,k-1) + v2(i,j-1,k-1) )
			  - d * ( f3(i,j,k)*dx_sqrd );
}

//#####################################################################

void StationaryNS3D::applyRelaxationGS( DoubleArray3D& v1, 
										DoubleArray3D& v2,
										DoubleArray3D& v3,
										const DoubleArray3D& f1,
										const DoubleArray3D& f2, 
										const DoubleArray3D& f3, 
										const double& dx )
{
	double a, b, c, d;
	a = mu / (8.0*mu + 2.0*lambda);
	b = (mu+lambda) / (8.0*mu + 2.0*lambda);
	c = (mu+lambda) / 4.0 / (8.0*mu + 2.0*lambda);
	d = 1.0 / (8.0*mu + 2.0*lambda);

	// cout << a << " " << b << " " << c << " " << d << endl << endl;

	long m = v1.getIndex1Size();
	long n = v1.getIndex2Size();
	long p = v1.getIndex3Size();

	double dx_sqrd = dx*dx;

	long i,j,k;
	
	//
	// Checkerboard updating:
	//
	long color;
	for( color = 0; color <= 1; color++ )
	{
		for( k=1; k < p-1; k++ )
		{	for( j=1; j < n-1; j++ )	// Do not update the boundary:
			{	for( i=1; i < m-1; i++ )
				{	
					if( (i+j+k)%2 == color )
					{	
						v1(i,j,k) = a * ( v1(i+1,j,k) + v1(i-1,j,k) + v1(i,j+1,k) + v1(i,j-1,k) + v1(i,j,k+1) + v1(i,j,k-1) ) 
								  + b * ( v1(i+1,j,k) + v1(i-1,j,k) )
								  + c * ( v2(i+1,j+1,k) - v2(i-1,j+1,k) - v2(i+1,j-1,k) + v2(i-1,j-1,k) )
								  + c * ( v3(i+1,j,k+1) - v3(i-1,j,k+1) - v3(i+1,j,k-1) + v3(i-1,j,k-1) )
								  - d * ( f1(i,j,k)*dx_sqrd );

						v2(i,j,k) = a * ( v2(i+1,j,k) + v2(i-1,j,k) + v2(i,j+1,k) + v2(i,j-1,k) + v2(i,j,k+1) + v2(i,j,k-1) ) 
								  + b * ( v2(i,j+1,k) + v2(i,j-1,k) )
								  + c * ( v1(i+1,j+1,k) - v1(i-1,j+1,k) - v1(i+1,j-1,k) + v1(i-1,j-1,k) )
								  + c * ( v3(i,j+1,k+1) - v3(i,j+1,k-1) - v3(i,j-1,k+1) + v3(i,j-1,k-1) )
								  - d * ( f2(i,j,k)*dx_sqrd );

						v3(i,j,k) = a * ( v3(i+1,j,k) + v3(i-1,j,k) + v3(i,j+1,k) + v3(i,j-1,k) + v3(i,j,k+1) + v3(i,j,k-1) ) 
								  + b * ( v3(i,j,k+1) + v3(i,j,k-1) )
								  + c * ( v1(i+1,j,k+1) - v1(i-1,j,k+1) - v1(i+1,j,k-1) + v1(i-1,j,k-1) )
								  + c * ( v2(i,j+1,k+1) - v2(i,j-1,k+1) - v2(i,j+1,k-1) + v2(i,j-1,k-1) )
								  - d * ( f3(i,j,k)*dx_sqrd );
					}
				}
			}
		}	
	}

}

//#####################################################################

void StationaryNS3D::findResidual( const DoubleArray3D& v1, 
								   const DoubleArray3D& v2,
								   const DoubleArray3D& v3,
								   const DoubleArray3D& f1,
								   const DoubleArray3D& f2,
								   const DoubleArray3D& f3,
								   DoubleArray3D& R1,
								   DoubleArray3D& R2,
								   DoubleArray3D& R3,
								   const double& dx )
{
	long m = v1.getIndex1Size();
	long n = v1.getIndex2Size();
	long p = v1.getIndex3Size();

	double dx_sqrd = dx*dx;
	
	//
	// find the residual:
	//
	long i,j,k;
	for( i=1; i < m-1; i++ )
	{
		for( j=1; j < n-1; j++ )
		{
			for( k=1; k < p-1; k++ )
			{
				R1(i,j,k) = f1(i,j,k) -          mu *( v1(i+1,j,k) + v1(i-1,j,k) + v1(i,j+1,k) + v1(i,j-1,k) + v1(i,j,k+1) + v1(i,j,k-1) - 6.0*v1(i,j,k) ) / dx_sqrd
									  - (mu+lambda) *( v1(i+1,j,k) - 2.0*v1(i,j,k) + v1(i-1,j,k) ) / dx_sqrd
									  - (mu+lambda) *( v2(i+1,j+1,k) - v2(i-1,j+1,k) - v2(i+1,j-1,k) + v2(i-1,j-1,k) ) / (4.0*dx_sqrd)
									  - (mu+lambda) *( v3(i+1,j,k+1) - v3(i-1,j,k+1) - v3(i+1,j,k-1) + v3(i-1,j,k-1) ) / (4.0*dx_sqrd);
					

				R2(i,j,k) = f2(i,j,k) -          mu *( v2(i+1,j,k) + v2(i-1,j,k) + v2(i,j+1,k) + v2(i,j-1,k) + v2(i,j,k+1) + v2(i,j,k-1) - 6.0*v2(i,j,k) ) / dx_sqrd
									  - (mu+lambda) *( v2(i,j+1,k) - 2.0*v2(i,j,k) + v2(i,j-1,k) ) / dx_sqrd
									  - (mu+lambda) *( v1(i+1,j+1,k) - v1(i-1,j+1,k) - v1(i+1,j-1,k) + v1(i-1,j-1,k) ) / (4.0*dx_sqrd)
									  - (mu+lambda) *( v3(i,j+1,k+1) - v3(i,j+1,k-1) - v3(i,j-1,k+1) + v3(i,j-1,k-1) ) / (4.0*dx_sqrd);

				R3(i,j,k) = f3(i,j,k) -          mu *( v3(i+1,j,k) + v3(i-1,j,k) + v3(i,j+1,k) + v3(i,j-1,k) + v3(i,j,k+1) + v3(i,j,k-1) - 6.0*v3(i,j,k) ) / dx_sqrd
									  - (mu+lambda) *( v3(i,j,k+1) - 2.0*v3(i,j,k) + v3(i,j,k-1) ) / dx_sqrd
									  - (mu+lambda) *( v1(i+1,j,k+1) - v1(i-1,j,k+1) - v1(i+1,j,k-1) + v1(i-1,j,k-1) ) / (4.0*dx_sqrd)
									  - (mu+lambda) *( v2(i,j+1,k+1) - v2(i,j-1,k+1) - v2(i,j+1,k-1) + v2(i,j-1,k-1) ) / (4.0*dx_sqrd);
			}
		}
	}

}
