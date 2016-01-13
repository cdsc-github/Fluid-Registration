//
//#####################################################################
//					 Convolution3D.cpp
//#####################################################################
//
// Igor Yanovsky (C) UCLA
//
//#####################################################################
//
#include <iostream>
using namespace std;

#include <math.h>
#include "DoubleArray3D.h"
#include "Grid3D.h"
//#include "grid3d.h"

#include <fftw3.h>

#define PI 4.0*atan2(1.0,1.0)

//#####################################################################

DoubleArray3D fftshift3D( const DoubleArray3D& A )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	DoubleArray3D B(m,n,p);

	long i,j,k;
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			for( k = 0; k < p; k++ )
			{
				B( (i+m/2)%m, (j+n/2)%n, (k+p/2)%p ) = A(i,j,k);	
				// if m-odd, then m/2 is integer division;
			}
		}
	}

	return B;
}

//#####################################################################

DoubleArray3D Gaussian3D( const Grid3D& grid, const double std)
{
	// NOTE: Center is defined to be at (m/2,n/2,p/2).
	// This produces no shift of the blurred image wrt original image.
	long   m    = grid.m;
	long   n    = grid.n;
	long   p    = grid.p;
	double dx   = grid.dx;
	double dy   = grid.dy;
	double dz   = grid.dz;
	double xMin = grid.xMin;
	double xMax = grid.xMax;
	double yMin = grid.yMin;
	double yMax = grid.yMax;
	double zMin = grid.zMin;
	double zMax = grid.zMax;

	double xcent, ycent, zcent;
	if(m%2 == 0)	{	xcent =  m/2.0;		}
	else			{	xcent = (m+1.0)/2.0;}
	if(n%2 == 0)	{	ycent =  n/2.0;		}
	else			{	ycent = (n+1.0)/2.0;}
	if(p%2 == 0)	{	zcent =  p/2.0;		}
	else			{	zcent = (p+1.0)/2.0;}
	
	double TotalWeight = 0.0;

	double alpha = 1.0 / pow( 2.0*PI*std*std, 1.5 );
	double beta = 2.0*std*std;

	double xx, yy, zz;
	
	DoubleArray3D K(m,n,p);
	
	long i,j,k;
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			for( k = 0; k < p; k++ )
			{
				xx = double(i)*dx + xMin;
				yy = double(j)*dy + yMin;
				zz = double(k)*dz + zMin;
				K(i,j,k) = alpha * exp( -((xx-xcent)*(xx-xcent)+(yy-ycent)*(yy-ycent)+(zz-zcent)*(zz-zcent)) / beta );
			}
		}
	}

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			for( k = 0; k < p; k++ )
			{
				TotalWeight = TotalWeight + K(i,j,k);
			}
		}
	}
	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			for( k = 0; k < p; k++ )
			{
				K(i,j,k) = K(i,j,k)/TotalWeight;
			}
		}
	}

	return K;
}

//#####################################################################

void convolution3D( DoubleArray3D& F, const DoubleArray3D& K )
{
	long m = K.getIndex1Size();
	long n = K.getIndex2Size();
	long p = K.getIndex3Size();

	fftw_complex *inK, *outK;
	fftw_complex *inF, *outF;
	fftw_complex *fftK_fftF;
	fftw_complex *ifft_fftK_fftF;
	fftw_plan p1, p2, p3;

	inK  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);
	outK = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);

	inF  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);
	outF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);

	fftK_fftF      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);
	ifft_fftK_fftF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n*p);

	long i,j,k, ind;
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	for( k = 0; k < p; k++ )
			{	ind = k + p*(j+i*n);
				inK[ind][0] = K(i,j,k);
				inK[ind][1] = 0.0;
				inF[ind][0] = F(i,j,k);
				inF[ind][1] = 0.0;
			}
		}
	}

	p1 = fftw_plan_dft_3d( m, n, p, inK, outK, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p1);

	p2 = fftw_plan_dft_3d( m, n, p, inF, outF, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p2);
	
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )	// complex multiplication:  fft(K)*fft(F)
		{	for( k = 0; k < p; k++ )
			{	ind = k + p*(j+i*n);
				fftK_fftF[ind][0] = outK[ind][0] * outF[ind][0] - outK[ind][1] * outF[ind][1];
				fftK_fftF[ind][1] = outK[ind][0] * outF[ind][1] + outK[ind][1] * outF[ind][0];
			}
		}
	}

	p3 = fftw_plan_dft_3d( m, n, p, fftK_fftF, ifft_fftK_fftF, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute(p3);

	// ifft( fft(K)*fft(F) ):
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	for( k = 0; k < p; k++ )
			{	ind = k + p*(j+i*n);
				F(i,j,k) = ifft_fftK_fftF[ind][0]/m/n/p;	// F is returned back
			}
		}
	}

	fftw_destroy_plan(p1); fftw_destroy_plan(p2); fftw_destroy_plan(p3);
	fftw_free(inK); fftw_free(outK);
	fftw_free(inF); fftw_free(outF);
	fftw_free(fftK_fftF);
	fftw_free(ifft_fftK_fftF);
}

//#####################################################################

void convGaussian3D( DoubleArray3D& F, const Grid3D& grid, const double& std )
{
	long m = F.getIndex1Size();
	long n = F.getIndex2Size();
	long p = F.getIndex3Size();

	DoubleArray3D K(m,n,p);		// define Kernel K
	K = Gaussian3D( grid, std );

	K = fftshift3D(K);			// FFTSHIFT

	convolution3D( F, K );
}
