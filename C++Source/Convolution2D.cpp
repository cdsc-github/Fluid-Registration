//
//#####################################################################
//					 Convolution2D.cpp
//#####################################################################
//
// Igor Yanovsky (C) UCLA
//
//#####################################################################
//

#include <iostream>
using namespace std;

#include <math.h>
#include "DoubleArray2D.h"
#include "Grid2D.h"
#include "FileIO_2D.h"
#include "common_routines.h"

#include <fftw3.h>

#define PI 4.0*atan2(1.0,1.0)

//#####################################################################

DoubleArray2D fftshift2D( const DoubleArray2D& A )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();

	DoubleArray2D B(m,n);

	long i, j;
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			B( (i+m/2)%m, (j+n/2)%n ) = A(i,j);	
			// if m-odd, then m/2 is integer division;
		}
	}

	return B;
}

//#####################################################################

/*
% This function takes the high-resolution MER MI 1D kernel (there are
% more samples than pixels) and radially reconstructs the 2D version of 
% the kernel. The kernel is a function of the radius and not of angle.
% Thus, K2D(x,y) = K1D(r), where r is the radius, i.e. r = sqrt(x^2+y^2).
%
%
% There are different ways of redefining the high-resolution 1D kernel
% onto 2D pixel space:
%     'sample ' option - defines the value of the kernel K2D at
%              pixel (x,y) using sampling, i.e. finding
%              the closest point defined on high-resolution 1D physical 
%              domain to the value of r and assisning the value at that 
%              location to the pixel (x,y).
%     'average' option - defines the value of the kernel K2D at
%              pixel (x,y) by averaging values of K1D defined at those
%              points in the 1D physical domain that are closest to the
%              value r.

% Also see corresponding Matlab code (process_PSF.m) inside \Images\MER\
*/ 

DoubleArray2D MER_MI_kernel( const Grid2D& grid, const char kernel_choice )
{
	// kernel_choice: 'a' - averaging, 's' - sampling
	
	long m = grid.m;
	long n = grid.n;

	long length = 257;

	DoubleArray2D data(3, length);
	
	FileIO_2D IO;
	//IO.readDAT2D( data, "C:\\Igor\\Images\\MER\\MER_MI_PSF_array1_1D.TXT" );		// Windows
	//IO.readDAT2D( data, "/home/yanovsky/MER/images/MER_MI_PSF_array1_1D.TXT" );	// Mipldevlinux
	IO.readDAT2D( data, "/Users/yanovsky/Igor/Images/MER/MER_MI_PSF_array1_1D.TXT" );		// Mac
	
	DoubleArray1D x(length);		// highly-resolved physical/spatial domain
	DoubleArray1D K1D_hres(length);	// highly-resolved value


	long i;
	for( i = 0; i < length; i++ )
	{								// Convert physical domain from [-165,165] to [0,28];
		x(i)        = data(1,i)/12.0;	// a pixel is 12 units wide in a physical domain
		K1D_hres(i) = data(2,i);
	}

	const long d = 14;

	DoubleArray2D Kc(29,29);		// the center of the actual kernel (the actual kernel will
								// be of the size of the image)

	DoubleArray2D number(29,29); // number of high resolution points in a pixel

	double r;
	long xx, yy;

	double min_distance_from_pixel = 10000.0;

	if( kernel_choice == 'a' )
	{
		for( xx = -d; xx <= d; xx++ )
		{
			for( yy = -d; yy <= d; yy++ )
			{
				r = sqrt(double(xx*xx + yy*yy));

				if( r <= double(d) )	// if r > 14, K = 0
				{
					for( i = 0; i < length; i++ )
					{
						if( fabs( x(i) - r  ) < 0.5 )
						{
							Kc(xx+d,yy+d) = Kc(xx+d,yy+d) + K1D_hres(i);
							number(xx+d,yy+d) = number(xx+d,yy+d) + 1;
						}
					}
				}
				else
				{
					number(xx+d,yy+d) = 1;	// to avoid later devision by 0
				}
			}
		}
		for( xx = -d; xx <= d; xx++ )
		{
			for( yy = -d; yy <= d; yy++ )
			{
				Kc(xx+d,yy+d) = Kc(xx+d,yy+d)/number(xx+d,yy+d);
			}
		}
		IO.writePGMsc( Kc, "Kc_average" );
	}

	else if( kernel_choice == 's' )
	{
		for( xx = -d; xx <= d; xx++ )
		{
			for( yy = -d; yy <= d; yy++ )
			{
				min_distance_from_pixel = 10000.0;
				r = sqrt(double(xx*xx + yy*yy));
			
				if( r <= double(d) )		// if r > 14, K = 0
				{
					for( i = 0; i < length; i++ )
					{
						if( fabs( x(i) - r  ) < min_distance_from_pixel )
						{
							Kc(xx+d,yy+d) = K1D_hres(i);	// set K1D to value which is closest to the pixel
							min_distance_from_pixel = fabs( x(i) - r  );
						}
					}
				}
			}
		}
		IO.writePGMsc( Kc, "Kc_sample" );
	}

	IO.write_bin_float( Kc, "Kc" );


	long a = (m - 29)/2;
	long b = (n - 29)/2;
	cout << a << endl << b << endl;

	DoubleArray2D K(m,n);

	long j;
	for(j = 0; j < 29; j++)
	{	
		for(i = 0; i < 29; i++)
		{
			K(i+a,j+b) = Kc(i,j);
		}
	}
	
	double TotalWeight = 0.0;
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	TotalWeight = TotalWeight + K(i,j);	}	}

	findMaxMin( K );

	cout << "Total weight of the PSF = " << TotalWeight << endl << endl;
	
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	K(i,j) = K(i,j)/TotalWeight;	}	}

	TotalWeight = 0.0;
	
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	TotalWeight = TotalWeight + K(i,j);	}	}

	findMaxMin( K );

	cout << "Total weight of the PSF = " << TotalWeight << endl << endl;
	
	return K;
}

//#####################################################################

DoubleArray2D GeoSTAR( const Grid2D& grid )
{
	FileIO_2D IO;
	long m = grid.m;
	long n = grid.n;
	
	DoubleArray2D Kc(101,101);
	IO.readDAT2D( Kc, "C:\\Igor\\Images\\GeoSTAR\\afk.txt" );

	DoubleArray2D K(m,n);

	long i, j;
	for(j = 0; j < 101; j++)
	{	
		for(i = 0; i < 101; i++)
		{
			K(i+77,j+77) = Kc(i,j);
		}
	}

	/*
	if( m == 256 && n == 256 )
	{	IO.readDAT2D( K, "C:\\Igor\\Images\\GeoSTAR\\afk256.txt" );		}
	else if( m == 258 && n == 258 )
	{	IO.readDAT2D( K, "C:\\Igor\\Images\\GeoSTAR\\afk258.txt" );		}
	else
	{	cout << "Generate aftXXX.txt file for appropriate dimension XXX!" << endl;
		exit(1);	}
	*/
	
	double TotalWeight = 0.0;
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	TotalWeight = TotalWeight + K(i,j);	}	}

	findMaxMin( K );
	
	cout << "Total weight of the PSF = " << TotalWeight << endl << endl;
	
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	K(i,j) = K(i,j)/TotalWeight;	}	}

	TotalWeight = 0.0;
	
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{	TotalWeight = TotalWeight + K(i,j);	}	}

	findMaxMin( K );

	cout << "Total weight of the PSF = " << TotalWeight << endl << endl;
	
	
	return K;
}

//#####################################################################

DoubleArray2D OutOfFocus2D( const Grid2D& grid, const double& radius )
{
	long m = grid.m;
	long n = grid.n;

	double xcent = (m-1.0)/2.0;
	double ycent = (n-1.0)/2.0;
	
	double TotalWeight = 0.0;
	
	DoubleArray2D K(m,n);

	long i, j;
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{
			//cout << i << " " << j << " " << sqrt( (double(i)-xcent)*(double(i)-xcent)+(double(j)-ycent)*(double(j)-ycent) ) << endl;
			if( sqrt( (double(i)-xcent)*(double(i)-xcent)+(double(j)-ycent)*(double(j)-ycent) ) < radius )
			{	K(i,j) = 1.0;
			}
		}
	}

	//cout << K << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K(i,j) = K(i,j)/TotalWeight;
		}
	}

	return K;
}

//#####################################################################

DoubleArray2D OutOfFocus2D_directional( const Grid2D& grid, const double& rx, const double& ry )
{
	// Out-of-Focus blur, having ellipsoid shape

	long m = grid.m;
	long n = grid.n;

	double xcent = (m-1.0)/2.0;
	double ycent = (n-1.0)/2.0;
	
	double TotalWeight = 0.0;
	
	DoubleArray2D K(m,n);

	long i, j;
	for( i = 0; i < m; i++ )
	{	for( j = 0; j < n; j++ )
		{
			if( sqrt( (double(i)-xcent)*(double(i)-xcent)/(rx*rx)+(double(j)-ycent)*(double(j)-ycent)/(ry*ry) ) < 1.0 )
			{	K(i,j) = 1.0;
			}
		}
	}

	//cout << K << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K(i,j) = K(i,j)/TotalWeight;
		}
	}

	return K;
}

//#####################################################################

DoubleArray2D RectBlur( const Grid2D& grid, const double& width, const double& height )
{
	// Out-of-Focus blur, having rectangular shape

	long m = grid.m;
	long n = grid.n;

	double xcent, ycent;
	if(m%2 == 0)	{	xcent = (m-1.0)/2.0;	}
	else			{	xcent = (m+1.0)/2.0;	}
	if(n%2 == 0)	{	ycent = (n-1.0)/2.0;	}
	else			{	ycent = (n+1.0)/2.0;	}
	
	double TotalWeight = 0.0;
	
	DoubleArray2D K(m,n);

	long i, j;
	for(j = 0; j < n; j++)		// Square
	{	
		for(i = 0; i < m; i++)
		{	
			if( fabs(double(i)-xcent) < (width)/2.0 && fabs(double(j)-ycent) < (height)/2.0 )
			{
				K(i,j) = 1.0;
			}
		}
	}

	//cout << K << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K(i,j) = K(i,j)/TotalWeight;
		}
	}

	return K;
}

//#####################################################################

DoubleArray2D ComplexShapeBlur( const Grid2D& grid, const double& w1, const double& h1, const double& w2, const double& h2 )
{
	long m = grid.m;
	long n = grid.n;

	long xcent, ycent;

	//if(m%2 == 0)	{	xcent = (m-1)/2;	}
	if(m%2 == 0)	{	xcent = (m)/2;	}
	else			{	xcent = (m+1)/2;	}
	//if(n%2 == 0)	{	ycent = (n-1)/2;	}
	if(n%2 == 0)	{	ycent = (n)/2;	}
	else			{	ycent = (n+1)/2;	}
	
	double TotalWeight = 0.0;
	
	DoubleArray2D delta(m,n);
	DoubleArray2D K1(m,n);
	DoubleArray2D K2(m,n);
	DoubleArray2D K_total(m,n);

	delta(xcent,ycent) = 1.0;

	long i, j;
	for(j = 0; j < n; j++)		// Square 1
	{	for(i = 0; i < m; i++)
		{	
			if( abs(i-xcent) < (w1)/2.0 && abs(j-ycent) < (h1)/2.0 )
			{
				K1(i,j) = 1.0;
			}
		}
	}
	
	for(j = 0; j < n; j++)		// Square 2
	{	
		for(i = 0; i < m; i++)
		{	
			if( abs(i-xcent) < w2/2.0 && abs(j-ycent) < h2/2.0 )
			{
				K2(i,j) = 1.0;
			}
		}
	}

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + delta(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			delta(i,j) = delta(i,j)/TotalWeight;
		}
	}

	TotalWeight = 0.0;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K1(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K1(i,j) = K1(i,j)/TotalWeight;
		}
	}

	TotalWeight = 0.0;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K2(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K2(i,j) = K2(i,j)/TotalWeight;
		}
	}



	cout << endl << delta(320,512) << " " << K1(320,512) << " " << K2(320,512) << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K_total(i,j) = delta(i,j) + K2(i,j) - K1(i,j);
			//if( K_total(i,j) < 0.0 )
			//{	cout<< "Exiting..." << endl;  exit(1);}
		}
	}	

	TotalWeight = 0.0;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K_total(i,j);
		}
	}

	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K_total(i,j) = K_total(i,j)/TotalWeight;
		}
	}

	return K_total;
}


//#####################################################################

DoubleArray2D Gaussian2D( const Grid2D& grid, const double& std )
{
	// NOTE: Center is defined to be at (m/2,n/2).
	// This produces no shift of the blurred image wrt original image.
	long   m    = grid.m;
	long   n    = grid.n;
	double dx   = grid.dx;
	double dy   = grid.dy;
	double xMin = grid.xMin;
	double xMax = grid.xMax;
	double yMin = grid.yMin;
	double yMax = grid.yMax;

	double xcent, ycent;
	//if(m%2 == 0)	{	xcent =  m/2.0;		}
	if(m%2 == 0)	{	xcent = (m-1.0)/2.0;		}
	else			{	xcent = (m+1.0)/2.0;}
	//if(n%2 == 0)	{	ycent =  n/2.0;		}
	if(n%2 == 0)	{	ycent = (n-1.0)/2.0;		}
	else			{	ycent = (n+1.0)/2.0;}
	
	double TotalWeight = 0.0;

	double alpha = 1.0/(2.0*PI*std*std);
	double beta = 2.0*std*std;

	double xx, yy;
	
	DoubleArray2D K(m,n);
	
	long i, j;
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			xx = double(i)*dx + xMin;
			yy = double(j)*dy + yMin;
			K(i,j) = alpha * exp( -((xx-xcent)*(xx-xcent)+(yy-ycent)*(yy-ycent)) / beta );
		}
	}

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			TotalWeight = TotalWeight + K(i,j);
		}
	}
	cout << TotalWeight << endl << endl;

	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
			K(i,j) = K(i,j)/TotalWeight;
		}
	}

	cout << "Kernel K: ";
	findMaxMin( K );

	return K;
}

//#####################################################################

DoubleArray2D Gaussian2D( const long m, const double std)
{
	// NOTE: Center is defined to be at (m/2,n/2).
	// This produces no shift of the blurred image wrt original image.

	double cent;
	if(m%2 == 0)	{	cent =  m/2.0;		}
	else			{	cent = (m+1.0)/2.0;	}
	
	double TotalWeight = 0.0;

	double alpha = 1.0/(2.0*PI*std*std);
	double beta = 2.0*std*std;

	double ii, jj;
	
	DoubleArray2D K(m,m);
	
	long i1, i2;
	for( i1 = 0; i1 < m; i1++ )
	{
		for( i2 = 0; i2 < m; i2++ )
		{
			ii = double(i1);
			jj = double(i2);
			K(i1,i2) = alpha * exp( -((ii-cent)*(ii-cent)+(jj-cent)*(jj-cent)) / beta );
		}
	}

	for( i1 = 0; i1 < m; i1++ )
	{
		for( i2 = 0; i2 < m; i2++ )
		{
			TotalWeight = TotalWeight + K(i1,i2);
		}
	}
	cout << TotalWeight << endl << endl;

	for( i1 = 0; i1 < m; i1++ )
	{
		for( i2 = 0; i2 < m; i2++ )
		{
			K(i1,i2) = K(i1,i2)/TotalWeight;
		}
	}

	return K;
}

//#####################################################################

DoubleArray2D dGdy( const long I, const double std)
{
	// Derivative of the Gaussian G(x,y) with respect to y.
	
	// NOTE: Center is defined to be at (m/2,n/2).
	// This produces no shift of the blurred image wrt original image.
	
	double cent;
	if(I%2 == 0)	{	cent =  I/2.0;		}
	else			{	cent = (I+1.0)/2.0;}

	double alpha = 1.0/(2.0*PI*pow(std,4.0));
	double beta = 2.0*std*std;

	double ii, jj;
	
	DoubleArray2D K(I,I);
	
	long i1, i2;
	for( i1 = 0; i1 < I; i1++ )
	{
		for( i2 = 0; i2 < I; i2++ )
		{
			ii = double(i1);
			jj = double(i2);
			K(i1,i2) = -(jj-cent) * alpha * exp( -((ii-cent)*(ii-cent)+(jj-cent)*(jj-cent)) / beta );
		}
	}

	return K;
}

//#####################################################################

void convolution2D( DoubleArray2D& F, const DoubleArray2D& K )
{
	// Input: F is the clean image, K is the PSF
	// Output: F is the convolution (F = K*F)

	long m = K.getIndex1Size();
	long n = K.getIndex2Size();

	fftw_complex *inK, *outK;
	fftw_complex *inF, *outF;
	fftw_complex *fftK_fftF;
	fftw_complex *ifft_fftK_fftF;
	fftw_plan p1, p2, p3;

	inK  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outK = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	inF  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	fftK_fftF      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	ifft_fftK_fftF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	long i, j;
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	inK[i*n + j][0] = K(i,j);
			inK[i*n + j][1] = 0.0;
			inF[i*n + j][0] = F(i,j);
			inF[i*n + j][1] = 0.0;
		}
	}

	p1 = fftw_plan_dft_2d( m, n, inK, outK, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p1);

	p2 = fftw_plan_dft_2d( m, n, inF, outF, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p2);
	
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )	// complex multiplication:  fft(K)*fft(F)
		{	fftK_fftF[i*n + j][0] = outK[i*n + j][0] * outF[i*n + j][0] - outK[i*n + j][1] * outF[i*n + j][1];
			fftK_fftF[i*n + j][1] = outK[i*n + j][0] * outF[i*n + j][1] + outK[i*n + j][1] * outF[i*n + j][0];
		}
	}

	p3 = fftw_plan_dft_2d( m, n, fftK_fftF, ifft_fftK_fftF, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute(p3);

	// ifft( fft(K)*fft(F) ):
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	F(i,j) = ifft_fftK_fftF[i*n + j][0]/m/n;	// F is returned back
			//ifft_fftK_fftF[i*n + j][1] = ifft_fftK_fftF[i*n + j][1]/m/n;
		}
	}

	fftw_destroy_plan(p1); fftw_destroy_plan(p2); fftw_destroy_plan(p3);
	fftw_free(inK); fftw_free(outK);
	fftw_free(inF); fftw_free(outF);
	fftw_free(fftK_fftF);
	fftw_free(ifft_fftK_fftF);
}

//#####################################################################

void deconvFidelity2D( DoubleArray2D& result, const DoubleArray2D& F, const DoubleArray2D& K, const DoubleArray2D& Krot, const DoubleArray2D& U )
{
	long m = K.getIndex1Size();
	long n = K.getIndex2Size();

	fftw_complex *inK,    *outK;
	fftw_complex *inKrot, *outKrot;
	fftw_complex *inF,    *outF;
	fftw_complex *inU,    *outU;

	fftw_complex *fftKrot_fftF;						// fft(Kr)*fft(F)
	fftw_complex *fftKrot_fftK;						// fft(Kr)*fft(K)
	fftw_complex *fftKrot_fftK_fftU;				// fft(Kr)*fft(K)*fft(u)
	fftw_complex *fftKrot_fftK_fftU__fftKrot_fftF;	// fft(Kr)*fft(K)*fft(u) - fft(Kr)*fft(F)

	fftw_complex *ifft_result;
	fftw_plan p1, p2, p3, p4, p5;

	inK     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outK    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	inKrot  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outKrot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	inF     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outF    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	inU     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	outU    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	fftKrot_fftF                    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	fftKrot_fftK                    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	fftKrot_fftK_fftU               = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	fftKrot_fftK_fftU__fftKrot_fftF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);
	ifft_result                     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m*n);

	long i, j;
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	inKrot[i*n + j][0] = Krot(i,j);
			inKrot[i*n + j][1] = 0.0;
			inK[i*n + j][0] = K(i,j);
			inK[i*n + j][1] = 0.0;

			inF[i*n + j][0] = F(i,j);
			inF[i*n + j][1] = 0.0;
			inU[i*n + j][0] = U(i,j);
			inU[i*n + j][1] = 0.0;
		}
	}

	p1 = fftw_plan_dft_2d( m, n, inK, outK, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p1);

	p2 = fftw_plan_dft_2d( m, n, inF, outF, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p2);

	p3 = fftw_plan_dft_2d( m, n, inU, outU, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p3);

	p4 = fftw_plan_dft_2d( m, n, inKrot, outKrot, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p4);
	
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	// complex multiplication:  fft(Kr)*fft(F)     NOT fft(K)*fft(F)
			fftKrot_fftF[i*n + j][0] = outKrot[i*n + j][0] * outF[i*n + j][0] - outKrot[i*n + j][1] * outF[i*n + j][1];
			fftKrot_fftF[i*n + j][1] = outKrot[i*n + j][0] * outF[i*n + j][1] + outKrot[i*n + j][1] * outF[i*n + j][0];
			
			// fft(Kr)*fft(K)
			fftKrot_fftK[i*n + j][0] = outKrot[i*n + j][0] * outK[i*n + j][0] - outKrot[i*n + j][1] * outK[i*n + j][1];
			fftKrot_fftK[i*n + j][1] = outKrot[i*n + j][0] * outK[i*n + j][1] + outKrot[i*n + j][1] * outK[i*n + j][0];
			
			// fft(Kr)*fft(K)*fft(u)
			fftKrot_fftK_fftU[i*n + j][0] = fftKrot_fftK[i*n + j][0] * outU[i*n + j][0] - fftKrot_fftK[i*n + j][1] * outU[i*n + j][1];
			fftKrot_fftK_fftU[i*n + j][1] = fftKrot_fftK[i*n + j][0] * outU[i*n + j][1] + fftKrot_fftK[i*n + j][1] * outU[i*n + j][0];

			// fft(Kr)*fft(K)*fft(u) - fft(Kr)*fft(F)
			fftKrot_fftK_fftU__fftKrot_fftF[i*n + j][0] = fftKrot_fftK_fftU[i*n + j][0] - fftKrot_fftF[i*n + j][0];
			fftKrot_fftK_fftU__fftKrot_fftF[i*n + j][1] = fftKrot_fftK_fftU[i*n + j][1] - fftKrot_fftF[i*n + j][1];
		}
	}

	p5 = fftw_plan_dft_2d( m, n, fftKrot_fftK_fftU__fftKrot_fftF, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute(p5);

	// ifft( fft(Kr)*fft(K)*fft(u) - fft(Kr)*fft(F) ):
	for(i = 0; i < m; i++)
	{	for( j = 0; j < n; j++ )
		{	result(i,j) = ifft_result[i*n + j][0]/m/n;	// F is returned back
			//ifft_result[i*n + j][1] = ifft_result[i*n + j][1]/m/n;
		}
	}

	fftw_destroy_plan(p1); fftw_destroy_plan(p2); fftw_destroy_plan(p3); fftw_destroy_plan(p4); fftw_destroy_plan(p5);
	fftw_free(inK); fftw_free(outK);
	fftw_free(inF); fftw_free(outF);
	fftw_free(inU); fftw_free(outU);
	fftw_free(inKrot); fftw_free(outKrot);

	fftw_free(fftKrot_fftF);
	fftw_free(fftKrot_fftK);
	fftw_free(fftKrot_fftK_fftU);
	fftw_free(fftKrot_fftK_fftU__fftKrot_fftF);

	fftw_free(ifft_result);
}

//#####################################################################

void convGaussian2D( DoubleArray2D& F, const Grid2D& grid, const double& std )
{
	long m = F.getIndex1Size();
	long n = F.getIndex2Size();

	DoubleArray2D K(m,n);		// define Kernel K
	//K = OutOfFocus2D( grid, 2.0 );
	K = Gaussian2D( grid, std );

	K = fftshift2D(K);			// FFTSHIFT

	//FileIO_2D IO;

	convolution2D( F, K );

	//IO.write_ascii( K, "kernel" );
}
