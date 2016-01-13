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
#include "DoubleArray1D.h"

#include "J_Reg3D.h"
#include "BC_3D.h"
#include "Norms.h"
#include "FileIO_3D.h"

#include "Convolution2D.h"
#include "Convolution3D.h"
#include "MGsystem3D.h"
#include "StationaryNS3D.h"


#define PI 4.0*atan2(1.0,1.0)

//#####################################################################

J_Reg3D::J_Reg3D()
{
	Jmin = 10.0;
	Jmax = 0.0;

	SSD = 0.0;			// L2-norm value
	MI  = 0.0;			// Mutual Information value

	BCtype = 5;

	totalTime = 0.0;

	triangle_u = 0.1;	// max displacement allowed for adaptive timestepping
	fidelity = 1;		// 1 = L2, 2 = MI

	smoothing = 1;		// 1 for Gaussian, 2 for Heat eqn, 3 for NS
	StD = 4.5;
	Parzen_StD = 6.0;
	heat_eqn_TimeSteps   = 100;
}

//#####################################################################

void J_Reg3D::calculateEnergy_L2( const Grid3D& grid, double& energy_L2 )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;

	energy_L2    = SSD /(m-2)/(n-2)/(p-2);
}

//#####################################################################

void J_Reg3D::calculateEnergy_MI( const Grid3D& grid, double& energy_MI )
{
	long m  = grid.m;
	long n  = grid.n;
	long p =  grid.p;

	energy_MI    = -MI /m/n/p;
}

//#####################################################################

DoubleArray3D J_Reg3D::linearInterpolation( const DoubleArray3D& T, 
											const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
											const Grid3D& grid )
{
	long   m    = grid.m;
	long   n    = grid.n;
	long   p    = grid.p;
	long   w    = grid.w;

	double dx   = grid.dx;
	double dy   = grid.dy;
	double dz   = grid.dz;

	double xMin = grid.xMin;
	double xMax = grid.xMax;
	double yMin = grid.yMin;
	double yMax = grid.yMax;
	double zMin = grid.zMin;
	double zMax = grid.zMax;

	DoubleArray3D Tx      = T;		// interpolation in x-direction
	DoubleArray3D Ty      = T;
	DoubleArray3D interpT = T;

	DoubleArray3D X(m,n,p), Y(m,n,p), Z(m,n,p);
	long i,j,k;
	for( k = 1; k < p-1; k++ )
	{
		for( j = 1; j < n-1; j++ )
		{
			for( i = 1; i < m-1; i++ )
			{
				X(i,j,k) = i - u1(i,j,k);		// X = x-u1
				Y(i,j,k) = j - u2(i,j,k);		// Y = y-u2
				Z(i,j,k) = k - u3(i,j,k);
			}
		}
	}

	long newX, newY, newZ;
	for( k = 1; k < p-1; k++ )
	{
		for( j = 1; j < n-1; j++ )
		{
			for( i = 1; i < m-1; i++ )
			{
				newX = long(floor(X(i,j,k)));
				Tx(i,j,k)      = T(newX,j,k) + ((T(newX+1,j,k)-T(newX,j,k))/dx) * (X(i,j,k)-newX);		// T(x-u1)
			}
		}
	}
	
	for( k = 1; k < p-1; k++ )
	{
		for( j = 1; j < n-1; j++ )
		{
			for( i = 1; i < m-1; i++ )
			{
				newY = long(floor(Y(i,j,k)));
				Ty(i,j,k)      = Tx(i,newY,k) + ((Tx(i,newY+1,k)-Tx(i,newY,k))/dy) * (Y(i,j,k)-newY);		// T(y-u2)
			}
		}
	}

	for( k = 1; k < p-1; k++ )
	{
		for( j = 1; j < n-1; j++ )
		{
			for( i = 1; i < m-1; i++ )
			{
				newZ = long(floor(Z(i,j,k)));
				interpT(i,j,k) = Ty(i,j,newZ) + ((Ty(i,j,newZ+1)-Ty(i,j,newZ))/dz) * (Z(i,j,k)-newZ);		// T(y-u2)
			}
		}
	}

	applyBC( interpT, w, BCtype );

	return interpT;
}

//#####################################################################

void J_Reg3D::update_U( DoubleArray3D& u1, DoubleArray3D& u2, DoubleArray3D& u3,
						const DoubleArray3D& v1, const DoubleArray3D& v2, const DoubleArray3D& v3,
						const Grid3D& grid )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;
	long w  = grid.w;

	double dx = grid.dx;
	double dy = grid.dy;
	double dz = grid.dz;

	double dx2 = dx*2.0;
	double dy2 = dy*2.0;
	double dz2 = dz*2.0;
	
	double du1_dx, du2_dx, du3_dx;
	double du1_dy, du2_dy, du3_dy;
	double du1_dz, du2_dz, du3_dz;

	DoubleArray3D R1(m,n,p), R2(m,n,p), R3(m,n,p);
	double Rmax_L2 = 0.0, R_L2 = 0.0;

	long i,j,k;
	for( i = w; i < m-w; i++ )
	{
		for( j = w; j < n-w; j++ )
		{
			for( k = w; k < p-w; k++ )
			{
				du1_dx = (u1(i+1,j,k) - u1(i-1,j,k)) / dx2;
				du2_dx = (u2(i+1,j,k) - u2(i-1,j,k)) / dx2;
				du3_dx = (u3(i+1,j,k) - u3(i-1,j,k)) / dx2;

				du1_dy = (u1(i,j+1,k) - u1(i,j-1,k)) / dy2;
				du2_dy = (u2(i,j+1,k) - u2(i,j-1,k)) / dy2;
				du3_dy = (u3(i,j+1,k) - u3(i,j-1,k)) / dy2;

				du1_dz = (u1(i,j,k+1) - u1(i,j,k-1)) / dz2;
				du2_dz = (u2(i,j,k+1) - u2(i,j,k-1)) / dz2;
				du3_dz = (u3(i,j,k+1) - u3(i,j,k-1)) / dz2;

				R1(i,j,k) = v1(i,j,k) - v1(i,j,k)*du1_dx - v2(i,j,k)*du1_dy - v3(i,j,k)*du1_dz;
				R2(i,j,k) = v2(i,j,k) - v1(i,j,k)*du2_dx - v2(i,j,k)*du2_dy - v3(i,j,k)*du2_dz;
				R3(i,j,k) = v3(i,j,k) - v1(i,j,k)*du3_dx - v2(i,j,k)*du3_dy - v3(i,j,k)*du3_dz;

				R_L2 = sqrt( R1(i,j,k)*R1(i,j,k) + R2(i,j,k)*R2(i,j,k) + R3(i,j,k)*R3(i,j,k) );
				Rmax_L2 = (Rmax_L2 > R_L2 ) ? Rmax_L2 : R_L2;
			}
		}
	}

	double dt;
	dt = triangle_u/Rmax_L2;
	//dt = 0.0001;

	u1 = u1 + R1*dt;
	u2 = u2 + R2*dt;
	u3 = u3 + R3*dt;

	totalTime += dt;
	cout << " dt = " << dt << "  time = " << totalTime;

	applyBC( u1, w, BCtype );
	applyBC( u2, w, BCtype );
	applyBC( u3, w, BCtype );
}

//#####################################################################

void J_Reg3D::evaluate_f_L2( const DoubleArray3D& T, const DoubleArray3D& S,
							 const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
							 const Grid3D& grid,
							 DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3 )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;
	long w  = grid.w;

	double dx = grid.dx;
	double dy = grid.dy;
	double dz = grid.dz;

	double dx2 = dx*2.0;
	double dy2 = dy*2.0;
	double dz2 = dz*2.0;

	double dTdx, dTdy, dTdz;

	SSD = 0.0;

	long i, j, k;
	for( i = 1; i < m-1; i++ )
	{	for( j = 1; j < n-1; j++ )		// unscaled SSD;
		{	for( k = 1; k < p-1; k++ )	// Calculated at beginning of each iteration
	{	SSD += 0.5*(T(i,j,k)-S(i,j,k))*(T(i,j,k)-S(i,j,k));	}	}	}

	cout << "SSD = " << SSD /(m-2.0)/(n-2.0)/(p-2.0) << "  ";
																		
	for( i = 1; i < m-1; i++ )
	{
		for( j = 1; j < n-1; j++ )
		{
			for( k = 1; k < p-1; k++ )
			{
				dTdx = (T(i+1,j,k) - T(i-1,j,k)) / dx2;
				dTdy = (T(i,j+1,k) - T(i,j-1,k)) / dy2;
				dTdz = (T(i,j,k+1) - T(i,j,k-1)) / dz2;

				f1(i,j,k) = -( T(i,j,k) - S(i,j,k) ) * dTdx;
				f2(i,j,k) = -( T(i,j,k) - S(i,j,k) ) * dTdy;
				f3(i,j,k) = -( T(i,j,k) - S(i,j,k) ) * dTdz;
			}
		}
	}


	//////////////////////////////////////////////////

	DoubleArray3D d1_d1(m,n,p), d1_d2(m,n,p), d1_d3(m,n,p);
	DoubleArray3D d2_d1(m,n,p), d2_d2(m,n,p), d2_d3(m,n,p);
	DoubleArray3D d3_d1(m,n,p), d3_d2(m,n,p), d3_d3(m,n,p);

	Jmin = 10.0;
	Jmax = 0.0;

	for ( i = 1; i < m-1; i++ ) {
		for ( j = 1; j < n-1; j++ ) {
			for ( k = 1; k < p-1; k++ ) {
				d1_d1(i,j,k) = 1.0 - (u1(i+1,j,k)-u1(i-1,j,k))/2.0;
				d1_d2(i,j,k) =	   - (u1(i,j+1,k)-u1(i,j-1,k))/2.0;
				d1_d3(i,j,k) =	   - (u1(i,j,k+1)-u1(i,j,k-1))/2.0;

				d2_d1(i,j,k) =	   - (u2(i+1,j,k)-u2(i-1,j,k))/2.0;
				d2_d2(i,j,k) = 1.0 - (u2(i,j+1,k)-u2(i,j-1,k))/2.0;
				d2_d3(i,j,k) =     - (u2(i,j,k+1)-u2(i,j,k-1))/2.0;

				d3_d1(i,j,k) =	   - (u3(i+1,j,k)-u3(i-1,j,k))/2.0;
				d3_d2(i,j,k) =     - (u3(i,j+1,k)-u3(i,j-1,k))/2.0;
				d3_d3(i,j,k) = 1.0 - (u3(i,j,k+1)-u3(i,j,k-1))/2.0;

				J(i,j,k) =  d1_d1(i,j,k)*d2_d2(i,j,k)*d3_d3(i,j,k)
					      + d2_d1(i,j,k)*d3_d2(i,j,k)*d1_d3(i,j,k)
						  + d1_d2(i,j,k)*d2_d3(i,j,k)*d3_d1(i,j,k)
						  - d1_d3(i,j,k)*d2_d2(i,j,k)*d3_d1(i,j,k)
						  - d1_d2(i,j,k)*d2_d1(i,j,k)*d3_d3(i,j,k)
						  - d2_d3(i,j,k)*d3_d2(i,j,k)*d1_d1(i,j,k);

				if (J(i,j,k) < Jmin)  Jmin = J(i,j,k);
				if (J(i,j,k) > Jmax)  Jmax = J(i,j,k);
			}
		}
	}
	applyBC( J, 1, BCtype );	// for energy/StD calculation
}

//#####################################################################

void J_Reg3D::evaluate_f_MI( const DoubleArray3D& T, const DoubleArray3D& S,
							 const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
							 const Grid3D& grid,
							 DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3 )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;
	long w  = grid.w;

	double dx = grid.dx;
	double dy = grid.dy;
	double dz = grid.dz;

	double dx2 = dx*2.0;
	double dy2 = dy*2.0;
	double dz2 = dz*2.0;

	SSD = 0.0;
	long i, j, k;
	for( i = 1; i < m-1; i++ )
	{	for( j = 1; j < n-1; j++ )
		{	for( k = 1; k < p-1; k++ )
	{	SSD += 0.5*(T(i,j,k)-S(i,j,k))*(T(i,j,k)-S(i,j,k));	}	}	}

	double dTdx, dTdy, dTdz;

	double eps = 10e-5;
	const long size = m*n*p;
	
	long Imax = 128;

	DoubleArray2D p12(Imax,Imax);
	DoubleArray1D p1(Imax);
	DoubleArray1D p2(Imax);

	DoubleArray2D L(Imax,Imax);
	
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			for ( k = 0; k < p; k++ )
			{
				p12( (long)(0.5*S(i,j,k)), (long)(0.5*T(i,j,k)) ) += 1;
				p1(  (long)(0.5*S(i,j,k)) )						  += 1;
				p2(  (long)(0.5*T(i,j,k)) )						  += 1;
			}
		}
	}

	MI = 0.0;
	long i1, i2;
	for ( i1 = 0; i1 < Imax; i1++ )
	{
		for ( i2 = 0; i2 < Imax; i2++ )
		{
			if( p12(i1,i2) < eps || p1(i1)*p2(i2) < eps )
			{	L(i1,i2) = 1.0;	}
			else
			{	L(i1,i2) = 1.0 + log( size * p12(i1,i2)/(p1(i1)*p2(i2)) );	}
			
			MI += p12(i1,i2) * (L(i1,i2)-1.0);	// unscaled MI; 
		}										// Calculated at beginning of each iteration
	}
	cout << "MI = " << MI/size << "  ";

	DoubleArray2D K(Imax,Imax);
	K = dGdy( Imax, StD );
	K = fftshift2D(K);
	
	convolution2D( L, K );		// smooth L with dGaussian/dy

	double temp;
	for ( i = 1; i < m-1; i++ )
	{
		for ( j = 1; j < n-1; j++ )
		{	
			for( k = 1; k < p-1; k++ )
			{
				dTdx = (T(i+1,j,k) - T(i-1,j,k)) / dx2;
				dTdy = (T(i,j+1,k) - T(i,j-1,k)) / dy2;
				dTdz = (T(i,j,k+1) - T(i,j,k-1)) / dz2;

				temp = L( (long)(0.5*S(i,j,k)), (long)(0.5*T(i,j,k)) );
				
				f1(i,j,k) = temp * dTdx;
				f2(i,j,k) = temp * dTdy;
				f3(i,j,k) = temp * dTdz;

			}
		}
	}




	//////////////////////////////////////////////////

	DoubleArray3D d1_d1(m,n,p), d1_d2(m,n,p), d1_d3(m,n,p);
	DoubleArray3D d2_d1(m,n,p), d2_d2(m,n,p), d2_d3(m,n,p);
	DoubleArray3D d3_d1(m,n,p), d3_d2(m,n,p), d3_d3(m,n,p);

	Jmin = 10.0;
	Jmax = 0.0;

	for ( i = 1; i < m-1; i++ ) {
		for ( j = 1; j < n-1; j++ ) {
			for ( k = 1; k < p-1; k++ ) {
				d1_d1(i,j,k) = 1.0 - (u1(i+1,j,k)-u1(i-1,j,k))/2.0;
				d1_d2(i,j,k) =	   - (u1(i,j+1,k)-u1(i,j-1,k))/2.0;
				d1_d3(i,j,k) =	   - (u1(i,j,k+1)-u1(i,j,k-1))/2.0;

				d2_d1(i,j,k) =	   - (u2(i+1,j,k)-u2(i-1,j,k))/2.0;
				d2_d2(i,j,k) = 1.0 - (u2(i,j+1,k)-u2(i,j-1,k))/2.0;
				d2_d3(i,j,k) =     - (u2(i,j,k+1)-u2(i,j,k-1))/2.0;

				d3_d1(i,j,k) =	   - (u3(i+1,j,k)-u3(i-1,j,k))/2.0;
				d3_d2(i,j,k) =     - (u3(i,j+1,k)-u3(i,j-1,k))/2.0;
				d3_d3(i,j,k) = 1.0 - (u3(i,j,k+1)-u3(i,j,k-1))/2.0;

				J(i,j,k) =  d1_d1(i,j,k)*d2_d2(i,j,k)*d3_d3(i,j,k)
					      + d2_d1(i,j,k)*d3_d2(i,j,k)*d1_d3(i,j,k)
						  + d1_d2(i,j,k)*d2_d3(i,j,k)*d3_d1(i,j,k)
						  - d1_d3(i,j,k)*d2_d2(i,j,k)*d3_d1(i,j,k)
						  - d1_d2(i,j,k)*d2_d1(i,j,k)*d3_d3(i,j,k)
						  - d2_d3(i,j,k)*d3_d2(i,j,k)*d1_d1(i,j,k);

				if (J(i,j,k) < Jmin)  Jmin = J(i,j,k);
				if (J(i,j,k) > Jmax)  Jmax = J(i,j,k);
			}
		}
	}
	applyBC( J, 1, BCtype );	// for energy/StD calculation
}

//#####################################################################

void J_Reg3D::evaluate_v( DoubleArray3D& f1, DoubleArray3D& f2, DoubleArray3D& f3, 
						  const Grid3D& grid )
{
	long kk = 0;

	if( smoothing == 1 )			// GAUSSIAN SMOOTHING
	{
		f1 = (-1.0)*f1;
		f2 = (-1.0)*f2;
		f3 = (-1.0)*f3;

		convolution3D( f1, K );
		convolution3D( f2, K );
		convolution3D( f3, K );
	}

	else if( smoothing == 3 )	// STATIONARY NS (LINEAR ELASTICA)
	{
		long m  = grid.m;
		long n  = grid.n;
		long p  = grid.p;

		double dx = grid.dx;

		DoubleArray3D v1(m,n,p);
		DoubleArray3D v2(m,n,p);
		DoubleArray3D v3(m,n,p);

		double accuracy = 1.0e-5;

		MGsystem3D MG;
		StationaryNS3D PDEsystem;

		DoubleArray3D r1(m,n,p);			// residual matrix
		DoubleArray3D r2(m,n,p);			// residual matrix
		DoubleArray3D r3(m,n,p);			// residual matrix

		PDEsystem.findResidual( v1, v2, v3, f1, f2, f3, r1, r2, r3, dx );

		double L2residual1 = findL2( r1, 4 );
		double L2residual2 = findL2( r2, 4 );
		double L2residual3 = findL2( r3, 4 );

		cout << 0 << "\t" << L2residual1 << "\t" << L2residual2 << "\t" << L2residual3;

		double L2residual1_old = L2residual1;
		double L2residual2_old = L2residual2;
		double L2residual3_old = L2residual3;

		double convergenceFactorResidual1, convergenceFactorResidual2, convergenceFactorResidual3;

		MG.setPDE(&PDEsystem);	// passing in a reference (e.g. a pointer value)
								// to the derived class
		MG.solver = 6;

		bool condition = true;
		while(condition)
		{	kk = kk + 1;

			MG.MultiGridCycle( v1, v2, v3, f1, f2, f3, dx );	cout << "           *** MultiGrid Recursive ***     " << endl;

			PDEsystem.findResidual( v1, v2, v3, f1, f2, f3, r1, r2, r3, dx );
			L2residual1 = findL2( r1, 4 );
			L2residual2 = findL2( r2, 4 );
			L2residual3 = findL2( r3, 4 );
			convergenceFactorResidual1 = L2residual1 / L2residual1_old;
			convergenceFactorResidual2 = L2residual2 / L2residual2_old;
			convergenceFactorResidual3 = L2residual3 / L2residual3_old;

			cout << kk << "\t" << L2residual1 << "\t" << convergenceFactorResidual1 << "\t" << L2residual2 << "\t" << convergenceFactorResidual2 << "\t" << L2residual3 << "\t" << convergenceFactorResidual3 << endl;
			
			if( L2residual1 < accuracy && L2residual2 < accuracy && L2residual3 < accuracy )
				break;

			L2residual1_old = L2residual1;
			L2residual2_old = L2residual2;
			L2residual3_old = L2residual3;
		}

		f1 = v1;
		f2 = v2;
		f3 = v3;	// f1,f2 are returned by reference, holding v1,v2
	}

	else
	{
		cout << "Incorrect specification of FIELD FORCE SMOOTHING!" << endl;
		exit(1);
	}
}

//#####################################################################

DoubleArray3D J_Reg3D::getJacobian( const DoubleArray3D& u1, const DoubleArray3D& u2, const DoubleArray3D& u3,
								    const Grid3D& grid )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;

	DoubleArray3D d1_d1(m,n,p), d1_d2(m,n,p), d1_d3(m,n,p);
	DoubleArray3D d2_d1(m,n,p), d2_d2(m,n,p), d2_d3(m,n,p);
	DoubleArray3D d3_d1(m,n,p), d3_d2(m,n,p), d3_d3(m,n,p);
	
	DoubleArray3D Jac(m,n,p);
	
	long i,j,k;
	for ( i = 1; i < m-1; i++ ) {
		for ( j = 1; j < n-1; j++ ) {
			for ( k = 1; k < p-1; k++ ) {

				d1_d1(i,j,k) = 1.0 - (u1(i+1,j,k)-u1(i-1,j,k))/2.0;
				d1_d2(i,j,k) =	   - (u1(i,j+1,k)-u1(i,j-1,k))/2.0;
				d1_d3(i,j,k) =	   - (u1(i,j,k+1)-u1(i,j,k-1))/2.0;

				d2_d1(i,j,k) =	   - (u2(i+1,j,k)-u2(i-1,j,k))/2.0;
				d2_d2(i,j,k) = 1.0 - (u2(i,j+1,k)-u2(i,j-1,k))/2.0;
				d2_d3(i,j,k) =     - (u2(i,j,k+1)-u2(i,j,k-1))/2.0;

				d3_d1(i,j,k) =	   - (u3(i+1,j,k)-u3(i-1,j,k))/2.0;
				d3_d2(i,j,k) =     - (u3(i,j+1,k)-u3(i,j-1,k))/2.0;
				d3_d3(i,j,k) = 1.0 - (u3(i,j,k+1)-u3(i,j,k-1))/2.0;

				Jac(i,j,k) =  d1_d1(i,j,k)*d2_d2(i,j,k)*d3_d3(i,j,k)
					        + d2_d1(i,j,k)*d3_d2(i,j,k)*d1_d3(i,j,k)
						    + d1_d2(i,j,k)*d2_d3(i,j,k)*d3_d1(i,j,k)
						    - d1_d3(i,j,k)*d2_d2(i,j,k)*d3_d1(i,j,k)
						    - d1_d2(i,j,k)*d2_d1(i,j,k)*d3_d3(i,j,k)
						    - d2_d3(i,j,k)*d3_d2(i,j,k)*d1_d1(i,j,k);
			}
		}
	}

	return Jac;
}

//#####################################################################

void J_Reg3D::outputJinv( long kk, long factor, const Grid3D& grid )
{
	long m  = grid.m;
	long n  = grid.n;
	long p  = grid.p;

	FileIO_3D IO;

	DoubleArray3D Jinv(m,n,p);

	long i,j,k;
	for( k = 0; k < p; k++ )
	{	for( j = 0; j < n; j++ )
		{	for( i = 0; i < m; i++ )
			{	Jinv(i,j,k) = 1.0/J(i,j,k);	}	}	}
	IO.write_bin_usi( Jinv, "Jinv", kk, factor, grid.bxo, grid.byo, grid.bzo );

}

//#####################################################################

void J_Reg3D::outputParameters( const Grid3D& grid, const long TimeSteps, const long outputCount )
{
    ostringstream outs;
    char fileName[256];
    outs.str("");

	outs << "parameters.dat";

	strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );

	FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )
    {
		printf( "The file %s could not be  opened\n",fileName);
		return;
    }

	fprintf(dataFile, "%-10.5e \n", double(grid.m) );
	fprintf(dataFile, "%-10.5e \n", double(grid.n) );
	fprintf(dataFile, "%-10.5e \n", double(grid.p) );

	fprintf(dataFile, "%-10.5e \n", lambda );
	fprintf(dataFile, "%-10.5e \n", double(smoothing) );
	fprintf(dataFile, "%-10.5e \n", triangle_u );
	fprintf(dataFile, "%-10.5e \n", double(TimeSteps) );
	fprintf(dataFile, "%-10.5e \n", double(outputCount) );
	
	fprintf(dataFile, "%-10.5e \n", StD );
	fprintf(dataFile, "%-10.5e \n", double(heat_eqn_TimeSteps) );

	fclose(dataFile);
}
