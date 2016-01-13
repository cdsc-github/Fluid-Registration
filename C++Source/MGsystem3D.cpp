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

#include "MGsystem3D.h"
#include "Norms.h"

#define MAX( a, b ) (((a)>(b)) ? (a) : (b) )
#define PI 4.0*atan2(1.0,1.0)

//#####################################################################

MGsystem3D::MGsystem3D()
{
	nr = 2;
	ns = 50000;		// ns is used only in Two-Grid solvers.
	np = 1;

	numVcyclesFMG = 1;

	R_RestrictionOrder = 2;		// restriction order for residual transfter
	U_RestrictionOrder = 0;		// restriction order for u transfer
}

//#####################################################################

void MGsystem3D::setPDE( PDE_System3D *pde )
{
	pPDE = pde;
	// & from obj to ptr
	// * from ptr to obj
}

//#####################################################################

void MGsystem3D::setNumberGSRelaxations( void )
{
	if( solver == 1)							// Gauss-Seidel
	{
		// 
	}
	else if( solver == 2 || solver == 3 )		// Two-Grid
	{
		cout << "Number of  pre-smoothing (nr) relaxations on Fine grid  : ";
		cin >> nr;
		cout << "Number of      smoothing (ns) relaxations on Coarse grid: ";
		cin >> ns;
		cout << "Number of post-smoothing (np) relaxations on Fine grid  : ";
		cin >> np;
	}
	else if( solver == 4 || solver == 5 || solver == 6 || solver == 7 )	// MultiGrid
	{		
		cout << "Number of  pre-smoothing (nr) relaxations : ";
		cin >> nr;
		cout << "Number of post-smoothing (np) relaxations : ";
		cin >> np;
	}
	else if( solver == 8 || solver == 9 || solver == 10 || solver == 11 )	// Full MultiGrid (FMG)
	{
		cout << "Number of  pre-smoothing (nr) relaxations : ";
		cin >> nr;
		cout << "Number of post-smoothing (np) relaxations : ";
		cin >> np;
		cout << "Number of V-cycles per level (default = 1): ";
		cin >> numVcyclesFMG;
	}
	else
	{
		cout << solver << " is not a valid choice of solver!" << endl;
		exit(1);
	}
}

//#####################################################################

void MGsystem3D::setNumberGSRelaxations( long preSmooth, long Smooth, long postSmooth )
{
	nr = preSmooth;
	ns = Smooth;
	np = postSmooth;
}

//#####################################################################

void MGsystem3D::initialize( DoubleArray3D& A, const Grid3D& grid, const long& choice )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();
	long p = A.getIndex3Size();

	long i,j,k;
	double xx = 0.0, yy = 0.0, zz = 0;

	switch( choice )
	{
	case 0:		// A = 0
		for( k=0; k <= p-1; k++ )
		{	for( j=0; j <= n-1; j++ )
			{	for( i=0; i <= m-1; i++ )
				{	A(i,j,k) = 0.0;	}	}	}
		cout << "\t 0 initialization" << endl;
		break;

	case 1:		// A = random
		for( k=0; k <= p-1; k++ )
		{	for( j=0; j <= n-1; j++ )
			{	for( i=0; i <= m-1; i++ )
				{	A(i,j,k) = -rand()/double(RAND_MAX);	}	}	}  // "-" => to make error positive
		cout << " \t Random initialization" << endl;
		break;

	case 2:		// A = 0
		for( k=0; k <= p-1; k++ )
		{	for( j=0; j <= n-1; j++ )
			{	for( i=0; i <= m-1; i++ )
				{	A(i,j,k) = 5.0;	}	}	}
		cout << "\t 5 initialization" << endl;
		break;

	default: cout << "Function has not been initialized." << endl;
		exit(1);
		break;
	}
}

//#####################################################################

void MGsystem3D::restriction( const DoubleArray3D& R, DoubleArray3D& R2, const long& order )
{
	//
	// Restrict quantity R to coarser grid.
	//
	long m = R.getIndex1Size();
	long n = R.getIndex2Size();
	long p = R.getIndex3Size();

	long m2 = (m-1)/2 + 1;
	long n2 = (n-1)/2 + 1;
	long p2 = (p-1)/2 + 1;

	long i, j, k;
	long I, J, K;

	switch( order )
	{
	case 0:		// INJECTION:
		for( I = 0; I <= m2-1; I++ )
		{	i = 2*I;
			for( J = 0; J <= n2-1; J++ )
			{	j = 2*J;
				for( K = 0; K <= p2-1; K++ )
				{	k = 2*K;
					R2(I,J,K) = R(i,j,k);		// injection
				}
			}
		}
		break;

	case 1:		// AVERAGING:
		for( I = 1; I < m2-1; I++ )
		{	i = 2*I;
			for( J = 1; J < n2-1; J++ )
			{	j = 2*J;
				for( K = 1; K < p2-1; K++ )
				{	k = 2*K;
					R2(I,J,K) = (1/27.0) * ( R(i,j,k) + R(i-1, j , k ) + R(i+1, j , k ) 
													  + R( i ,j-1, k ) + R( i ,j+1, k )
													  + R( i , j ,k-1) + R( i , j ,k+1)
													  + R(i-1,j-1, k ) + R(i-1,j+1, k ) 
													  + R(i+1,j-1, k ) + R(i+1,j+1, k ) 
													  + R(i-1, j ,k-1) + R(i-1, j ,k+1)
													  + R(i+1, j ,k-1) + R(i+1, j, k+1)
													  + R( i ,j-1,k-1) + R( i ,j-1,k+1)
													  + R( i ,j+1,k-1) + R( i ,j+1,k+1)
													  + R(i-1,j-1,k-1) + R(i-1,j-1,k+1)
													  + R(i-1,j+1,k-1) + R(i-1,j+1,k+1)
													  + R(i+1,j-1,k-1) + R(i+1,j-1,k+1)
													  + R(i+1,j+1,k-1) + R(i+1,j+1,k+1) );
				}
			}
		}
		break;

	case 2:		// FULL WEIGHTING:
		for( I = 1; I < m2-1; I++ )
		{	i = 2*I;
			for( J = 1; J < n2-1; J++ )
			{	j = 2*J;
				for( K = 1; K < p2-1; K++ )
				{	k = 2*K;
					R2(I,J,K) = (1/64.0) * ( 8.0*R(i,j,k) + 4.0*R(i-1, j , k ) + 4.0*R(i+1, j , k ) 
														  + 4.0*R( i ,j-1, k ) + 4.0*R( i ,j+1, k )
														  + 4.0*R( i , j ,k-1) + 4.0*R( i , j ,k+1)
														  + 2.0*R(i-1,j-1, k ) + 2.0*R(i-1,j+1, k ) 
														  + 2.0*R(i+1,j-1, k ) + 2.0*R(i+1,j+1, k ) 
														  + 2.0*R(i-1, j ,k-1) + 2.0*R(i-1, j ,k+1)
														  + 2.0*R(i+1, j ,k-1) + 2.0*R(i+1, j, k+1)
														  + 2.0*R( i ,j-1,k-1) + 2.0*R( i ,j-1,k+1)
														  + 2.0*R( i ,j+1,k-1) + 2.0*R( i ,j+1,k+1)
														  +     R(i-1,j-1,k-1) +     R(i-1,j-1,k+1)
														  +     R(i-1,j+1,k-1) +     R(i-1,j+1,k+1)
														  +     R(i+1,j-1,k-1) +     R(i+1,j-1,k+1)
														  +     R(i+1,j+1,k-1) +     R(i+1,j+1,k+1) );
				}
			}
		}
		break;
	
	default: cout << "Incorrect specification of the interpolation order of restriction." << endl;
		exit(1);
		break;
	} // end switch
}

//#####################################################################

void MGsystem3D::prolongation( const DoubleArray3D& V2, DoubleArray3D& V, const long& order )
{
	// This routine prolongates quantity V2 from coarse to fine grid
	// using either injection, linear, quadratic, cubic, or 4th order interpolation.
	//
	// The dimensions of V2 have to be odd.
	//

	long m2 = V2.getIndex1Size();
	long n2 = V2.getIndex2Size();
	long p2 = V2.getIndex3Size();

	long i; long j; long k;
	long I; long J; long K;

	switch( order )
	{
		case 1:		// LINEAR INTERPOLATION
			for( K = 0; K <= p2-1; K++)
			{	k = 2*K;
				for( J = 0; J <= n2-1; J++)
				{	j = 2*J;
					for( I = 0; I < m2-1; I++ )
					{	i = 2*I;
						V( i ,j,k)   =  V2(I,J,K);
						V(i+1,j,k)   = (V2(I,J,K)+V2(I+1,J,K)) / 2.0;
					}
				}
			}

			for( K = 0; K <= p2-1; K++)
			{	k = 2*K;
				for( I = 0; I <= m2-2; I++ )
				{	i = 2*I;
					for( J = 0; J <= n2-2; J++)
					{	j = 2*J;
						V( i ,j+1,k) = (V( i ,j,k)+V( i ,j+2,k)) / 2.0;
						V(i+1,j+1,k) = (V(i+1,j,k)+V(i+1,j+2,k)) / 2.0;
					}
				}
			}

			for( I = 0; I <= m2-2; I++ )
			{	i = 2*I;
				for( J = 0; J <= n2-2; J++)
				{	j = 2*J;
					for( K = 0; K < p2-1; K++)
					{	k = 2*K;
						V( i , j ,k+1) = (V( i , j ,k)+V( i , j ,k+2)) / 2.0;
						V(i+1, j ,k+1) = (V(i+1, j ,k)+V(i+1, j ,k+2)) / 2.0;
						V( i ,j+1,k+1) = (V( i ,j+1,k)+V( i ,j+1,k+2)) / 2.0;
						V(i+1,j+1,k+1) = (V(i+1,j+1,k)+V(i+1,j+1,k+2)) / 2.0;
					}
				}
			}

			// Right column:
			for( K = 0; K <= p2-1; K++)
			{	k = 2*K;
				for( J = 0; J < n2-1; J++ )
				{	j = 2*J;
					I = m2-1; i = 2*I;
					V(i, j ,k)   =    V2(I,J,K);
					V(i,j+1,k)   =   (V2(I,J,K)+V2(I,J+1,K)) / 2.0;
				}
				I = m2-1; J = n2-1; i = 2*I; j = 2*J;
				V(i,j,k) = V2(I,J,K);		// Upper-right corner
			}
			for( K = 0; K < p2-1; K++)
			{	k = 2*K;
				for( J = 0; J < n2-1; J++ )
				{	j = 2*J;
					I = m2-1; i = 2*I;
					V(i, j ,k+1) = (V(i, j ,k)+V(i, j ,k+2)) / 2.0;
					V(i,j+1,k+1) = (V(i,j+1,k)+V(i,j+1,k+2)) / 2.0;
				}
			}

			// Top row:
			for( K = 0; K < p2-1; K++)
			{	k = 2*K;
				for( I = 0; I < m2-1; I++ )
				{	i = 2*I;
					J = n2-1; j = 2*J;
					V( i ,j,k+1) = (V( i ,j,k)+V( i ,j,k+2)) / 2.0;
					V(i+1,j,k+1) = (V(i+1,j,k)+V(i+1,j,k+2)) / 2.0;
				}
				I = m2-1; J = n2-1; i = 2*I; j = 2*J;
				V(i,j,k+1) = (V(i,j,k)+V(i,j,k+2))/2.0;		// Upper-right corner
			}

		break;

		default: cout << "Incorrect specification of the interpolation order of prolongation." << endl;
			exit(1);
		break;
	}  // end switch
}

//#####################################################################

void MGsystem3D::correct( const DoubleArray3D& V, DoubleArray3D& u )
{
	long m = u.getIndex1Size();
	long n = u.getIndex2Size();
	long p = u.getIndex3Size();

	long i,j,k;

	for( i=1; i < m-1; i++ )
	{	for( j=1; j < n-1; j++ )
		{	for( k=1; k < p-1; k++ )
			{	u(i,j,k) = u(i,j,k) + V(i,j,k);	}	}	}

}

//#####################################################################

void MGsystem3D::MultiGridCycle( DoubleArray3D& u1, 
								 DoubleArray3D& u2,
								 DoubleArray3D& u3,
								 const DoubleArray3D& f1,
								 const DoubleArray3D& f2,
								 const DoubleArray3D& f3,
								 const double& dx )
{
    long m = u1.getIndex1Size();
	long n = u1.getIndex2Size();
	long p = u1.getIndex3Size();

	/*
	// solve EXACTLY if there is only one interior point:
	if (m == 3) {
		pPDE->solveExactly( u1, u2, u3, f1, f2, f3, dx );
		return;
    }
	*/
	long kk = 0;
	// solve ALMOST EXACTLY if there are 10 or less interior points:
	if (m <= 10) {
		for ( kk = 1; kk <= 12; kk++ )
		{	pPDE->applyRelaxationGS( u1, u2, u3, f1, f2, f3, dx );	}
		return;
    }

    // do a few pre-smoothing Gauss-Seidel steps
    for ( kk = 1; kk <= nr; kk++ )
		pPDE->applyRelaxationGS( u1, u2, u3, f1, f2, f3, dx );

	//
    // find the residual:
	//
    DoubleArray3D R1(m,n,p);
	DoubleArray3D R2(m,n,p);
	DoubleArray3D R3(m,n,p);
	pPDE->findResidual( u1, u2, u3, f1, f2, f3, R1, R2, R3, dx );

	long m2 = (m-1)/2 + 1;
	long n2 = (n-1)/2 + 1;
	long p2 = (p-1)/2 + 1;

	//
    // transfer residual to coarser grid:
    //
	DoubleArray3D R1_coarse(m2,n2,p2);
	DoubleArray3D R2_coarse(m2,n2,p2);
	DoubleArray3D R3_coarse(m2,n2,p2);
	restriction( R1, R1_coarse, R_RestrictionOrder );
	restriction( R2, R2_coarse, R_RestrictionOrder );
	restriction( R3, R3_coarse, R_RestrictionOrder );

    // initialize correction V on coarse grid to zero
    DoubleArray3D V1_coarse(m2,n2,p2);
	DoubleArray3D V2_coarse(m2,n2,p2);
	DoubleArray3D V3_coarse(m2,n2,p2);

    // call twoGrid recursively
    MultiGridCycle( V1_coarse, V2_coarse, V3_coarse,
					R1_coarse, R2_coarse, R3_coarse, 2*dx );

	//
    // prolongate V1,V2,V3 to fine grid:
	//
    DoubleArray3D V1(m,n,p);
	DoubleArray3D V2(m,n,p);
	DoubleArray3D V3(m,n,p);
	prolongation( V1_coarse, V1, 1 );
	prolongation( V2_coarse, V2, 1 );
	prolongation( V3_coarse, V3, 1 );

	//
	// correct u1, u2, u3
	//
	correct( V1, u1 );
	correct( V2, u2 );
	correct( V3, u3 );

	//
	// do a few post-smoothing Gauss-Seidel steps:
	//
	for( kk = 1; kk <= np; kk++ )
	{	pPDE->applyRelaxationGS( u1, u2, u3, f1, f2, f3, dx );	}
}

//#####################################################################

void MGsystem3D::MultiGridCycle_FAS( DoubleArray3D& u1, 
								 DoubleArray3D& u2,
								 DoubleArray3D& u3,
								 const DoubleArray3D& f1,
								 const DoubleArray3D& f2,
								 const DoubleArray3D& f3,
								 const double& dx )
{
    long m = u1.getIndex1Size();
	long n = u1.getIndex2Size();
	long p = u1.getIndex3Size();
	
	long m2 = (m-1)/2 + 1;
	long n2 = (n-1)/2 + 1;
	long p2 = (p-1)/2 + 1;

	DoubleArray3D u1_coarse(m2,n2,p2);
	DoubleArray3D u2_coarse(m2,n2,p2);
	DoubleArray3D u3_coarse(m2,n2,p2);
	DoubleArray3D V1(m,n,p);
	DoubleArray3D V2(m,n,p);
	DoubleArray3D V3(m,n,p);
	DoubleArray3D R1(m,n,p);
	DoubleArray3D R2(m,n,p);
	DoubleArray3D R3(m,n,p);
	DoubleArray3D R1_coarse(m2,n2,p2);
	DoubleArray3D R2_coarse(m2,n2,p2);
	DoubleArray3D R3_coarse(m2,n2,p2);
	DoubleArray3D L2u1_coarse(m2,n2,p2);
	DoubleArray3D L2u2_coarse(m2,n2,p2);
	DoubleArray3D L2u3_coarse(m2,n2,p2);
	DoubleArray3D f2bar1(m2,n2,p2);
	DoubleArray3D f2bar2(m2,n2,p2);
	DoubleArray3D f2bar3(m2,n2,p2);

	// solve exactly if there is only one interior point:
	if (m == 3) {
		pPDE->solveExactly( u1, u2, u3, f1, f2, f3, dx );
		return;
    }

	long kk = 0;

    // do a few pre-smoothing Gauss-Seidel steps
    for ( kk = 1; kk <= nr; kk++ )
        pPDE->applyRelaxationGS( u1, u2, u3, f1, f2, f3, dx );

	// restrict to coarse grid:  u -> u_coarse
	restriction( u1, u1_coarse, U_RestrictionOrder );
	restriction( u2, u2_coarse, U_RestrictionOrder );
	restriction( u3, u3_coarse, U_RestrictionOrder );

	// find Laplace operator for u
	pPDE->operatorL( u1_coarse, u2_coarse, u3_coarse, L2u1_coarse, L2u2_coarse, L2u3_coarse, 2.0*dx );
	
    // find the residual:
	pPDE->findResidual( u1, u2, u3, f1, f2, f3, R1, R2, R3, dx );

	// restrict to coarse grid
	restriction( R1, R1_coarse, R_RestrictionOrder );
	restriction( R2, R2_coarse, R_RestrictionOrder );
	restriction( R3, R3_coarse, R_RestrictionOrder );

	f2bar1 = R1_coarse + L2u1_coarse;
	f2bar2 = R2_coarse + L2u2_coarse;
	f2bar3 = R3_coarse + L2u3_coarse;

    DoubleArray3D V1_coarse(m2,n2,p2);
	DoubleArray3D V2_coarse(m2,n2,p2);
	DoubleArray3D V3_coarse(m2,n2,p2);
	DoubleArray3D U1_coarse = u1_coarse;
	DoubleArray3D U2_coarse = u2_coarse;
	DoubleArray3D U3_coarse = u3_coarse;

    // call MultiGridCycle_FAS recursively
    MultiGridCycle_FAS( U1_coarse, U2_coarse, U3_coarse, f2bar1, f2bar2, f2bar3, 2.0*dx );

	V1_coarse = U1_coarse - u1_coarse;
	V2_coarse = U2_coarse - u2_coarse;
	V3_coarse = U3_coarse - u3_coarse;

    // prolongate V1,V2,V3 to fine grid:
	prolongation( V1_coarse, V1, 1 );
	prolongation( V2_coarse, V2, 1 );
	prolongation( V3_coarse, V3, 1 );

	// correct u1, u2, u3
	correct( V1, u1 );
	correct( V2, u2 );
	correct( V3, u3 );

	//
	// do a few post-smoothing Gauss-Seidel steps:
	//
	for( kk = 1; kk <= np; kk++ )
	{	pPDE->applyRelaxationGS( u1, u2, u3, f1, f2, f3, dx );	}
}

//#####################################################################

void MGsystem3D::outputParameters( const double& accuracy, const Grid3D& grid, 
								   const long& TimeSteps, const long& outputCount, const double& timeTaken )
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
		printf( "The file %s could not be opened\n", fileName);
		return;
    }

	fprintf(dataFile, "%-10.5e \n", double(grid.m) );
	fprintf(dataFile, "%-10.5e \n", double(grid.n) );
	fprintf(dataFile, "%-10.5e \n", double(grid.p) );
	fprintf(dataFile, "%-10.5e \n", double(pPDE->problemNumber) );
	fprintf(dataFile, "%-10.5e \n", double(solver) );
	fprintf(dataFile, "%-10.5e \n", double(nr) );
	fprintf(dataFile, "%-10.5e \n", double(ns) );
	fprintf(dataFile, "%-10.5e \n", double(np) );
	fprintf(dataFile, "%-10.5e \n", double(numVcyclesFMG) );
	fprintf(dataFile, "%-10.5e \n", double(R_RestrictionOrder) );
	fprintf(dataFile, "%-10.5e \n", double(U_RestrictionOrder) );
	fprintf(dataFile, "%-10.5e \n", accuracy );
	fprintf(dataFile, "%-10.5e \n", double(TimeSteps) );
	fprintf(dataFile, "%-10.5e \n", double(outputCount) );
	fprintf(dataFile, "%-10.5e \n", timeTaken);

	fclose(dataFile);
}
