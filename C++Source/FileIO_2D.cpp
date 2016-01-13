#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <fstream>
#include <math.h>

//
//#####################################################################
//						FileIO_2D.cpp
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// version:  11/25/2008
//
//#####################################################################
//

#include "FileIO_2D.h"
#include "dbh.h"
#include "common_routines.h"

//#####################################################################

inline void endian_swap(unsigned short& x)
{
    x = (x>>8) | 
        (x<<8);
}

inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}

inline void endian_swap(signed short& x)
{		// Igor wrote this routine
    x = (x>>8) | 
        (x<<8);
	x = (unsigned short) x;
}

//
//#####################################################################
//						  GET IMAGE INFO
//#####################################################################
//
void FileIO_2D::getImageInfo( string filename, long& m, long& n )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	double intensity;
	char type[100];

	infile >> type;			// Read the Header part of PGM file that contains "P2"	
	infile >> m >> n;		// Read the dimension of the PGM image
	infile >> intensity;	// Read the maximum gray level. i.e. 255	
}
//
//#####################################################################
//							PGM ASCII READ
//#####################################################################
//
DoubleArray2D FileIO_2D::readPGM( string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long m; long n;

	char tmp[100];
	infile >> tmp;		// Read the Header part of PGM file that contains "P2"	
	infile >> m >> n;		// Read the dimension of the PGM image
	infile >> tmp;		// Read the maximum gray level. i.e. 255	

	DoubleArray2D A(m,n);

	long i;  long j;
	for(j = 0; j < n; j++)
    {
		for(i = 0; i < m; i++)
		{
			infile >> A(i,j);
		}
	}
	infile.close();

	return A;
}

//#####################################################################

DoubleArray2D FileIO_2D::readPGM( string filename, 
								  long& m, long& n, 
								  const long bx, const long by )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	char tmp[100];
	infile >> tmp;		// Read the Header part of PGM file that contains "P2"	
	infile >> m >> n;	// Read the dimension of the PGM image
	infile >> tmp;		// Read the maximum gray level. i.e. 255	

	m = m + 2*bx;
	n = n + 2*bx;
	DoubleArray2D A(m,n);

	long i, j;
	for(j = by; j < n-by; j++)
    {
		for(i = bx; i < m-bx; i++)
		{
			infile >> A(i,j);
		}
	}
	infile.close();

	return A;
}
//
//#####################################################################
//							ASCII READ
//#####################################################################
//
void FileIO_2D::readDAT2D( DoubleArray2D& A, long m, long n, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long i;  long j;
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < m; i++)
		{
			infile >> A(i,j);
		}
	}
	infile.close();
}

//#####################################################################

void FileIO_2D::readDAT2D( DoubleArray2D& A, string filename )
{
	long m = A.getIndex1Size();
	long n = A.getIndex2Size();

	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long i;  long j;
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < m; i++)
		{
			infile >> A(i,j);
		}
	}
	infile.close();
}

//#####################################################################

void FileIO_2D::readDAT1D( DoubleArray1D& A, string filename )
{
	long m = A.getIndex1Size();

	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long i;
	for(i = 0; i < m; i++)
	{
		infile >> A(i);
	}
	infile.close();
}

//
//#####################################################################
//							PGM ASCII WRITE
//#####################################################################
//
void FileIO_2D::writePGM( const DoubleArray2D& A, const string& str )
{
	ostringstream outs;
	outs << str << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << endl;
	outfile << 255 << endl;

    long i; long j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			//if( A(i,j) < 0 )
			//{	outfile << setw(5) << 0 << " ";	}
			//else if( A(i,j) > 255 )
			//{	outfile << setw(5) << 255 << " ";	}
			//else
			{	outfile <<  setw(5) << round_to_int( A(i,j) ) << " ";	}
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################

void FileIO_2D::writePGM( const DoubleArray2D& A, const string& str, 
						  const long& stepCount )
{
	ostringstream outs;
	outs << str << stepCount << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	findMaxMin( A );

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << endl;
	outfile << 255 << endl;

    long i; long j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			if( A(i,j) < 0 )
			{	outfile << setw(5) << 0 << " ";	}
			else if( A(i,j) > 255 )
			{	outfile << setw(5) << 255 << " ";	}
			else
			{	outfile << setw(5) << round_to_int( A(i,j) ) << " ";	}
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################

void FileIO_2D::writePGM( const DoubleArray2D& A, const string& str,
						  const long& stepCount,
						  const long bx, const long by )
{
	ostringstream outs;
	outs << str << stepCount << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() - 2*bx << " " << A.getIndex2Size() - 2*by << endl;
	outfile << 255 << endl;

    long i; long j;
	for(j = A.getIndex2Begin() + by; j <= A.getIndex2End() - by; j++)
    {
		for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
		{			
			if( A(i,j) < 0 )
			{	outfile << setw(5) << 0 << " ";	}
			else if( A(i,j) > 255 )
			{	outfile << setw(5) << 255 << " ";	}
			else
			{
				outfile <<  setw(5) << round_to_int( A(i,j) ) << " ";
			}
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################
void FileIO_2D::writePGM( const DoubleArray2D& A, const string& str, 
						  const long& stepCount, 
						  const long& b )
{
	ostringstream outs;
	outs << str << stepCount << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size()-2*b << " " << A.getIndex2Size()-2*b << endl;
	outfile << 255 << endl;

    long i; long j;
	for(j = A.getIndex2Begin()+b; j <= A.getIndex2End()-b; j++)
    {
		for(i = A.getIndex1Begin()+b; i <= A.getIndex1End()-b; i++)
		{
			outfile << setw(5) << round_to_int( A(i,j) ) << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//
//#####################################################################
//					WRITE SCALED PGM
//#####################################################################
//
void FileIO_2D::writePGMsc( const DoubleArray2D& AA, const string& str )
{

	ostringstream outs;
	outs << str << ".pgm";
	string filename = outs.str();

	DoubleArray2D A = AA;

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	double minA = A(1,1);

	long i; long j;
	
	// Find the minimum of A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			minA = (minA < A(i,j)) ? minA : A(i,j);
		}
	}
	
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) - minA;		// will have min(A) = 0
		}
	}

	double maxA = A(1,1);

	// Find the maximum of new A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			maxA = (maxA > A(i,j)) ? maxA : A(i,j);
		}
	}

	// Rescale:
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) * 255.0 / maxA;
		}
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << endl;
	outfile << 255 << endl;

    
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			outfile <<  setw(5) << (long)A(i,j) << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################

void FileIO_2D::writePGMsc( const DoubleArray2D& AA, const string& str, 
						    const long& stepCount )
{
	ostringstream outs;
	outs << str << stepCount << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}


	DoubleArray2D A = AA;

	double minA =  100000000;
	double maxA = -100000000;

	long i; long j;
	
	// Find the minimum of A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			minA = (minA < A(i,j)) ? minA : A(i,j);
		}
	}

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) - minA;		// will have min(A) = 0
		}
	}

	// Find the maximum of new A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			maxA = (maxA > A(i,j)) ? maxA : A(i,j);
		}
	}

	// Rescale:
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) * 255.0 / maxA;
		}
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << endl;
	outfile << 255 << endl;

    
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			outfile <<  setw(5) << (long)A(i,j) << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################

void FileIO_2D::writePGMsc( const DoubleArray2D& AA, const string& str, 
						    const long& stepCount,
							const long bx, const long by)
{
	ostringstream outs;
	outs << str << stepCount << ".pgm";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}


	DoubleArray2D A = AA;

	double minA =  100000000;
	double maxA = -100000000;

	long i; long j;
	
	// Find the minimum of A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			minA = (minA < A(i,j)) ? minA : A(i,j);
		}
	}

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) - minA;		// will have min(A) = 0
		}
	}

	// Find the maximum of new A
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			maxA = (maxA > A(i,j)) ? maxA : A(i,j);
		}
	}

	// Rescale:
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) * 255.0 / maxA;
		}
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() - 2*bx << " " << A.getIndex2Size() - 2*by << endl;
	outfile << 255 << endl;

    
	for(j = A.getIndex2Begin() + by ; j <= A.getIndex2End() - by; j++)
    {
		for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
		{
			outfile <<  setw(5) << (long)A(i,j) << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//
//#####################################################################
//					WRITE SCALED PGM
//
//          Consider pixels EXCEPT for: 
// r columns on the right,  l columns on the left,
// t rows on the top,       b rows on the bottom
//#####################################################################
//
void FileIO_2D::writePGMsc( const DoubleArray2D& AA, long r, long l, long t, long b, const string& str )
{
	ostringstream outs;
	outs << str << ".pgm";
	string filename = outs.str();

	DoubleArray2D A = AA;

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	double minA =  100000000;
	double maxA = -100000000;

	long i; long j;
	
	// Find the minimum of A inside the narrow boundary
	for(j = A.getIndex2Begin()+r; j <= A.getIndex2End()-l; j++)
    {
		for(i = A.getIndex1Begin()+t; i <= A.getIndex1End()-b; i++)
		{
			minA = (minA < A(i,j)) ? minA : A(i,j);
		}
	}

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) - minA;		// will have min(A) = 0 inside
		}
	}

	// Find the maximum of new A inside the narrow boundary
	for(j = A.getIndex2Begin()+r; j <= A.getIndex2End()-l; j++)
    {
		for(i = A.getIndex1Begin()+t; i <= A.getIndex1End()-b; i++)
		{
			maxA = (maxA > A(i,j)) ? maxA : A(i,j);
		}
	}

	// Rescale:
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			A(i,j) = A(i,j) * 255.0 / maxA;
		}
	}

	outfile << "P2" << endl;
	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << endl;
	outfile << 255 << endl;

    
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			outfile <<  setw(5) << (long)A(i,j) << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//
//#####################################################################
//							ASCII WRITE
//#####################################################################
//
void FileIO_2D::write_ascii( const vector2Dint& A, const string& str )
{
	ostringstream outs;
	outs << str << ".dat";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

    long i; long j;
	for(j = 0; j < long(A[0].size()); j++)
    {
		for(i = 0; i < long(A.size()); i++)
		{
			outfile <<  setw(5) << A[i][j] << " ";
		}
		outfile << endl;
    }
	outfile.close();
}

//#####################################################################

void FileIO_2D::write_ascii( const vector1Ddouble& A, const string& str )
{
	ostringstream outs;
	outs << str << ".dat";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	long i; 

		for(i = 0; i < long(A.size()); i++)
		{
			outfile <<  setw(5) << A[i] << " "; outfile << endl;
		}	
	outfile.close();
}

//#####################################################################

void FileIO_2D::write_ascii( const DoubleArray2D& A,
							 const string& str, 
							 const long& stepCount )
{
	//
	//  Create ostringstream for construction output file names
	// 
    ostringstream outs;
    char fileName[256];
    // 
    // Create file name of form strXXX.dat for step XXX 
    //
    outs.str("");
	outs << str << stepCount << ".dat";
    strcpy(fileName,(outs.str()).c_str());
	//strcpy_s( fileName, 256, (outs.str()).c_str() );	// 

    //
	//  Open and then write to a file
	//
    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		
    {
      printf( "The file %s could not be opened\n",fileName);
      return;
    }

	long i; long j;
	//
	//  Output the data.
	// 

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			fprintf( dataFile, "%-10.5e ", A(i,j));
		}
		fprintf( dataFile, "\n" );
	}

    fclose(dataFile);
}

//#####################################################################

void FileIO_2D::write_ascii( const DoubleArray2D& A, const string& str, 
							 const long& stepCount,
							 const long bx, const long by )
{
	//
	//  Create ostringstream for construction output file names
	// 
    ostringstream outs;
    char fileName[256];
    // 
    // Create file name of form strXXX.dat for step XXX 
    //
    outs.str("");
	outs << str << stepCount << ".dat";
    strcpy(fileName,(outs.str()).c_str());
	//strcpy_s( fileName, 256, (outs.str()).c_str() );	// 

    //
	//  Open and then write to a file
	//
    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		
    {
      printf( "The file %s could not be opened\n",fileName);
      return;
    }

	long i; long j;
	//
	//  Output the data.
	// 

	for(j = A.getIndex2Begin() + by; j <= A.getIndex2End() - by; j++)
    {
		for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
		{
			fprintf( dataFile, "%-10.5e ", A(i,j));
		}
		fprintf( dataFile, "\n" );
	}

    fclose(dataFile);
}

//#####################################################################

void FileIO_2D::write_ascii( const DoubleArray2D& A,
							 const string& str, 
							 const long& stepCount,
							 const long& w )
{
	//
	//  Create ostringstream for construction output file names
	// 
    ostringstream outs;
    char fileName[256];
    // 
    // Create file name of form strXXX.dat for step XXX 
    //
    outs.str("");
	outs << str << stepCount << ".dat";
	strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// 

    //
	//  Open and then write to a file
	//
    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
      printf( "The file %s could not be opened\n",fileName);
      return;
    }

	long i; long j;
	//
	//  Output the data.
	// 

	for(j = w; j < A.getIndex2Size()-w; j++)
    {
		for(i = w; i < A.getIndex1Size()-w; i++)
		{
			fprintf( dataFile, "%-10.5e ", A(i,j));
		}
		fprintf( dataFile, "\n" );
	}

    fclose(dataFile);
}

//#####################################################################

void FileIO_2D::write_ascii( const DoubleArray2D& A,
							 const string& str )
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
	outs << str << ".dat";
	strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

    FILE* dataFile;

    if( (dataFile = fopen(fileName, "w+" )) == NULL )
	//if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
      printf( "The file %s could not be opened\n",fileName);
      return;
    }

	long i; long j;

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
    {
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			fprintf( dataFile, "%-10.5e ", A(i,j));
		}
		fprintf( dataFile, "\n" );
	}

    fclose(dataFile);
}
//
//#####################################################################
//						WRITE 1D array
//#####################################################################
void FileIO_2D::write_ascii( const DoubleArray1D& A,
							 const string& str )
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
	outs << str << ".dat";
	strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
		printf( "The file %s could not be  opened\n",fileName);
		return;
    }

	long i;

	for(i = 0;  i < A.getSize();  i++)
	{	fprintf( dataFile, "%-10.5e \n", A(i) );	}
	
	fclose(dataFile);
}

//#####################################################################

void FileIO_2D::write_ascii( const DoubleArray1D& A,
							 const string& str,
							 const long& stepCount )
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
	outs << str << stepCount << ".dat";
	strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
		printf( "The file %s could not be  opened\n",fileName);
		return;
    }

	long i;

	for(i = 0;  i < A.getSize();  i++)
	{	fprintf( dataFile, "%-10.5e \n", A(i) );	}
	
	fclose(dataFile);
}

//#####################################################################
//					BINARY READ
//#####################################################################

void FileIO_2D::read_bin_uc( string filename, DoubleArray2D& A )
{
	//
	// This routine reads binary files that do NOT have 
	// a header info included.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}
	
	long i;  long j;
	unsigned char uc;

	cout << "READ: unsigned char: " << sizeof(uc) << " byte ( " << sizeof(uc) * 8 << " bit ) data reading." << endl;

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			// The following code retrieves an unsigned char from a file
			// and copies it to uc:
		    infile.read(reinterpret_cast <char*> (&uc),sizeof(unsigned char));
			
			// Convert unsigned char uc to double and assign to A(i,j):
			A(i,j) = uc;
		}
	}
}

//#####################################################################

void FileIO_2D::read_bin_usi( string filename, DoubleArray2D& A, bool byteswap )
{
	//
	// This routine reads binary files that do NOT have 
	// a header info included.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}
	
	long i;  long j;
	unsigned short int usi;
	
	cout << "READ: unsigned short int: " << sizeof(usi) << " byte ( " << sizeof(usi) * 8 << " bit ) data." << endl;

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			// The following code retrieves an unsigned short int from a file
			// and copies it to usi:
		    infile.read(reinterpret_cast <char*> (&usi),sizeof(unsigned short int));
			
			if( byteswap == true )
			{	// Byte-swapping:
				endian_swap(usi);	//cout << "Byte-Swapping" << endl;
			}
			// Convert unsigned short int usi to double and assign to A(i,j,k):
			A(i,j) = usi;
		}
	}
}

//#####################################################################

void FileIO_2D::read_bin_ssi( string filename, DoubleArray2D& A, bool byteswap )
{
	//
	// This routine reads binary files that do NOT have 
	// a header info included.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}
	
	long i;  long j;
	signed short int ssi;
	
	cout << "READ: signed short int: " << sizeof(ssi) << " byte ( " << sizeof(ssi) * 8 << " bit ) data." << endl;

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			// The following code retrieves an unsigned short int from a file
			// and copies it to usi:
		    infile.read(reinterpret_cast <char*> (&ssi),sizeof(signed short int));
			
			if( byteswap == true )
			{	// Byte-swapping:
				endian_swap(ssi);	//cout << "Byte-Swapping" << endl;
			}				

			// Convert signed short int ssi to double and assign to A(i,j,k):
			A(i,j) = ssi;

			/*if( A(i,j,k) != 0)	
			{	display( A, i, j, k );	
				exit(1);	
			}*/
		}
	}
}

//#####################################################################

void FileIO_2D::read_bin_float( string filename, DoubleArray2D& A )
{
	//
	// This routine reads binary files that do NOT have 
	// a header info included.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}
	
	long i;  long j;
	float flt;
	
	cout << "READ: float: " << sizeof(float) << " byte ( " << sizeof(float) * 8 << " bit ) data." << endl;

	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			// The following code retrieves an unsigned short int from a file
			// and copies it to usi:
		    infile.read(reinterpret_cast <char*> (&flt),sizeof(float));		

			// Convert float to double and assign to A(i,j,k):
			A(i,j) = flt;

			/*if( A(i,j,k) != 0)	
			{	display( A, i, j, k );	
				exit(1);	
			}*/
		}
	}
}

//
//#####################################################################
//							BINARY WRITE
//#####################################################################
//
void FileIO_2D::write_bin_usi( const DoubleArray2D& A, 
							   const string& str, 
							   long factor )
{	//
	// This routine WRITES .IMG and .HDR.
	//
	ostringstream outs;
	outs << str << ".img";
	string filename = outs.str();

	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	unsigned short num;

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i,j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			num = (unsigned short)( A(i,j)*factor );
			outfile.write((char*)&num,sizeof(unsigned short));
		}
	}
	make_header( A, str );
}

//#####################################################################

void FileIO_2D::write_bin_usi( const DoubleArray2D& A, 
							   const string& str, 
							   const long& stepCount, 
							   long factor )
{	//
	// This routine WRITES .IMG and .HDR.
	//
	ostringstream outs, outs2;
	outs << str << stepCount << ".img";
	string filename = outs.str();
	
	outs2 << str << stepCount;
	string basename = outs2.str();
	
	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	unsigned short num;

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i,j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			num = (unsigned short)( A(i,j)*factor );
			outfile.write((char*)&num,sizeof(unsigned short));
		}
	}

	make_header( A, basename );
}

//#####################################################################

void FileIO_2D::write_bin_usi( const DoubleArray2D& A, 
							   const string& str, 
							   const long& stepCount, 
							   long factor,
							   const long bx, const long by )
{	//
	// This routine WRITES .IMG and .HDR.
	//
	ostringstream outs, outs2;
	outs << str << stepCount << ".img";
	string filename = outs.str();
	
	outs2 << str << stepCount;
	string basename = outs2.str();
	
	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	unsigned short num;

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i,j;
	for(j = A.getIndex2Begin() + by; j <= A.getIndex2End() - by; j++)
	{
		for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
		{
			num = (unsigned short)( A(i,j)*factor );
			outfile.write((char*)&num,sizeof(unsigned short));
		}
	}

	make_header( A, basename, bx, by );
}

//#####################################################################

void FileIO_2D::make_header( const DoubleArray2D& A,
							 const string& str ) /* file x y z t datatype max min */
{
	ostringstream outs;
    char fileName[256];

    outs.str("");
	outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

	int m = A.getIndex1Size();
	int n = A.getIndex2Size();


	int i;
	struct dsr hdr;
	FILE *fp;
	static char DataTypes[9][12] = {"UNKNOWN", "BINARY", "CHAR", "SHORT", "INT","FLOAT", "COMPLEX", "DOUBLE", "RGB"};
	static int DataTypeSizes[9] = {0,1,8,16,32,32,64,64,24};

	memset(&hdr,0, sizeof(struct dsr));
	
	for(i=0;i<8;i++)
		hdr.dime.pixdim[i] = 0.0;

	hdr.dime.vox_offset = 0.0;
	hdr.dime.funused1 = 0.0;
	hdr.dime.funused2 = 0.0;
	hdr.dime.funused3 = 0.0;
	hdr.dime.cal_max = 0.0;
	hdr.dime.cal_min = 0.0;
	hdr.dime.datatype = -1;

	for(i=1;i<=8;i++)
		if(!strcmp("SHORT",DataTypes[i]))	// short int (16-bit, 0-65536)
		{
			hdr.dime.datatype = (1<<(i-1));
			hdr.dime.bitpix = DataTypeSizes[i];
			break;
		}

	if((fp=fopen(fileName,"w"))==0)
	{
		printf("unable to create: %s\n",fileName);
		exit(0);
	}
	hdr.dime.dim[0] = 4; /* all Analyze images are taken as 4 dimensional */
	hdr.hk.regular = 'r';
	hdr.hk.sizeof_hdr = sizeof(struct dsr);
	hdr.dime.dim[1] = m; /* slice width in pixels */
	hdr.dime.dim[2] = n; /* slice height in pixels */
	hdr.dime.dim[3] = 1; /* volume depth in slices */
	hdr.dime.dim[4] = 1; /* number of volumes per file */
	hdr.dime.glmax = 256*256; /* maximum voxel value */
	hdr.dime.glmin = 0; /* minimum voxel value */
/* Set the voxel dimension fields:
A value of 0.0 for these fields implies that the value is unknown.
Change these values to what is appropriate for your data
or pass additional command line arguments */
hdr.dime.pixdim[1] = 1.0; /* voxel x dimension */
hdr.dime.pixdim[2] = 1.0; /* voxel y dimension */
hdr.dime.pixdim[3] = 1.0; /* pixel z dimension, slice thickness */
/* Assume zero offset in .img file, byte at which pixel
data starts in the image file */
hdr.dime.vox_offset = 0.0;
/* Planar Orientation; */
/* Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
Movie flag ON: 3 = transverse, 4 = coronal, 5 = sagittal */
hdr.hist.orient = 0;
/* up to 3 characters for the voxels units label; i.e. mm., um., cm. */
//strcpy(hdr.dime.vox_units," ");		// commented by Igor
/* up to 7 characters for the calibration units label; i.e. HU */
// strcpy(hdr.dime.cal_units," ");		// commented by Igor
/* Calibration maximum and minimum values;
values of 0.0 for both fields imply that no
calibration max and min values are used */
hdr.dime.cal_max = 0.0;
hdr.dime.cal_min = 0.0;
fwrite(&hdr,sizeof(struct dsr),1,fp);
fclose(fp);
}

//#####################################################################

void FileIO_2D::make_header( const DoubleArray2D& A,
							 const string& str,
							 const long& bx, const long& by ) /* file x y z t datatype max min */
{
	ostringstream outs;
    char fileName[256];

    outs.str("");
	outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

	int m = A.getIndex1Size() - 2*bx;
	int n = A.getIndex2Size() - 2*by;


	int i;
	struct dsr hdr;
	FILE *fp;
	static char DataTypes[9][12] = {"UNKNOWN", "BINARY", "CHAR", "SHORT", "INT","FLOAT", "COMPLEX", "DOUBLE", "RGB"};
	static int DataTypeSizes[9] = {0,1,8,16,32,32,64,64,24};

	memset(&hdr,0, sizeof(struct dsr));
	
	for(i=0;i<8;i++)
		hdr.dime.pixdim[i] = 0.0;

	hdr.dime.vox_offset = 0.0;
	hdr.dime.funused1 = 0.0;
	hdr.dime.funused2 = 0.0;
	hdr.dime.funused3 = 0.0;
	hdr.dime.cal_max = 0.0;
	hdr.dime.cal_min = 0.0;
	hdr.dime.datatype = -1;

	for(i=1;i<=8;i++)
		if(!strcmp("SHORT",DataTypes[i]))	// short int (16-bit, 0-65536)
		{
			hdr.dime.datatype = (1<<(i-1));
			hdr.dime.bitpix = DataTypeSizes[i];
			break;
		}

	if((fp=fopen(fileName,"w"))==0)
	{
		printf("unable to create: %s\n",fileName);
		exit(0);
	}
	hdr.dime.dim[0] = 4; /* all Analyze images are taken as 4 dimensional */
	hdr.hk.regular = 'r';
	hdr.hk.sizeof_hdr = sizeof(struct dsr);
	hdr.dime.dim[1] = m; /* slice width in pixels */
	hdr.dime.dim[2] = n; /* slice height in pixels */
	hdr.dime.dim[3] = 1; /* volume depth in slices */
	hdr.dime.dim[4] = 1; /* number of volumes per file */
	hdr.dime.glmax = 256*256; /* maximum voxel value */
	hdr.dime.glmin = 0; /* minimum voxel value */
/* Set the voxel dimension fields:
A value of 0.0 for these fields implies that the value is unknown.
Change these values to what is appropriate for your data
or pass additional command line arguments */
hdr.dime.pixdim[1] = 1.0; /* voxel x dimension */
hdr.dime.pixdim[2] = 1.0; /* voxel y dimension */
hdr.dime.pixdim[3] = 1.0; /* pixel z dimension, slice thickness */
/* Assume zero offset in .img file, byte at which pixel
data starts in the image file */
hdr.dime.vox_offset = 0.0;
/* Planar Orientation; */
/* Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
Movie flag ON: 3 = transverse, 4 = coronal, 5 = sagittal */
hdr.hist.orient = 0;
/* up to 3 characters for the voxels units label; i.e. mm., um., cm. */
//strcpy(hdr.dime.vox_units," ");		// commented by Igor
/* up to 7 characters for the calibration units label; i.e. HU */
// strcpy(hdr.dime.cal_units," ");		// commented by Igor
/* Calibration maximum and minimum values;
values of 0.0 for both fields imply that no
calibration max and min values are used */
hdr.dime.cal_max = 0.0;
hdr.dime.cal_min = 0.0;
fwrite(&hdr,sizeof(struct dsr),1,fp);
fclose(fp);
}

//#####################################################################

//#####################################################################

void FileIO_2D::write_bin_ssi( const DoubleArray2D& A, 
							   const string& str,
							   long factor )
{	//
	// This routine does NOT write a header file.
	//
	ostringstream outs;
	outs << str << ".img";
	string filename = outs.str();

	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	signed short num;
	
	cout << endl << "write_bin_ssi() : ";
	findMaxMin( A );

	cout << "WRITE: signed short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			num = (signed short)( A(i,j)*factor );
			outfile.write((char*)&num,sizeof(signed short));
		}
	}
}

//#####################################################################

void FileIO_2D::write_bin_float( const DoubleArray2D& A, 
							   const string& str )
{	//
	// This routine does NOT write a header file.
	//
	ostringstream outs;
	outs << str << ".raw";
	string filename = outs.str();

	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	float num;
	
	cout << endl << "write_bin_float() : ";
	findMaxMin( A );

	cout << "WRITE: float: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j;
	for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
	{
		for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
		{
			num = (float)( A(i,j) );
			outfile.write((char*)&num,sizeof(float));
		}
	}
}
