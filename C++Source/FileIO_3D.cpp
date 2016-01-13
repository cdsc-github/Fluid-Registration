#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <fstream>
#include <math.h>

//
//#####################################################################
//						FileIO_3D.cpp
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: Nov. 20, 2006
//
//#####################################################################
//

#include "FileIO_3D.h"
#include "dbh.h"
#include "common_routines.h"


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

//#####################################################################

	// To read a binary value from an input file stream 
	//    (ifstream) object, use the read() member function. 
	// read() takes two parameters: char * and long, which hold 
	//    the address of the buffer to which the value is written 
	//    and its size, respectively. 
//
//#####################################################################
//							BINARY READ
//#####################################################################
//
void FileIO_3D::read_bin_header_info( string filename, DoubleArray3D& A )
{
	//
	// This routine reads binary files that have a header info included.
	// It checks if the dimensions of A are compatible with
	// the dimensions provided by the header portion.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	int x0,y0,z0,xf,yf,zf;
	infile.read(reinterpret_cast <char*> (&x0),sizeof(int));
	infile.read(reinterpret_cast <char*> (&y0),sizeof(int));
	infile.read(reinterpret_cast <char*> (&z0),sizeof(int));
	cout << "x0 = " << x0 << " y0 = " << y0 << " z0 = " << z0 << endl;
	
	infile.read(reinterpret_cast <char*> (&xf),sizeof(int));
	infile.read(reinterpret_cast <char*> (&yf),sizeof(int));
	infile.read(reinterpret_cast <char*> (&zf),sizeof(int));
	cout << "xf = " << xf << " yf = " << yf << " zf = " << zf << endl;
	
	if( A.getIndex1End()+1 != xf || A.getIndex2End()+1 != yf || A.getIndex3End()+1 != zf )
	{
		cout << "Grid Dimensions mismatch!!!" << endl;
		exit(1);
	}
	
	cout << "READ: double: " << sizeof(double) << " byte ( " << sizeof(double) * 8 << " bit ) data reading." << endl;

	long i;  long j;  long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
			    infile.read(reinterpret_cast <char*> (&A(i,j,k)),sizeof(double));
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin( string filename, DoubleArray3D& A )
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
	
	cout << "READ: double: " << sizeof(double) << " byte ( " << sizeof(double) * 8 << " bit ) data reading." << endl;

	long i;  long j;  long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
			    infile.read(reinterpret_cast <char*> (&A(i,j,k)),sizeof(double));
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin_float( string filename, DoubleArray3D& A )
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

	long i;  long j;  long k;
	float fl;

	cout << "READ: float: " << sizeof(float) << " byte ( " << sizeof(float) * 8 << " bit ) data reading." << endl;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				// The following code retrieves a float from a file
				// and copies it to fl:
			    infile.read(reinterpret_cast <char*> (&fl),sizeof(float));
				// Convert float fl to double and assign to A(i,j,k):
				A(i,j,k) = (double)fl;
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin_unsigned_char( string filename, DoubleArray3D& A )
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
	
	long i;  long j;  long k;
	unsigned char uc;

	cout << "READ: unsigned char: " << sizeof(uc) << " byte ( " << sizeof(uc) * 8 << " bit ) data reading." << endl;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				// The following code retrieves an unsigned char from a file
				// and copies it to uc:
			    infile.read(reinterpret_cast <char*> (&uc),sizeof(unsigned char));
				
				// Convert unsigned char uc to double and assign to A(i,j,k):
				A(i,j,k) = uc;
				//if( A(i,j,k) > 0)	{	display( A, i, j, k );	}
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin_unsigned_short_int( string filename, DoubleArray3D& A, bool byteswap )
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
	
	long i;  long j;  long k;
	unsigned short int usi;
	
	cout << "READ: unsigned short int: " << sizeof(usi) << " byte ( " << sizeof(usi) * 8 << " bit ) data." << endl;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
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
				A(i,j,k) = usi;
				//if( A(i,j,k) > 0)	{	display( A, i, j, k );	}
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin_signed_short_int( string filename, DoubleArray3D& A, bool byteswap )
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
	
	long i;  long j;  long k;
	signed short int ssi;
	
	cout << "READ: signed short int: " << sizeof(ssi) << " byte ( " << sizeof(ssi) * 8 << " bit ) data." << endl;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
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
				A(i,j,k) = ssi;

				/*if( A(i,j,k) != 0)	
				{	display( A, i, j, k );	
					exit(1);	
				}*/
			}
		}
	}
}

//#####################################################################

void FileIO_3D::read_bin_int( string filename, DoubleArray3D& A )
{
	//
	// This routine reads binary files that do NOT have 
	// a header info included.
	//

	fstream infile( filename.c_str(), ios::in | ios::binary );
	if(infile.fail())
	{
		cout << "Error reading " << filename << "!!!" << endl;
		exit(-1);
	}
	
	long i;  long j;  long k;
	int num;

	cout << "READ: int: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				// The following code retrieves an int from a file
				// and copies it to num:
			    infile.read(reinterpret_cast <char*> (&num),sizeof(int));

				// Convert integer num to double and assign to A(i,j,k):
				A(i,j,k) = num;
				//if( A(i,j,k) > 0)	{	display( A, i, j, k );	}
			}
		}
	}
}
//
//#####################################################################
//							BINARY WRITE
//#####################################################################
//
void FileIO_3D::write_bin( const DoubleArray3D& A, 
						   const string& str )
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

	/*
	distw_file.write(reinterpret_cast <char*> (&grid_x),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_y),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_z),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_x),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_y),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_z),sizeof(int));
	*/

	double num;

	cout << "WRITE: double: " << sizeof(double) << " byte ( " << sizeof(double) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = A(i,j,k);
				outfile.write(reinterpret_cast <char*> (&num),sizeof(double));
			}
		}
	}
	make_header_double( A, str );
}

void FileIO_3D::write_bin_uc( const DoubleArray3D& A, const string& str )
{
	//
	// This routine WRITES a binary file with NO header.
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

	unsigned char num;
	
	cout << endl << "In write_bin_uc() : ";
	findMaxMin( A );

	cout << "WRITE: unsigned char: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = (unsigned char)( A(i,j,k) );
				outfile.write((char*)&num,sizeof(unsigned char));
			}
		}
	}
	make_header_uc( A, str );
}

//#####################################################################

void FileIO_3D::write_bin_int( DoubleArray3D& A, string filename )
{
	//
	// This routine WRITES a binary file with NO header.
	//

	fstream outfile( filename.c_str(), ios::out | ios::binary );
	if(outfile.fail())
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	/*
	distw_file.write(reinterpret_cast <char*> (&grid_x),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_y),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_z),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_x),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_y),sizeof(int));
	distw_file.write(reinterpret_cast <char*> (&grid_z),sizeof(int));
	*/

	int num;

	cout << "WRITE: int: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = round_to_int( A(i,j,k) );
				outfile.write(reinterpret_cast <char*> (&num),sizeof(int));
			}
		}
	}
}

//#####################################################################

void FileIO_3D::write_bin_usi( const DoubleArray3D& A, 
							   const string& str )
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
	
	cout << endl << "In write_bin_usi() : ";
	findMaxMin( A );

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = (unsigned short)( A(i,j,k) );
				//if (A(i,j,k) > 260.0 )	{ cout << A(i,j,k) << " " << num << endl;	}
				outfile.write((char*)&num,sizeof(unsigned short));
			}
		}
	}

	make_header( A, str );
}

//#####################################################################

void FileIO_3D::write_bin_usi( const DoubleArray3D& A, 
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
	
	cout << endl << "write_bin_usi() : ";
	findMaxMin( A );

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = (unsigned short)( A(i,j,k)*factor );
				outfile.write((char*)&num,sizeof(unsigned short));
			}
		}
	}

	make_header( A, str );
}

//#####################################################################

void FileIO_3D::write_bin_usi( const DoubleArray3D& A, 
							   const string& str,
							   long factor,
							   const long bx, const long by, const long bz )
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
	
	//cout << endl << "In write_bin_usi() : ";
	//findMaxMin( A );

	cout << "WRITE: unsigned short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin() + bz; k <= A.getIndex3End() - bz; k++)
	{
		for(j = A.getIndex2Begin() + by; j <= A.getIndex2End() - by; j++)
		{
			for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
			{
				num = (unsigned short)( A(i,j,k)*factor );
				/*if( i == 111 && j == 111 && k == 111 )
				{
					display( A, i, j, k );
					cout << A(i,j,k)*factor << " " << num << endl << endl;
				}*/
				outfile.write((char*)&num,sizeof(unsigned short));
			}
		}
	}

	make_header( A, str, bx, by, bz );
}

//#####################################################################

void FileIO_3D::write_bin_usi( const DoubleArray3D& A, 
							   const string& str,
							   const long& stepCount,
							   long factor,
							   const long bx, const long by, const long bz )
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

	long i; long j; long k;
	for(k = A.getIndex3Begin() + bz; k <= A.getIndex3End() - bz; k++)
	{
		for(j = A.getIndex2Begin() + by; j <= A.getIndex2End() - by; j++)
		{
			for(i = A.getIndex1Begin() + bx; i <= A.getIndex1End() - bx; i++)
			{
				num = (unsigned short)( A(i,j,k)*factor );
				outfile.write((char*)&num,sizeof(unsigned short));
			}
		}
	}

	make_header( A, basename, bx, by, bz );
}

//#####################################################################

void FileIO_3D::write_bin_usi( const DoubleArray3D& A, 
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

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = (unsigned short)( A(i,j,k)*factor );
				outfile.write((char*)&num,sizeof(unsigned short));
			}
		}
	}

	make_header( A, basename );
}

//#####################################################################

void FileIO_3D::write_bin_ssi( const DoubleArray3D& A, 
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

	signed short num;
	
	cout << endl << "write_bin_ssi() : ";
	findMaxMin( A );

	cout << "WRITE: signed short: " << sizeof(num) << " byte ( " << sizeof(num) * 8 << " bit ) data." << endl;

	long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				num = (signed short)( A(i,j,k)*factor );
				outfile.write((char*)&num,sizeof(signed short));
			}
		}
	}

	make_header_ssi( A, str );
}

//
//#####################################################################
//						  GET IMAGE INFO
//#####################################################################
//
void FileIO_3D::getDAT3DInfo( string filename, long& m, long& n, long& p)
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	infile >> m >> n >> p;		// Read the dimension of the PGM image
}
//
//#####################################################################
//							ASCII READ
//#####################################################################
//
void FileIO_3D::readDAT3D( DoubleArray3D& A, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long m; long n; long p;

	infile >> m >> n >> p;		// Read the dimension of the DAT image

	long i;  long j;  long k;
	for(k = 0; k < p; k++)
	{
		for(j = 0; j < n; j++)
		{
			for(i = 0; i < m; i++)
			{
				infile >> A(i,j,k);
			}
		}
	}
	infile.close();
}
//
//#####################################################################
//
void FileIO_3D::readDAT1D( DoubleArray1D& A, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long m = A.getIndex1Size();
	long i;
	for(i = 0; i < m; i++)
	{
		infile >> A(i);
	}

	infile.close();
}
//
//#####################################################################
//
void FileIO_3D::readDAT1D( vector<int>& A, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long m = long(A.size());
	long i;
	for(i = 0; i < m; i++)
	{
		infile >> A[i];
	}

	infile.close();
}

//#####################################################################

/*void FileIO_3D::readDAT3D( DoubleArray3D& A, long m, long n, long p, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	long i;  long j;  long k;
	for(k = 0; k < p; k++)
	{
		for(j = 0; j < n; j++)
		{
			for(i = 0; i < m; i++)
			{
				infile >> A(i,j,k);
			}
		}
	}
	infile.close();
}*/

void FileIO_3D::readDAT3D( DoubleArray3D& A, long m, long n, long p, string filename )
{
	ifstream infile( filename.c_str() );
	if(!infile)
	{
		cout << "Error reading " << filename<< "!!!" << endl;
		exit(-1);
	}

	//float value;
	long count = 0;
	long i;  long j;  long k;
	for(k = 0; k < p; k++)
	{
		for(j = 0; j < n; j++)
		{
			for(i = 0; i < m; i++)
			{
				infile >> A(i,j,k);
			}
		}
	}
	cout << endl << endl << count;
	infile.close();
}
//
//#####################################################################
//							ASCII WRITE
//#####################################################################
//
void FileIO_3D::write_ascii( const vector3Dint& A, const string& str )
{
	cout << "write_ascii( const vector3Dint& A, const string& str )" << endl;
	cout << "Written by James but not tested by me" << endl;
	ostringstream outs;
	outs << str << ".dat";
	string filename = outs.str();

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

    long i; long j; long k;
	for(k = 0; k < long(A[0][0].size()); k++)
	{
		for(j = 0; j < long(A[0].size()); j++)
		{
			for(i = 0; i < long(A.size()); i++)
			{
				outfile <<  setw(5) << A[i][j][k] << " ";
				cout << A[i][j][k] << endl;
			}
			outfile << endl;
		}
		outfile << endl;
	}
	outfile.close();
}
//
//#####################################################################
//
void FileIO_3D::write_ascii( const DoubleArray3D& A,
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
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

	//
	//  Open and then write to a file
	//
    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
		printf( "The file %s could not be  opened\n",fileName);
		return;
    }

    long i; long j; long k;
	//
	//  Output the data.
	//

	// x is the fastest-running variable:
	//
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++) 
    {
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				fprintf(dataFile,"%-10.5e ",A(i,j,k));
			}
			fprintf(dataFile,"\n");
		}
		fprintf(dataFile,"\n");
	}
	
	fclose(dataFile);
}

//#####################################################################

void FileIO_3D::write_ascii( const DoubleArray3D& A,
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
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

	//
	//  Open and then write to a file
	//
    FILE* dataFile;

	if( (dataFile = fopen(fileName, "w+" )) == NULL )
    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
		printf( "The file %s could not be  opened\n",fileName);
		return;
    }

    long i; long j; long k;
	//
	//  Output the data.
	//

	// x is the fastest-running variable:
	//
	for(k = w; k <= A.getIndex3End()-w; k++) 
    {
		for(j = w; j <= A.getIndex2End()-w; j++)
		{
			for(i = w; i <= A.getIndex1End()-w; i++)
			{
				fprintf(dataFile,"%-10.5e ",A(i,j,k));
			}
			fprintf(dataFile,"\n");
		}
		fprintf(dataFile,"\n");
	}
	
	fclose(dataFile);
}

//#####################################################################

void FileIO_3D::write_ascii( const DoubleArray3D& A,
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

    long i; long j; long k;

	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++) 
    {
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				fprintf(dataFile,"%-10.5e ",A(i,j,k));
			}
			fprintf(dataFile,"\n");
		}
		fprintf(dataFile,"\n");
	}
	
	fclose(dataFile);
}

//#####################################################################

void FileIO_3D::writeDAT3D( const DoubleArray3D& A, string filename )
{
	cout << "The preferred and faster way of writing is to use write_ascii() routines." << endl << endl;

	ofstream outfile( filename.c_str() );
	if(!outfile)
	{
		cout << "Error opening  " << filename << "!!!" << endl;
		exit(-1);
	}

	outfile << A.getIndex1Size() << " " << A.getIndex2Size() << " " << A.getIndex3Size() << endl;

    long i; long j; long k;
	for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
	{
		for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
		{
			for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
			{
				outfile <<  setw(5) << A(i,j,k) << " ";
			}
			outfile << endl;
		}
    }
	outfile.close();
}

//#####################################################################

void FileIO_3D::writeImageSeg3D( const DoubleArray3D& A, long n )
{
	cout << "Use  write_ascii( A, str ) or write_ascii( A, str, number ) with str = \"u_seg\" " << endl;
	cout << "Exiting..." << endl;
	exit(1);
}
//
//#####################################################################
//						WRITE 1D array
//#####################################################################
void FileIO_3D::write_ascii( const vector<int> & A,
							 const string& str )
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
	for( i = 0; i < long(A.size()); i++ )
	{
		outfile << A[i] << endl;;
	}

	outfile.close();
}

//#####################################################################

void FileIO_3D::write_ascii( const vector<double>& A,
							 const string& str, const long& stepCount )
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

	for(i = 0;  i < long(A.size());  i++)
	{	fprintf( dataFile, "%-10.5e \n", A[i] );	}
	
	fclose(dataFile);
}

//#####################################################################

void FileIO_3D::write_ascii( const DoubleArray1D& A,
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

void FileIO_3D::write_ascii( const DoubleArray1D& A,
        const string& str, const long& stepCount )
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << stepCount << ".dat";
    strcpy(fileName,(outs.str()).c_str());
    //strcpy_s( fileName, 256, (outs.str()).c_str() );	// strcpy(fileName,(outs.str()).c_str());

    FILE* dataFile;

    //if( fopen_s( &dataFile, fileName, "w+" ) != 0 )		// if( (dataFile = fopen(fileName, "w+" )) == NULL )
    if( (dataFile = fopen(fileName, "w+" )) == NULL )
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

void FileIO_3D::make_header( const DoubleArray3D& A,
        const string& str ) /* file x y z t datatype max min */
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

    int m = A.getIndex1Size();
    int n = A.getIndex2Size();
    int p = A.getIndex3Size();

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
    hdr.dime.dim[3] = p; /* volume depth in slices */
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

void FileIO_3D::make_header_uc( const DoubleArray3D& A,
        const string& str ) /* file x y z t datatype max min */
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

    int m = A.getIndex1Size();
    int n = A.getIndex2Size();
    int p = A.getIndex3Size();

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
        if(!strcmp("CHAR",DataTypes[i]))	// char (8-bit, 0-255)
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
    hdr.dime.dim[3] = p; /* volume depth in slices */
    hdr.dime.dim[4] = 1; /* number of volumes per file */
    hdr.dime.glmax = 255; /* maximum voxel value */
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

void FileIO_3D::make_header( const DoubleArray3D& A,
        const string& str,
        const long& bx, const long& by, const long& bz ) /* file x y z t datatype max min */
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

    int m = A.getIndex1Size() - 2*bx;
    int n = A.getIndex2Size() - 2*by;
    int p = A.getIndex3Size() - 2*bz;

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
    hdr.dime.dim[3] = p; /* volume depth in slices */
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

void FileIO_3D::make_header_ssi( const DoubleArray3D& A,
        const string& str ) /* file x y z t datatype max min */
{
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

    int m = A.getIndex1Size();
    int n = A.getIndex2Size();
    int p = A.getIndex3Size();

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
    hdr.dime.dim[3] = p; /* volume depth in slices */
    hdr.dime.dim[4] = 1; /* number of volumes per file */
    hdr.dime.glmax = 256*256/2; /* maximum voxel value */
    hdr.dime.glmin = -256*256/2+1; /* minimum voxel value */
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

void FileIO_3D::make_header_double( const DoubleArray3D& A,
        const string& str )
{
    cout << "THIS routine is not functional yet! " << endl;
    ostringstream outs;
    char fileName[256];

    outs.str("");
    outs << str << ".hdr";
    strcpy(fileName,(outs.str()).c_str());

    int m = A.getIndex1Size();
    int n = A.getIndex2Size();
    int p = A.getIndex3Size();


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

    /*for(i=1;i<=8;i++)					// Commented by Igor
      if(!strcmp("DOUBLE",DataTypes[i]))
      {
      hdr.dime.datatype = (1<<(i-1));
      hdr.dime.bitpix = DataTypeSizes[i];
      break;
      }
      */
    hdr.dime.datatype = 8;	// Added by Igor
    hdr.dime.bitpix = 64;	// Added by Igor

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
    hdr.dime.dim[3] = p; /* volume depth in slices */
    hdr.dime.dim[4] = 1; /* number of volumes per file */
    hdr.dime.glmax = 1000; /* maximum voxel value */
    hdr.dime.glmin = -1000; /* minimum voxel value */
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
