#ifndef __FileIO_3D__
#define __FileIO_3D__

#include <vector>
#include "DoubleArray1D.h"
#include "DoubleArray3D.h"
#include "Grid3D.h"

typedef vector<    int    > vector1Dint;
typedef vector<vector1Dint> vector2Dint;
typedef vector<vector2Dint> vector3Dint;

//
//#####################################################################
//						FileIO_3D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: Nov. 20, 2006
//
//#####################################################################
//


class FileIO_3D
{
public:

	// ASCII Info:
	void getDAT3DInfo( string filename, long& m, long& n, long& p );

	// ASCII READ:
	void readDAT3D( DoubleArray3D& A, string filename );
	void readDAT3D( DoubleArray3D& A, long m, long n, long p, string filename );
	void readDAT1D( DoubleArray1D& A, string filename );
	void readDAT1D( vector<int>&   A, string filename );

	// ASCII WRITE:
	void write_ascii( const vector<int>&    A, const string& str );
	void write_ascii( const vector<double>& A, const string& str, const long& stepCount );
	void write_ascii( const vector3Dint&    A, const string& str );
	void write_ascii( const DoubleArray3D&  A, const string& str, const long& stepCount );
	void write_ascii( const DoubleArray3D&  A, const string& str, const long& stepCount, const long& w );
	void write_ascii( const DoubleArray3D&  A, const string& str );
	void write_ascii( const DoubleArray1D&  A, const string& str );	
	void write_ascii( const DoubleArray1D&  A, const string& str, const long& stepCount );

	void writeDAT3D( const DoubleArray3D& A, string filename );	// slow write, write_ascii() is preferred.

	void writeImageSeg3D( const DoubleArray3D& u, long n );


	// BINARY READ:
	void read_Analyze_hdr(            string filename, long& m, long& n, long& p );

	void read_bin_header_info(        string filename, DoubleArray3D& A );	// if file contains the header info in certain form; reads in doubles
	void read_bin(                    string filename, DoubleArray3D& A );	// double
	void read_bin_float(              string filename, DoubleArray3D& A );	// float
	void read_bin_unsigned_char(      string filename, DoubleArray3D& A );	// unsigned char
	void read_bin_signed_short_int(   string filename, DoubleArray3D& A, bool byteswap );	// signed short int
	void read_bin_unsigned_short_int( string filename, DoubleArray3D& A, bool byteswap );	// unsigned short int
	void read_bin_int(                string filename, DoubleArray3D& A );	// int

	// BINARY WRITE:
	void write_bin(     const DoubleArray3D& A, const string& str );				// double
	void write_bin_int( DoubleArray3D& A, string filename );						// int
	void write_bin_uc(  const DoubleArray3D& A, const string& str );				// unsigned char
	void write_bin_usi( const DoubleArray3D& A, const string& str );				// unsigned short int
	void write_bin_usi( const DoubleArray3D& A, const string& str, long factor );	// unsigned short int
	void write_bin_usi( const DoubleArray3D& A, const string& str, const long& stepCount, long factor );
	void write_bin_usi( const DoubleArray3D& A, const string& str, long factor,
					    const long bx, const long by, const long bz );
	void write_bin_usi( const DoubleArray3D& A, const string& str, const long& stepCount, long factor,
						const long bx, const long by, const long bz );
	void write_bin_ssi( const DoubleArray3D& A, const string& str, long factor );	// signed  short int
	
	// HEADER:
	void make_header_double( const DoubleArray3D& A, const string& str );
	void make_header( const DoubleArray3D& A, const string& str );
	void make_header( const DoubleArray3D& A, const string& str,
					  const long& bx, const long& by, const long& bz );
	void make_header(    const DoubleArray3D& A, const string& str, char *dtype );
	void make_header_uc( const DoubleArray3D& A, const string& str );
	void make_header_ssi( const DoubleArray3D& A, const string& str );
};

#endif

