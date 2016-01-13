//
//#####################################################################
//						Grid3D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// Version: March 15, 2005
//
//#####################################################################
//

#ifndef __Grid3D__
#define __Grid3D__


class Grid3D
{
public :

	double xMin;	// Computational Region is:
	double xMax;	//  [xMin,xMax]x[yMin,yMax]x[zMin,zMax]	
	double yMin;
	double yMax;
	double zMin;
	double zMax;
	long   m;       // Number of Points in the x direction
	long   n;       // Number of Points in the y direction
	long   p;		// Number of Points in the z direction

	long   bxi;		// Number of points the image is extended at the time of input.
	long   byi;		// These points will be cropped before output.
	long   bzi;		// These are usually set to 0 (no extension or cropping).
	long   bxo;		// Number of points the image is cropped at the time of output.
	long   byo;		// These points will be cropped before output.
	long   bzo;		// These are usually set to 0 (no extension or cropping).


	long   w;		// Number of Points used for the width of the boundary of the domain
	long   gamma;	// the Width of the boundary of the tube (for Local Level Set)

	double dx;      // mesh width in x direction
	double dy;      // mesh width in y direction
	double dz;      // mesh width in z direction

	double deltaWidth;

	Grid3D() {	bxi = 0; byi = 0; bzi = 0;
				bxo = 0; byo = 0; bzo = 0;	}
};

#endif
