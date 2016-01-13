
# makefile
#
# Created at:    August 7, 2007
# Last modified: August 7, 2007

# BEGINNING OF MAKEFILE
TARGET  = fluid3d 

# change the following 'HOMEDIR' to your local path
#HOMEDIR = /home1/.....
HOMEDIR = /curr/boyuan/Yanovsky_Vese_NonlinFluid3D

OBJDIR  = $(HOMEDIR)
LIBDIR  = $(HOMEDIR)
CURRDIR = $(HOMEDIR)
OUTPUT  = $(CURRDIR)

SRCDIR  = $(CURRDIR):$(OBJDIR)

vpath %.cpp $(SRCDIR)

# put all the used .cpp file here:
ALLSRC = BC_2D.cpp BC_3D.cpp  common_routines.cpp  Convolution2D.cpp  Convolution3D.cpp  FileIO_2D.cpp  FileIO_3D.cpp  INIT_img3D.cpp  J_Reg3D.cpp  MAIN.cpp  MGsystem3D.cpp  Norms.cpp  StationaryNS3D.cpp  timer.cpp  
ALLOBJ = BC_2D.o BC_3D.o  common_routines.o  Convolution2D.o  Convolution3D.o  FileIO_2D.o  FileIO_3D.o  INIT_img3D.o  J_Reg3D.o  MAIN.o  MGsystem3D.o  Norms.o  StationaryNS3D.o  timer.o  

#ALLOBJ = $(ALLSRC:%.cpp=$(OUTPUT)/%.o)

#FFTW3 = /u/local/lib/fftw-3.1.2-i386-gcc
#FFTW3 = /u/local/lib/fftw-3.1.2-x86_64-gcc
#FFTW3 = /curr/boyuan/Yanovsky_Vese_NonlinFluid3D/lib/fftw-3.2.2-x86_64-gcc
FFTW3 = /curr/boyuan/Yanovsky_Vese_NonlinFluid3D


INCLUDEDIR = -I$(CURRDIR) -I$(OBJDIR) -I$(FFTW3)/include
CPPFLAGS   = -O3 $(INCLUDEDIR)

LDFLAG = -lm -L$(FFTW3)/lib -lfftw3  # -L$(LIBDIR)  
#/usr/lib/libfftw.a

CL = g++

%.o : %.cpp 
	$(CL) $(CPPFLAGS) -c $< -o $@

$(TARGET) : $(ALLOBJ)
	$(CL) $(CPPFLAGS) -o $@ $^ $(LDFLAG)

depend : $(ALLSRC)
	makedepend $(CPPFLAGS) $^


.PHONY : clean
clean:
	rm $(TARGET) $(ALLOBJ)


# DO NOT DELETE



