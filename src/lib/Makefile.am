# src/lib/Makefile.am
# by Se-Heon Oh (KASI/ICRAR)

# cfitsio
DIR_CFITSIO = /mnt/g0/research/packages/cfitsio/cfitsio
# gsl
DIR_GSL = /mnt/g0/research/packages/gsl/gsl-1.16/.libs

# additional include paths necessary to compile 2dbat
AM_CFLAGS= -O3 -DMPI -fopenmp -m64 -g -Bstatic -I$(top_srcdir)/src/lib -I$(DIR_CFITSIO) -I$(DIR_GSL)

# the library names to build (note we are building static libs only)
lib_LIBRARIES = lib2dbat.a

# where to install the headers on the system
lib2dbat_adir = $(includedir)/lib

# the list of header files that belong to the library (to be installed later)
lib2dbat_a_HEADERS = 2dbat.2dmaps.h 2dbat.global_params.h 2dbat.multinest.h 2dbat.cfitsio.h 2dbat.gsl.h 2dbat.priors.h 2dbat.einastofit.h 2dbat.main.h 2dbat.sort.h 2dbat.ellipsefit.h 2dbat.memory.h 2dbat.trfit.h 2dbat.etc.h 2dbat.mpi.h 2dbat.gfit.h 2dbat.mpi_datatype.h

# the sources to add to the library and to add to the source distribution
lib2dbat_a_SOURCES= $(libmylib_a_HEADERS) 2dbat.2dmaps.c 2dbat.etc.c 2dbat.mpi.c 2dbat.trfit.c 2dbat.gfit.c 2dbat.multinest.c 2dbat.einastofit.c 2dbat.gsl.c 2dbat.priors.c 2dbat.ellipsefit.c 2dbat.memory.c 2dbat.sort.c 

# end of line
