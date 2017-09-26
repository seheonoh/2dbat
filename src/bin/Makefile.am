# src/bin/Makefile.am
# by Se-Heon Oh (KASI/ICRAR)

#-----------------------------------------------------
# additional libraries linked to 2dbat
#-----------------------------------------------------

#++++++++++
# 1/7. multinest
DIR_MULTINEST = /mnt/g0/research/packages/multinest/MultiNest_v3.7

#++++++++++
# 2/7. cfitsio
DIR_CFITSIO = /mnt/g0/research/packages/cfitsio/cfitsio

#++++++++++
# 3/7. gsl
DIR_GSL = /mnt/g0/research/packages/gsl/gsl-1.16/.libs

#++++++++++
# 4/7. lapack
DIR_LAPACK = /usr/lib/lapack

#++++++++++
# 5/7. blas
DIR_BLAS = /usr/lib/libblas

#++++++++++
# 6/7. cblas
DIR_CBLAS = /usr/lib/atlas-base

#++++++++++
# 7/7. atlas
DIR_ATLAS = /usr/lib/atlas-base


AM_CFLAGS= -O3 -DMPI -fopenmp -m64 -g -Bstatic -I$(top_srcdir)/src/lib -I$(DIR_CFITSIO) -I$(DIR_GSL)

# the program to build (the names of the final binaries)
bin_PROGRAMS= 2dbat

# list of sources for the '2dbat' binary
2dbat_SOURCES= 2dbat.c

# the additional libraries needed to link 2dbat
2dbat_LDADD= -L$(top_srcdir)/src/lib -l2dbat -L$(DIR_MULTINEST) -lnest3 -L$(DIR_LAPACK) -llapack -L$(DIR_CFITSIO) -lcfitsio -L$(DIR_GSL) -lgsl -L$(DIR_BLAS) -lblas -L$(DIR_CBLAS) -lcblas -L$(DIR_ATLAS) -latlas -lm

# end of line
