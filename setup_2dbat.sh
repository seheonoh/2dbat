#!/bin/bash
# setup_2dbat.sh 
# by Se-Heon Oh (KASI/ICRAR)
# Version v1.0.0

#......................
# PRE-STEP FOR AUTOCONF
#......................

#--------------------------------------------------------
# EDIT ./src/bin/Makefile.am
#--------------------------------------------------------
# 1. MULTINEST DIRECTORY from the user
echo ""
echo "------------------------------------------------------------------------"
echo "... Setup the paths to the libraries required for 2DBAT installation ..."
echo "------------------------------------------------------------------------"
echo ""
default_path="/mnt/g0/research/packages/multinest/MultiNest_v3.7"
read -p "--> 1/7. Enter the absolute path to the directory where MULTINEST is installed <--
    [Check with 'locate libnest3.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/multinest/MultiNest_v3.7' 
    or enter another directory and press 'Enter']: " dir_multinest
dir_multinest="${dir_multinest:-$default_path}"
sed -i "s|^DIR_MULTINEST.*|DIR_MULTINEST = ${dir_multinest}|g" ./src/bin/Makefile.am
echo "    >> DIR_MULTINEST= ${dir_multinest=}"


# 2. CFITSIO DIRECTORY from the user
echo ""
default_path="/mnt/g0/research/packages/cfitsio/cfitsio"
read -p "--> 2/7. Enter the absolute path to the directory where CFITSIO is installed <--
    [Check with 'locate libcfitsio.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/cfitsio/cfitsio' 
    or enter another directory and press 'Enter']: " dir_cfitsio
dir_cfitsio="${dir_cfitsio:-$default_path}"
sed -i "s|^DIR_CFITSIO.*|DIR_CFITSIO = ${dir_cfitsio}|g" ./src/bin/Makefile.am
echo "    >> DIR_CFITSIO= ${dir_cfitsio=}"


# 3. GSL DIRECTORY from the user
echo ""
default_path="/mnt/g0/research/packages/gsl/gsl-1.16/.libs"
read -p "--> 3/7. Enter the absolute path to the directory where GSL is installed <--
    [Check with 'locate libgsl.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/gsl/gsl-1.16/.libs' 
    or enter another directory and press 'Enter']: " dir_gsl
dir_gsl="${dir_gsl:-$default_path}"
sed -i "s|^DIR_GSL.*|DIR_GSL = ${dir_gsl}|g" ./src/bin/Makefile.am
echo "    >> DIR_GSL= ${dir_gsl=}"


# 4. LAPACK DIRECTORY from the user
echo ""
default_path="/usr/lib/lapack"
read -p "--> 4/7. Enter the absolute path to the directory where LAPACK is installed <--
    [Check with 'locate liblapack.a']
    [Press 'Enter' to use the default directory '/usr/lib/lapack' 
    or enter another directory and press 'Enter']: " dir_lapack
dir_lapack="${dir_lapack:-$default_path}"
sed -i "s|^DIR_LAPACK.*|DIR_LAPACK = ${dir_lapack}|g" ./src/bin/Makefile.am
echo "    >> DIR_LAPACK= ${dir_lapack=}"

# 5. BLAS DIRECTORY from the user
echo ""
default_path="/usr/lib/libblas"
read -p "--> 5/7. Enter the absolute path to the directory where BLAS is installed <--
    [Check with 'locate libblas.a']
    [Press 'Enter' to use the default directory '/usr/lib/libblas' 
    or enter another directory and press 'Enter']: " dir_blas
dir_blas="${dir_blas:-$default_path}"
sed -i "s|^DIR_BLAS.*|DIR_BLAS = ${dir_blas}|g" ./src/bin/Makefile.am
echo "    >> DIR_BLAS= ${dir_blas=}"

# 6. CBLAS DIRECTORY from the user
echo ""
default_path="/usr/lib/atlas-base"
read -p "--> 6/7. Enter the absolute path to the directory where CBLAS is installed <--
    [Check with 'locate libcblas.a']
    [Press 'Enter' to use the default directory '/usr/lib/atlas-base' 
    or enter another directory and press 'Enter']: " dir_cblas
dir_cblas="${dir_cblas:-$default_path}"
sed -i "s|^DIR_CBLAS.*|DIR_CBLAS = ${dir_cblas}|g" ./src/bin/Makefile.am
echo "    >> DIR_CBLAS= ${dir_cblas=}"

# 7. ATLAS DIRECTORY from the user
echo ""
default_path="/usr/lib/atlas-base"
read -p "--> 7/7. Enter the absolute path to the directory where ATLAS is installed <--
    [Check with 'locate libatlas.a']
    [Press 'Enter' to use the default directory '/usr/lib/atlas-base' 
    or enter another directory and press 'Enter']: " dir_atlas
dir_atlas="${dir_atlas:-$default_path}"
sed -i "s|^DIR_ATLAS.*|DIR_ATLAS = ${dir_atlas}|g" ./src/bin/Makefile.am
echo "    >> DIR_ATLAS= ${dir_atlas=}"

#--------------------------------------------------------
# EDIT ./src/lib/Makefile.am
#--------------------------------------------------------
# Update the path of the CFITSIO DIRECTORY obtained from the user
sed -i "s|^DIR_CFITSIO.*|DIR_CFITSIO = ${dir_cfitsio}|g" ./src/lib/Makefile.am
# Update the path of the GSL DIRECTORY obtained from the user
sed -i "s|^DIR_GSL.*|DIR_GSL = ${dir_gsl}|g" ./src/lib/Makefile.am

echo ""
echo "-------------------------------------------------"
echo "... The paths to the libraries are configured ..."
echo "-------------------------------------------------"
echo ""

echo ""
echo "----------------------"
echo "... Start AUTOCONF ..."
echo "----------------------"
echo ""
#......................
# AUTOCONF STEP
#......................
aclocal
autoheader
autoconf
#touch NEWS README AUTHORS ChangeLog
#libtoolize --automake --copy --force
automake --foreign --copy --add-missing
./configure

echo ""
echo "---------------------------------------------------------------------------"
echo "... Makefile is configured... Now, type 'make' in the current directory ..."
echo "---------------------------------------------------------------------------"
echo ""


