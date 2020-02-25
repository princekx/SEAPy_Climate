
export PATH=${PATH}:.
export CC=/usr/bin/gcc
export FC=/usr/bin/gfortran
export ARFLAGS=
export NETCDF=/project/ukmo/rhel7/fortran/opt/gfortran/packages/gnu/8.1.0/netcdf/4.6.1/

module purge
module load gcc
master -build -i=linux -f=linux
make utils


Documentation at 
http://www.met.reading.ac.uk/~dispersion/track/docs.html
