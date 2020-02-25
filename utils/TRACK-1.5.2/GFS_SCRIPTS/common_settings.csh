#!/bin/csh

# this is used as the base output directory, with GFS appended.
set HOMEDIR=${HOME}
set GFSDIR=/users/kih/TRACK-hg-working/GFS
set VERSION=hg-working

setenv PATH ${PATH}:${HOMEDIR}/TRACK-${VERSION}/GFS_SCRIPTS
setenv LD_LIBRARY_PATH $HOMEDIR/local/LINUX/NETCDF-4.1.2-beta3-snapshot2011020522/lib
# set REMOTE="kih@atlas.nerc-essc.ac.uk"
# set REMOTEDIR=/data/kevin/GFS

# define base thredds URLS to download from
# these are nomads, could also use kryten in theory.
set BASE_D=http://nomads.ncep.noaa.gov:9090/dods/gfs_0p50
set BASE_E=http://nomads.ncep.noaa.gov:9090/dods/gens

set track=bin/track.linux
set master=./master
