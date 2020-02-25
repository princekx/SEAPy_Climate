#!/bin/csh

set me = `readlink -f ${0}` ; set GFS_SCRIPTS = `dirname $me`
source ${GFS_SCRIPTS}/common_settings.csh

set DATE=$1
set TIME=$2
set HEMI=$3
set FIELD=$4

set TRACK=${HOMEDIR}/TRACK-1.4.3
set GFSIN=${TRACK}/GFS
set BINDR=${TRACK}/utils/bin

set DATETM=${DATE}${TIME}
set DDIR=${DATETM}_$FIELD

set GFSDR=${GFSDIR}/${DATETM}/$HEMI

mkdir ${GFSDR}/${DDIR}
cd ${GFSDR}/$DDIR
ln -s ${TRACK}/GFS/hist.dat hist.dat

\rm eps.in

if( "$FIELD" == "mslp" )then
    sed -e "s/DATETM/${DATETM}/g;s/NH/${HEMI}/" ${GFSIN}/eps_cluster_mslp.in > eps.in
else if( "$FIELD" == "vor" )then
   if( "$HEMI" == "NH" )then
      sed -e "s/DATETM/${DATETM}/g;s/NH/${HEMI}/" ${GFSIN}/eps_cluster_vor.in > eps.in
   else if ( "$HEMI" == "SH" )then
      sed -e "s/DATETM/${DATETM}/g;s/NH/${HEMI}/;s/pos/neg/" ${GFSIN}/eps_cluster_vor.in > eps.in
   endif
endif

${BINDR}/eps < eps.in
