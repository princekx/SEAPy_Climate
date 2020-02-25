#!/usr/bin/csh

#source setup.TRACK
#setenv LD_LIBRARY_PATH ${HOME}/local/lib/SUNWspro:${LD_LIBRARY_PATH}

set FIELD = $1

set SEASON = $2

set MODEL = $3

set EXT=$4

set OUT=$5

set TRACK = /users/nutis/kih/R_TRACK

set DIR = /data/cadfael2/atmos/kih/DATA/ECHAM5/${SEASON}

cd $DIR

#set LIST=`ls ${MODEL}_${FIELD}_${SEASON}*.utf`
#set LIST=`ls EH5_AMIP*${FIELD}*.nc`
set LIST=`ls EH5*SVO*.nc`
#set LIST="EH5_OMmam1974.nc"

#set LIST=`ls  *fill.dat`

echo $LIST


cd ${TRACK}



foreach filenm ($LIST)

   \rm $OUT

   set STUB=`echo $filenm | sed -e 's/\.utf//'`

   echo $STUB

   if !( -f ${DIR}/${STUB}_filt.dat) then

      sed -e "s/mslp\.dat/${filenm}/;s/DJF/${SEASON}/" indat/SPACE_FILTER.echam5.in > $OUT

#      nice +19 bin/track.$EXT < $OUT > /dev/null
      bin/track.$EXT < $OUT > /dev/null

      set STUB=`echo $filenm | sed -e 's/\.nc//'`

      echo ${STUB}_filt.dat

      \mv outdat/specfil.${EXT}_band000 ${DIR}/${STUB}_plw.dat
      \mv outdat/specfil.${EXT}_band001 ${DIR}/${STUB}_filt.dat

      endif

end
