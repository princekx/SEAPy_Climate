#!/usr/bin/csh

#source setup.TRACK

set MODEL = $1
set SEASON = $2
set FIELD = $3

set EXT = time

set DIRO=/data/tethys/kih/AMIP/${MODEL}/NH_${SEASON}/TSERIES/YEARS
set TRACK=/users/nutis/kih/TRACK27

cd $DIRO

cat ${TRACK}/indat/grid.dat user_tavg.${FIELD}_* > tot_tavg.${FIELD}.dat
cat ${TRACK}/indat/grid.dat user_var.${FIELD}_* > tot_var.${FIELD}.dat
cat ${TRACK}/indat/grid.dat user_varfil.${FIELD}_* > tot_varfil.${FIELD}.dat

set LIST=`ls tot_*.${FIELD}.dat`

cd ${TRACK}

foreach var ($LIST)

   \rm mean.ss

   echo $var

   sed -e "s@user\.dat@${DIRO}\/${var}@" mean.in > mean.ss
   bin/track.${EXT} < mean.ss

   cat indat/grid1.dat outdat/user_tavg.${EXT} > ${DIRO}/avg_$var
   \rm outdat/user_tavg.${EXT}

   \rm mean.ss

   set NEW = `echo $var | sed -e 's/\.dat//'`

   converters/bin2nc ${DIRO}/avg_$var ${DIRO}/avg_$NEW
end

