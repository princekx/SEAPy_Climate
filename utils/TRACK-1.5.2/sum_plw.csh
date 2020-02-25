#!/usr/bin/csh

#source setup.TRACK

set MODEL = $1
set SEASON = $2
set FIELD = $3

set EXT = $4

set DIRO=/data/atlas/kih/AMIP/${MODEL}/NH_${SEASON}/PLW

set TRACK=${HOME}/S_TRACK

mkdir -p $DIRO

set DATA=/data/scratch03/kih/GISS

cd $DATA

set LOWES=`echo ${SEASON} | tr "[:upper:]" "[:lower:]"`
set LIST=`ls ${MODEL}_${FIELD}_${LOWES}*plw.dat`

#set LIST=`ls ${FIELD}_${LOWES}*plw.dat`


echo $LIST

cd ${TRACK}

foreach var ($LIST)

   \rm mean.plw

   echo $var

   sed -e "s@user\.dat@${DATA}\/${var}@" mean_plw.in > mean.plw
   bin/track.${EXT} < mean.plw 

   mv outdat/user_tavg.${EXT} ${DIRO}/avg_${var}

   \rm mean.plw

end

cat indat/grid_giss.dat ${DIRO}/avg_*${FIELD}_${LOWES}*.dat > ${DIRO}/tot_avg_${FIELD}_${LOWES}_plw.dat

sed -e "s@user\.dat@${DIRO}/tot_avg_${FIELD}_${LOWES}_plw.dat@" mean_plw.in > mean.plw
bin/track.${EXT} < mean.plw
cat indat/grid1_giss.dat outdat/user_tavg.${EXT} > ${DIRO}/avg_tot_avg_${FIELD}_${LOWES}_plw.dat
\rm outdat/user_tavg.${EXT}
cd ${DIRO}
${TRACK}/converters/bin2nc tot_avg_${FIELD}_${LOWES}_plw.dat tot_avg_${FIELD}_${LOWES}_plw
${TRACK}/converters/bin2nc avg_tot_avg_${FIELD}_${LOWES}_plw.dat avg_tot_avg_${FIELD}_${LOWES}_plw

