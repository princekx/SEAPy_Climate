#!/usr/bin/csh

#source setup.TRACK
setenv LD_LIBRARY_PATH ${HOME}/local/lib/SUNWspro:${LD_LIBRARY_PATH}

@ START = $1
@ END = $2

@ STEP = $3

set YEAR = $4

set FIELD = $5

set LEVEL = $6

set SEASON = $7
set CDIR = $8
set MODEL = $9

set DIR = /data/tethys/kih/AMIP/${MODEL}/${CDIR}/TSERIES/YEARS

mkdir -p $DIR

set TYPE=${FIELD}${LEVEL}_${YEAR}

echo $TYPE

set EXT=time
set INFILE=RUNDATIN.spectra.in

\rm outdat/RUN.${YEAR}_${FIELD}${LEVEL}


sed -e "s/^[0-9]*\?/${START}/;s/^[0-9]*\~/${STEP}/;s/^[0-9]*\!/${END}/;s/NCAR/${MODEL}/;s/PSL/${FIELD}${LEVEL}/;s/7980/${YEAR}/;s/djf/${SEASON}/" < indat/$INFILE > outdat/RUN.${YEAR}_${FIELD}${LEVEL}

bin/track.$EXT < outdat/RUN.${YEAR}_${FIELD}${LEVEL} > /dev/null

mv outdat/user_tavg.$EXT $DIR/user_tavg.${TYPE}
mv outdat/user_tavg.${EXT}_var $DIR/user_var.${TYPE}

mv outdat/user_tavg.${EXT}_varfil $DIR/user_varfil.${TYPE}


echo "$START $STEP $END"

\rm outdat/RUN.${YEAR}_${FIELD}${LEVEL}

