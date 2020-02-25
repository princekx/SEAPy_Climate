#!/usr/bin/csh

set MODEL=$1
set TELE=$2
set MODE=$3
set HEMI=$4
set SEASON=$5

#set FIELDS=(MSLP O PV T250 T850 V250 V850 VOR250 VOR850 Z250 Z850 theta_pv2)

set FIELDS=(VOR850 VOR700 VOR600)

set TELED=${TELE}.TANH
set TYPE=`echo $TELE | tr "[:upper:]" "[:lower:]"`

set BASE=/data/atlas/kih/AMIP/${MODEL}/${HEMI}_${SEASON}

foreach FIELD ($FIELDS)

  set POS=${MODE}_${TYPE}_pos_${FIELD}_pos

  bb_stats.job ${POS} $TELED $BASE

  set POS=${MODE}_${TYPE}_pos_${FIELD}_neg


  bb_stats.job ${POS} $TELED $BASE

  set NEG=${MODE}_${TYPE}_neg_${FIELD}_pos

  bb_stats.job ${NEG} $TELED $BASE

  set NEG=${MODE}_${TYPE}_neg_${FIELD}_neg


  bb_stats.job ${NEG} $TELED $BASE


end
