#!/usr/bin/csh

set FIELD=$1
set HEMI=$2
set MODEL=$3
set SEASON=$4
set MODE=$5

set DIR=/data/medusa/kih/AMIP

bb_stats.job ${MODE}_tot_${FIELD}_pos TOTAL ${DIR}/${MODEL}/${HEMI}_${SEASON}
bb_stats.job ${MODE}_tot_${FIELD}_neg TOTAL ${DIR}/${MODEL}/${HEMI}_${SEASON}
