#!/usr/bin/csh


set HEMI=$1
set MODEL=$2
set SEASON=$3

set FIELD=(MSLP O PV T250 T850 V250 V850 VOR250 VOR850 Z250 Z850 theta_pv2)
#set FIELD=(O PV T250 T850 V250 V850 VOR250 VOR850 Z250 Z850 theta_pv2)

foreach var ($FIELD)

   stats_season.csh $var $HEMI $MODEL $SEASON ERAOP
   stats_season.csh $var $HEMI $MODEL $SEASON ERA
   stats_season.csh $var $HEMI $MODEL $SEASON OP

end
