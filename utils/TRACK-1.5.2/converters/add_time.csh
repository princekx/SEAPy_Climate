#!/usr/bin/csh

set FIELD=$1
set SEASON=$2
set MODEL=$3

set LIST=`ls ${MODEL}_${FIELD}_${SEASON}*.utf`

\rm ttt

foreach var ($LIST)
   sed -e 's/                                                  \$/RUN    1.000     DAYS       0.00 TO     91.75     \$/' $var > ttt
   mv ttt $var
end
