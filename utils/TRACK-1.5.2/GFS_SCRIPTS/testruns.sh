#!/bin/bash

DATE=$1
TIME=$2

FILEIN=log.${DATE}${TIME}

NEWCN=condor.tryagain.${DATE}${TIME}

\rm $NEWCN 

cat GFS/condor.gfs.head > $NEWCN

STATUS=0

while read LINE
do
  grep -q "$LINE" $FILEIN 
  if [ $? -ne 0 ] 
  then

     STATUS=1
     echo "$LINE"
     TYP=`echo $LINE | cut -d " " -f 1`
     PRT=`echo $LINE | cut -d " " -f 2`
     HMI=`echo $LINE | cut -d " " -f 3`
     FLD=`echo $LINE | cut -d " " -f 4`

     echo $TYP $PRT $HMI $FLD

     FFLD=`echo $FLD | tr "[:lower:]" "[:upper:]"` 

     if [ "$TYPE" = "0" ]
     then
        EXT=gfs_det_${DATE}${TIME}_${HMI}_${FLD}
        \rm -r RUN_AT_gfs_det_${DATE}${TIME}_${HMI}_${FLD}
        \rm ../GFS/${HMI}/GFS_DET_${DATE}${TIME}_$FFLD
     elif [ "$TYPE" = "1" -a "$PRT" = "00" ]
     then
        EXT=gfs_cntl_${DATE}${TIME}_${HMI}_${FLD}
        \rm -r RUN_AT_gfs_cntl_${DATE}${TIME}_${HMI}_${FLD}
        \rm ../GFS/${HMI}/GFS_CNTL_${DATE}${TIME}_$FFLD
     else
        EXT=gfs_pb${PRT}_${DATE}${TIME}_${HMI}_${FLD}
        \rm -r RUN_AT_gfs_pb${PRT}_${DATE}${TIME}_${HMI}_${FLD}
        \rm ../GFS/${HMI}/GFS_PB${PRT}_${DATE}${TIME}_$FFLD
     fi

     \rm outdat/*${EXT}*
     echo " " >> $NEWCN
     echo "Arguments       = $DATE $TIME $TYP $PRT $HMI $FLD" >> $NEWCN
     echo 'Output          = CONDOR.OUT/out.$(Process)' >> $NEWCN
     echo 'Queue' >> $NEWCN
  fi
done <  GFS/runs.log.test

exit "$STATUS"
