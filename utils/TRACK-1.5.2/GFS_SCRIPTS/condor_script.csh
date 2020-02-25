#!/bin/csh 

set DATE=$1
set TIME=$2

set REG=$3   # region switch, '0' is extra-tropics, '1' is tropics

\rm data/zone.dat0
if($REG == 0)then
   set SUBMIT=condor.gfs_${DATE}${TIME}_ET
   ln -s data/zone_st.dat data/zone.dat0
else if($REG == 1)then
   set SUBMIT=condor.gfs_${DATE}${TIME}_TR
   ln -s data/zone_tr.dat data/zone.dat0 
else
   echo "Region switch incorrect, use '0' or '1'"
   exit(1)
endif

@ ST = 0
@ END = 20

@ NTRY = 5

# create log file

\rm log.${DATE}${TIME}
touch log.${DATE}${TIME}

\rm $SUBMIT 

cat GFS/condor.gfs.head > $SUBMIT

echo "Arguments       = $DATE $TIME 0 0 NH vor" >> $SUBMIT
echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
echo 'Queue' >> $SUBMIT
echo ' ' >> $SUBMIT
echo "Arguments       = $DATE $TIME 0 0 NH mslp" >> $SUBMIT
echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
echo 'Queue' >> $SUBMIT
echo ' ' >> $SUBMIT
echo "Arguments       = $DATE $TIME 0 0 SH vor" >> $SUBMIT
echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
echo 'Queue' >> $SUBMIT
echo ' ' >> $SUBMIT
echo "Arguments       = $DATE $TIME 0 0 SH mslp" >> $SUBMIT
echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
echo 'Queue' >> $SUBMIT
echo ' ' >> $SUBMIT

while($ST <= $END)

   echo "Arguments       = $DATE $TIME 1 $ST NH vor" >> $SUBMIT
   echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
   echo 'Queue' >> $SUBMIT
   echo ' ' >> $SUBMIT

   echo "Arguments       = $DATE $TIME 1 $ST NH mslp" >> $SUBMIT
   echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
   echo 'Queue' >> $SUBMIT
   echo ' ' >> $SUBMIT

   echo "Arguments       = $DATE $TIME 1 $ST SH vor" >> $SUBMIT
   echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
   echo 'Queue' >> $SUBMIT
   echo ' ' >> $SUBMIT

   echo "Arguments       = $DATE $TIME 1 $ST SH mslp" >> $SUBMIT
   echo 'Output          = CONDOR.OUT/out.$(Process)' >> $SUBMIT
   echo 'Queue' >> $SUBMIT

   @ ST ++
end


set CLUSTER=`condor_submit $SUBMIT | grep "submitted to cluster" |  sed -e "s/^.*cluster//;s/\.//"`
if( $? == 1 )then
   echo "**** condor_submit error ****"
   exit(1)
endif
echo $CLUSTER

@ RET = 0
while( ! $RET )
   condor_q | grep -q $CLUSTER
   @ RET = $?
   sleep 10m
end


# Retry's

if($REG == 0)then
   set NEWCN=condor.tryagain.${DATE}${TIME}_ET
else if($REG == 1)then
   set NEWCN=condor.tryagain.${DATE}${TIME}_TR
else
   echo "Region switch incorrect, use '0' or '1'"
   exit(1)
endif

@ NT = 0
while( $NT < $NTRY )
   testruns.sh $DATE $TIME
   if( $? == 0 )then
      echo "**** no further re-try's required ****"
      break
   endif

   set CLUSTER=`condor_submit $NEWCN | grep "submitted to cluster" |  sed -e "s/^.*cluster//;s/\.//"`
   echo $CLUSTER
   @ RET = 0
   while( ! $RET )
      condor_q | grep -q $CLUSTER
      @ RET = $?
      sleep 10m
   end


   @ NT ++
end

eps_match.csh $DATE $TIME NH vor >> /dev/null
eps_match.csh $DATE $TIME NH mslp >> /dev/null
eps_match.csh $DATE $TIME SH vor >> /dev/null
eps_match.csh $DATE $TIME SH mslp >> /dev/null


