#!/bin/csh 
#$ -S /bin/csh
##$ -o gfs.log
#$ -cwd
#$ -V

source /etc/profile.d/sge-binaries.csh

set DATE=$1
set TIME=$2

@ ST = 0
@ END = 20

qsub rungfs.csh $DATE $TIME 0 0 NH vor
qsub rungfs.csh $DATE $TIME 0 0 NH mslp
qsub rungfs.csh $DATE $TIME 0 0 SH vor
qsub rungfs.csh $DATE $TIME 0 0 SH mslp

while($ST <= $END)

   qsub rungfs.csh $DATE $TIME 1 $ST NH vor
   qsub rungfs.csh $DATE $TIME 1 $ST NH mslp
   qsub rungfs.csh $DATE $TIME 1 $ST SH vor
   qsub rungfs.csh $DATE $TIME 1 $ST SH mslp

   @ ST ++
end

#eps_match.csh $DATE $TIME NH vor >> /dev/null
#eps_match.csh $DATE $TIME NH mslp >> /dev/null
#eps_match.csh $DATE $TIME SH vor >> /dev/null
#eps_match.csh $DATE $TIME SH mslp >> /dev/null


