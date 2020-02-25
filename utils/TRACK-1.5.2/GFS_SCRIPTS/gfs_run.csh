#!/bin/csh -x
#$ -S /bin/csh
##$ -o gfs.log
#$ -cwd
#$ -V

set me = `readlink -f ${0}` ; set GFS_SCRIPTS = `dirname $me`
source ${GFS_SCRIPTS}/common_settings.csh

set DATE=$1      # YYYYMMDD e.g. 20100106
set TIME=$2      # HH, e.g. 06
set TYPE=$3      # deteriministic (0) or ensemble (1)

@ PRT=$4         # perturbation, 0 for control and deterministic

set HEMI=$5
set FIELD=$6     # vor or mslp

set SUBDIR=${DATE}${TIME}

# check /tmp/vxx05162 exists

#if( ! -d /tmp/vxx05162) mkdir /tmp/vxx05162

# check tmpin exists

if( ! -d tmpin ) mkdir tmpin

if( "$TYPE" == "1" && ($PRT > 20 || $PRT < 0 )) then
  echo "ERROR: perturbation ID $PERT is not valid."
endif

if($PRT < 10) then
  set PERT=0$PRT
else
  set PERT=$PRT
endif

if( "$TYPE" == "0" )then

  if( "$FIELD" == "vor" ) then
     set HTML="${BASE_D}/gfs${DATE}/gfs_0p50_${TIME}z?time,lev[6:1:6],lat,lon,absvprs[0:1:80][6:1:6][0:1:360][0:1:719]"
     set DIROUT=GFS_DET_${DATE}${TIME}_VOR850
     set EXT=gfs_det_${DATE}${TIME}_${HEMI}_vor

  else if( "$FIELD" == "mslp" ) then
     set HTML="${BASE_D}/gfs${DATE}/gfs_0p50_${TIME}z?time,lat,lon,prmslmsl[0:1:80][0:1:360][0:1:719]"
     set DIROUT=GFS_DET_${DATE}${TIME}_MSLP
     set EXT=gfs_det_${DATE}${TIME}_${HEMI}_mslp

  endif

else if( "$TYPE" == "1" )then

  if( "$PERT" == "00" )then 
    if( "$FIELD" == "vor" ) then
       set HTML="${BASE_E}/gens${DATE}/gec00_${TIME}z?ens,time,lev[5:1:5],lat,lon,absvprs[0:1:0][0:1:64][5:1:5][0:1:180][0:1:359]"
       set DIROUT=GFS_EPS_CNTL_${DATE}${TIME}_VOR850
       set EXT=gfs_cntl_${DATE}${TIME}_${HEMI}_vor

    else if( "$FIELD" == "mslp" ) then
       set HTML="${BASE_E}/gens${DATE}/gec00_${TIME}z?ens,time,lat,lon,prmslmsl[0:1:0][0:1:64][0:1:180][0:1:359]"
       set DIROUT=GFS_EPS_CNTL_${DATE}${TIME}_MSLP
       set EXT=gfs_cntl_${DATE}${TIME}_${HEMI}_mslp

    endif
  else
    if( "$FIELD" == "vor" ) then
       set HTML="${BASE_E}/gens${DATE}/gep${PERT}_${TIME}z?ens,time,lev[5:1:5],lat,lon,absvprs[0:1:0][0:1:64][5:1:5][0:1:180][0:1:359]"
       set DIROUT=GFS_EPS_PB${PERT}_${DATE}${TIME}_VOR850
       set EXT=gfs_pb${PERT}_${DATE}${TIME}_${HEMI}_vor

    else if( "$FIELD" == "mslp" ) then
       set HTML="${BASE_E}/gens${DATE}/gep${PERT}_${TIME}z?ens,time,lat,lon,prmslmsl[0:1:0][0:1:64][0:1:180][0:1:359]"
       set DIROUT=GFS_EPS_PB${PERT}_${DATE}${TIME}_MSLP
       set EXT=gfs_pb${PERT}_${DATE}${TIME}_${HEMI}_mslp

    endif 
  endif

else
  echo "ERROR: forecast type $TYPE not valid."
  exit(1)
endif


if( "$FIELD" == "vor" ) then
   if( "$HEMI" == "NH" ) then
     set INI=initial_NH_T42.dat
     set FFF=ff_trs_pos.addvor_addmslp_add925winds_add10mwinds_addprecip
   else if( "$HEMI" == "SH" ) then
     set INI=initial_SH_T42.dat
     set FFF=ff_trs_neg.addvor_addmslp_add925winds_add10mwinds_addprecip
   endif

   if( "$TYPE" == "0" ) then
      $track -f $EXT -w "$HTML" < GFS/spec_filt_T42_det.in 
   else
      $track -f $EXT -w "$HTML" < GFS/spec_filt_T42_eps.in 
   endif

   \rm indat/specfil.${EXT}_band000.dat
   ln -s `pwd`/outdat//specfil.${EXT}_band000 indat/specfil.${EXT}_band000.dat
   $master -c=$DIROUT -d=now -e=`basename $track` -f=${EXT} -i=specfil.${EXT}_band000.dat -j=RUN_AT_GFS.in -k=$INI -n=1,100,1 -o=${GFSDIR}/${SUBDIR}/${HEMI} -r=RUN_AT_ -s=RUNDATIN.GFS_VOR

else if( "$FIELD" == "mslp" ) then
   set FFF=ff_trs_neg.addvor_addmslp_add925winds_add10mwinds_addprecip
   if( "$HEMI" == "NH" ) then
     set INI=initial_NH_T63.dat
   else if( "$HEMI" == "SH" ) then
     set INI=initial_SH_T63.dat
   endif

   if( "$TYPE" == "0" ) then
      $track -f $EXT -w "$HTML" < GFS/spec_filt_T63_det.in 
   else
      $track -f $EXT -w "$HTML" < GFS/spec_filt_T63_eps.in 
   endif

   \rm indat/specfil.${EXT}_band000.dat
   ln -s `pwd`/outdat/specfil.${EXT}_band000 indat/specfil.${EXT}_band000.dat
   $master -c=$DIROUT -d=now -e=`basename $track` -f=${EXT} -i=specfil.${EXT}_band000.dat -j=RUN_AT_GFS.in -k=$INI -n=1,100,1 -o=${GFSDIR}/${SUBDIR}/${HEMI} -r=RUN_AT_ -s=RUNDATIN.GFS_MSLP
endif

\rm outdat/*${EXT}* indat/specfil.${EXT}_band000.dat
gfs_addvor.csh $DATE $TIME $TYPE $PRT $HEMI $FIELD 
gfs_addmslp.csh $DATE $TIME $TYPE $PRT $HEMI $FIELD 
gfs_add925winds.csh $DATE $TIME $TYPE $PRT $HEMI $FIELD
gfs_add10mwinds.csh $DATE $TIME $TYPE $PRT $HEMI $FIELD 
gfs_addprecip.csh $DATE $TIME $TYPE $PRT $HEMI $FIELD 

set DIROUT2=${GFSDIR}/${SUBDIR}/${HEMI}/$DIROUT

utils/bin/count ${DIROUT2}/${FFF} 0 0 5 4 0 ${DATE}${TIME} 6
mv ${DIROUT2}/${FFF}.new ${DIROUT2}/${FFF}.time

if( -f ${DIROUT2}/${FFF}.time ) then
   echo "$TYPE $PERT $HEMI $FIELD" >> log.${DATE}${TIME}
endif
