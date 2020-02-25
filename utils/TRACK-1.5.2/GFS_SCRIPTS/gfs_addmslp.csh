#!/bin/csh

set me = `readlink -f ${0}` ; set GFS_SCRIPTS = `dirname $me`
source ${GFS_SCRIPTS}/common_settings.csh

# Current setup adds MSLP

set DATE=$1      # YYYYMMDD e.g. 20100106
set TIME=$2      # HH, e.g. 06
set TYPE=$3      # deteriministic (0) or ensemble (1)

@ PRT=$4      # perturbation, 0 for control and deterministic

set HEMI=$5
set FIELD=$6     # vor or mslp

set SUBDIR=${DATE}${TIME}

if( "$TYPE" == "1" && ($PRT > 20 || $PRT < 0 )) then
  echo "ERROR: perturbation ID $PERT is not valid."
endif

if($PRT < 10) then
  set PERT=0$PRT
else
  set PERT=$PRT
endif

if( "$TYPE" == "0" ) then
  set HTML="${BASE_D}/gfs_hd${DATE}/gfs_hd_${TIME}z?time,lat,lon,prmslmsl[0:1:64][0:1:360][0:1:719]"
  if( "$FIELD" == "vor" ) then
     set DIROUT=GFS_DET_${DATE}${TIME}_VOR850
     set EXT=gfs_det_${DATE}${TIME}_${HEMI}_vor

  else if( "$FIELD" == "mslp" ) then
     set DIROUT=GFS_DET_${DATE}${TIME}_MSLP
     set EXT=gfs_det_${DATE}${TIME}_${HEMI}_mslp

  endif

else if( "$TYPE" == "1" ) then
  if( "$PERT" == "00" ) then 
    set HTML="${BASE_E}/gens${DATE}/gec00_${TIME}z?ens,time,lat,lon,prmslmsl[0:1:0][0:1:64][0:1:180][0:1:359]"
    if( "$FIELD" == "vor" ) then
       set DIROUT=GFS_EPS_CNTL_${DATE}${TIME}_VOR850
       set EXT=gfs_cntl_${DATE}${TIME}_${HEMI}_vor
    else if( "$FIELD" == "mslp" ) then
       set DIROUT=GFS_EPS_CNTL_${DATE}${TIME}_MSLP
       set EXT=gfs_cntl_${DATE}${TIME}_${HEMI}_mslp
    endif

  else
    set HTML="${BASE_E}/gens${DATE}/gep${PERT}_${TIME}z?ens,time,lat,lon,prmslmsl[0:1:0][0:1:64][0:1:180][0:1:359]"
    if( "$FIELD" == "vor" ) then
       set DIROUT=GFS_EPS_PB${PERT}_${DATE}${TIME}_VOR850
       set EXT=gfs_pb${PERT}_${DATE}${TIME}_${HEMI}_vor
    else if( "$FIELD" == "mslp" ) then
       set DIROUT=GFS_EPS_PB${PERT}_${DATE}${TIME}_MSLP
       set EXT=gfs_pb${PERT}_${DATE}${TIME}_${HEMI}_mslp
    endif

  endif 

else
  echo "ERROR: forecast type $TYPE not valid."
  exit(1)
endif

\rm tmpin/addmslp.$EXT
if("$TYPE" == "0" )then
  set AVFILE=GFS/det_addmslp.in
else
  set AVFILE=GFS/ens_addmslp.in
endif


if( "$FIELD" == "vor" ) then

   if( "$HEMI" == "NH" ) then
      set DIROUT2=${GFSDIR}/${SUBDIR}/NH/${DIROUT}
      sed -e "s/EXT/${EXT}/;s@HTML@${HTML}@;s@GFSDET@${DIROUT2}@" $AVFILE > tmpin/addmslp.$EXT
      set FFF=ff_trs_pos.addvor
   else
      set DIROUT2=${GFSDIR}/${SUBDIR}/SH/${DIROUT}
      sed -e "s/EXT/${EXT}/;s/pos/neg/;s@HTML@${HTML}@;s@GFSDET@${DIROUT2}@" $AVFILE > tmpin/addmslp.$EXT
      set FFF=ff_trs_neg.addvor
   endif

else if( "$FIELD" == "mslp" ) then
   if( "$HEMI" == "NH" ) then
      set DIROUT2=${GFSDIR}/${SUBDIR}/NH/${DIROUT}
      sed -e "s/EXT/${EXT}/;s/pos/neg/;s@HTML@${HTML}@;s@GFSDET@${DIROUT2}@" $AVFILE > tmpin/addmslp.$EXT
      set FFF=ff_trs_neg.addvor
   else
      set DIROUT2=${GFSDIR}/${SUBDIR}/SH/${DIROUT}
      sed -e "s/EXT/${EXT}/;s/pos/neg/;s@HTML@${HTML}@;s@GFSDET@${DIROUT2}@" $AVFILE > tmpin/addmslp.$EXT
      set FFF=ff_trs_neg.addvor
   endif
endif

$track -w "${HTML}" -f $EXT < tmpin/addmslp.$EXT
mv outdat/ff_trs.${EXT}_addfld ${DIROUT2}/${FFF}_addmslp

\rm tmpin/addmslp.$EXT
\rm outdat/*${EXT}*

