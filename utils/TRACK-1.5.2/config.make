#!/bin/sh
#
# Configure shell script for editing the Makefile
#

if [ "$EXEC_EXIST" ]; then exit 0; fi

usage="[MAKEFILE.in][-execn=EXEC] [-fext=FEXT] [-inf=INFILE] [-tim=TIME_AV] [-user=USER]"

configname=$0

MAKEIN=$1

echo "config.make Usage: $usage"

execname=
filextns=
userpath=`(cd ..; pwd)`
track_dir=`pwd | sed -e "s@${userpath}/@@"`
indatfil=
timavfil=TIME_AVG

INCLUDE=./include
SRC=./src
P77=./f77_param

name=`echo $configname | sed 's|/[^/]*$||'`

for arg in $*;
do
  case $arg in
#
# exec name
#
   -execn=* | -exec=* | -exe=* | -ex=* | -e=*)
	execname=`echo $arg | sed 's/[+-]*e[a-z]*=//'`
	cat ${MAKEIN} | sed "s@track@${execname}@" > Makefile 
	;;
#
# output files extension
#
   -fext=* | -fex=* | -fe=* | -f=*)
	filextns=`echo $arg | sed 's/[+-]*f[a-z]*=//'`
        cat $INCLUDE/files_out.in | sed "s@\?@.${filextns}@;s@\!@${filextns}@" > $INCLUDE/files_out.h
	cat $INCLUDE/file_cat_out.in | sed "s@\?@.${filextns}@;s@\!@${filextns}@" > $INCLUDE/file_cat_out.h
	;;
#
# input data file name
#
   -inf=* | -in=* | -i=*)
	indatfil=`echo $arg | sed 's/[+-]*i[a-z]*=//'`
        cat $INCLUDE/files_in.in | sed "s@inputt@${indatfil}@;s@?@.${filextns}@" \
 > $INCLUDE/files_in.temp
	cp $INCLUDE/files_in.temp $INCLUDE/files_in.h
	;;
#
# another input file name
#
   -tim=* | -ti=* | -t=*)
	timavfil=`echo $arg | sed 's/[+-]*t[a-z]*=//'`
	if test -r $INCLUDE/files_in.temp
	  then
	    cat $INCLUDE/files_in.temp | sed "s@DTI.*\".*\"@DTIMAVG PATHI##\"${timavfil}\"@" \
 > $INCLUDE/files_in.h
	else
	    cat $INCLUDE/files_in.in | sed "s@DTI.*\".*\"@DTIMAVG PATHI##\"${timavfil}\"@" \
 > $INCLUDE/files_in.h
	fi
	;;
#
# users path
#
   -user=* | -use=* | -us=* | -u=*)
	userpath=`echo $arg | sed 's/[+-]*u[a-z]*=//'`
        track_dir=`pwd | sed -e "s@${userpath}@@"`

	cat $INCLUDE/file_path.in | sed "s@USER .*@USER       ${userpath}@;s@S_TRACK@${track_dir}@" \
 > $INCLUDE/file_path.h

	CN=0
	outpath=${userpath}${track_dir}/outdat/
	opath=$outpath

	while test -n "$opath"
	do
	  CN=`expr $CN + 1`
	  opath=`echo $opath | sed -e 's/.//'`
	done

	cat $P77/path.ptr.in | sed "s/\*[0-9]*/\*${CN}/;s@=['a-zA-Z0-9/._]*@=\'${outpath}\'@" > $P77/path.ptr


	;;
   *)
        if [ ! "$arg" = "$MAKEIN" ]
        then 
	  echo Incorrect option for $configname
          echo "Usage: $name $usage"
	  exit 1
        fi
	;;

   esac

done


if test -r $INCLUDE/files_in.temp
then
  rm $INCLUDE/files_in.temp
fi

exit 0
