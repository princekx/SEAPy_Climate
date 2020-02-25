#!/bin/csh -f

if !( -d bin) mkdir bin

#set FILES='meantrd.c read_tracks.c statdmp.c malloc_initl.c read_stats.c realloc_n.c sincos.c new_time.c convert_track.c netcdf_write_stats.c'
set FILES='meantrd.c read_tracks.c statdmp.c malloc_initl.c read_stats.c realloc_n.c new_time.c convert_track.c netcdf_write_stats.c'

echo $FILES

set LIST='COMBINE CONVOLVE COUNT DIFFSTAT ENSEMBLE HART IBT-REF LANDFALL LIFECYCLE LOCATION MISC MULTI-IDENT SAMPLE SPLIT TELE_WEIGHTS TC TILT TRANSADD TRTRANS TR2KML TR2NC TR2PGSL SIEGEL'

echo $LIST

foreach var ($LIST)
   foreach file ($FILES)
      \rm $var/$file
   end
end

# COMBINE


ln -s ../COMMON/meantrd.c COMBINE/meantrd.c
ln -s ../COMMON/read_tracks.c COMBINE/read_tracks.c

# COUNT

ln -s ../COMMON/meantrd.c COUNT/meantrd.c
ln -s ../COMMON/read_tracks.c COUNT/read_tracks.c
ln -s ../COMMON/new_time.c COUNT/new_time.c
#ln -s ../COMMON/sincos.c COUNT/sincos.c


# CONVOLVE

ln -s ../COMMON/meantrd.c CONVOLVE/meantrd.c
ln -s ../COMMON/read_tracks.c CONVOLVE/read_tracks.c

# DIFFSTAT

ln -s ../COMMON/read_stats.c DIFFSTAT/read_stats.c
ln -s ../COMMON/statdmp.c DIFFSTAT/statdmp.c
ln -s ../COMMON/malloc_initl.c DIFFSTAT/malloc_initl.c
ln -s ../COMMON/netcdf_write_stats.c DIFFSTAT/netcdf_write_stats.c

# ENSEMBLE

ln -s ../COMMON/read_tracks.c ENSEMBLE/read_tracks.c
ln -s ../COMMON/meantrd.c ENSEMBLE/meantrd.c
ln -s ../COMMON/realloc_n.c ENSEMBLE/realloc_n.c
ln -s ../COMMON/malloc_initl.c ENSEMBLE/malloc_initl.c
#ln -s ../COMMON/sincos.c ENSEMBLE/sincos.c
ln -s ../COMMON/new_time.c ENSEMBLE/new_time.c
ln -s ../COMMON/convert_track.c ENSEMBLE/convert_track.c

# HART

ln -s ../COMMON/read_tracks.c HART/read_tracks.c
ln -s ../COMMON/meantrd.c HART/meantrd.c
ln -s ../COMMON/realloc_n.c HART/realloc_n.c

# SAMPLE

ln -s ../COMMON/meantrd.c SAMPLE/meantrd.c
ln -s ../COMMON/read_tracks.c SAMPLE/read_tracks.c
ln -s ../COMMON/read_stats.c SAMPLE/read_stats.c
ln -s ../COMMON/malloc_initl.c SAMPLE/malloc_initl.c
ln -s ../COMMON/statdmp.c SAMPLE/statdmp.c
ln -s ../COMMON/realloc_n.c SAMPLE/realloc_n.c
ln -s ../COMMON/netcdf_write_stats.c SAMPLE/netcdf_write_stats.c

# SPLIT

ln -s ../COMMON/meantrd.c SPLIT/meantrd.c
ln -s ../COMMON/read_tracks.c SPLIT/read_tracks.c

# TELE_WEIGHTS

ln -s ../COMMON/read_tracks.c TELE_WEIGHTS/read_tracks.c
#ln -s ../COMMON/sincos.c TELE_WEIGHTS/sincos.c

# TILT

ln -s ../COMMON/read_tracks.c TILT/read_tracks.c
#ln -s ../COMMON/sincos.c TILT/sincos.c

# TC

ln -s ../COMMON/meantrd.c TC/meantrd.c
ln -s ../COMMON/read_tracks.c TC/read_tracks.c
#ln -s ../COMMON/sincos.c TC/sincos.c
ln -s ../COMMON/convert_track.c TC/convert_track.c
ln -s ../COMMON/realloc_n.c TC/realloc_n.c

# LIFECYCLE

ln -s ../COMMON/read_tracks.c LIFECYCLE/read_tracks.c
ln -s ../COMMON/meantrd.c LIFECYCLE/meantrd.c
ln -s ../COMMON/realloc_n.c LIFECYCLE/realloc_n.c
#ln -s ../COMMON/sincos.c LIFECYCLE/sincos.c

#LOCATION

ln -s ../COMMON/read_tracks.c LOCATION/read_tracks.c
ln -s ../COMMON/meantrd.c LOCATION/meantrd.c

# MISC

ln -s ../COMMON/read_tracks.c MISC/read_tracks.c
ln -s ../COMMON/meantrd.c MISC/meantrd.c
ln -s ../COMMON/realloc_n.c MISC/realloc_n.c

#TRANSADD

ln -s ../COMMON/read_tracks.c TRANSADD/read_tracks.c
ln -s ../COMMON/meantrd.c TRANSADD/meantrd.c

#TRTRANS

ln -s ../COMMON/read_tracks.c TRTRANS/read_tracks.c
ln -s ../COMMON/meantrd.c TRTRANS/meantrd.c

#TR2KML

ln -s ../COMMON/read_tracks.c TR2KML/read_tracks.c
ln -s ../COMMON/meantrd.c TR2KML/meantrd.c
ln -s ../COMMON/new_time.c TR2KML/new_time.c
ln -s ../COMMON/realloc_n.c TR2KML/realloc_n.c
ln -s ../COMMON/malloc_initl.c TR2KML/malloc_initl.c

#TR2PGSL

ln -s ../COMMON/read_tracks.c TR2PGSL/read_tracks.c
ln -s ../COMMON/meantrd.c TR2PGSL/meantrd.c
ln -s ../COMMON/new_time.c TR2PGSL/new_time.c

#SIEGEL

ln -s ../COMMON/new_time.c SIEGEL/new_time.c


#LANDFALL

ln -s ../COMMON/read_tracks.c LANDFALL/read_tracks.c
ln -s ../COMMON/new_time.c LANDFALL/new_time.c
ln -s ../COMMON/convert_track_add.c LANDFALL/convert_track_add.c
ln -s ../COMMON/realloc_n.c LANDFALL/realloc_n.c
#ln -s ../COMMON/sincos.c LANDFALL/sincos.c

#MULTI-IDENT

ln -s ../COMMON/read_tracks.c MULTI-IDENT/read_tracks.c
ln -s ../COMMON/meantrd.c MULTI-IDENT/meantrd.c

#IBT-REF

ln -s ../COMMON/read_tracks.c IBT-REF/read_tracks.c
ln -s ../COMMON/meantrd.c IBT-REF/meantrd.c
ln -s ../COMMON/realloc_n.c IBT-REF/realloc_n.c

#TR2NC
ln -s ../COMMON/read_tracks.c TR2NC/read_tracks.c
ln -s ../COMMON/new_time.c TR2NC/new_time.c
