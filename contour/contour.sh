#!/bin/tcsh -f
#
# example:  contour.sh    VOR250  files
#
set field = $1
shift

set testfile = $1
set ncols = `ncdump -h $testfile  | grep "ncol = " `
if ( $status ) then
   echo "looks like lat/lon grid: running contour_latlon.ncl"
   set cmd = ~/ncl/contour_latlon.ncl
else
   echo "looks like native grid: running contour_native.ncl"
   set cmd = ~/ncl/contour_native.ncl
endif

set files = "$*"

set arg0 = "projection=0"
#set arg0 = "projection=1"    #orthographic over N pole
#set arg0 = "projection=2"    # north america
#set arg0 = "projection=3"    # indonesia
#set arg0 = "projection=4"    # Andes
#set arg0 = "projection=5"    # Himalyas
set arg1 = field=\"$field\"
set arg2 = `echo files=\"$files\"`


ncl  "$arg0" "$arg1"  "$arg2" $cmd

