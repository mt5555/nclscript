#!/bin/tcsh  
if ($#argv < 1) then
  echo "climoplot.sh casename"
  exit
endif

conda activate nglmpl2


if ( $1 =~ sd.nc ) then 
   #
   # standalone HOMME Held Suarez test case data
   #
   set name = `basename $1 .nc`
   set wdir = `pwd`
   set dirname = `basename $wdir`
   ls $1
   if ( $status == 0 ) then
   foreach opt ( andes himalaya latlon )
      contour.py -i $name.nc -c -.15,.15 -r 300x600 -m $opt omega500
      mv -f $name.omega500.pdf $dirname.omega500.$opt.sd.pdf
   end
   endif
   # snapshots
   ls movies/held_suarez0*3.nc
   if ( $status == 0 ) then
      set filename = `\ls movies/held_suarez0*3.nc  | tail -1`
      echo snapshot plot: $filename
      #contour.py -i $filename -t -1  -c -1.0,1.0 -r 300x600 -m andes -p 500  omega
      #mv -f movies/held_suarez0*omega.pdf omega500-t1.pdf
   endif
endif



if ( $1 =~ hs* ) then 
   #
   # standalone HOMME Held Suarez test case data
   #
   set name = `basename $1 .nc`
   set wdir = `pwd`
   set dirname = `basename $wdir`
   ls $1
   if ( $status == 0 ) then
   foreach opt ( andes himalaya latlon )
      contour.py -i $name.nc -c -.15,.15 -r 300x600 -m $opt omega500
      mv -f $name.omega500.pdf $dirname.omega500.$opt.pdf
   end
   endif
   # snapshots
   ls movies/held_suarez0*3.nc
   if ( $status == 0 ) then
      set filename = `\ls movies/held_suarez0*3.nc  | tail -1`
      echo snapshot plot: $filename
      #contour.py -i $filename -t -1  -c -1.0,1.0 -r 300x600 -m andes -p 500  omega
      #mv -f movies/held_suarez0*omega.pdf omega500-t1.pdf
   endif
else
   #
   # EAM data
   #
   set base = $1
   set Y1 = 20    # 2
   set Y2 = 24    # 6
   set Y1M = ${Y1}01     # for climos with out the jfd option, these might run
   set Y2M = ${Y2}12     # from ${Y1-1}12 -> ${Y2}11
   
   \ls climo/${base}_ANN*_0*${Y1M}_0*${Y2M}_climo.nc
   if ( $status ) then
      echo no climo files, running nco...

      \ls ${base}.eam.h0.*.nc
      if ( $status == 0 ) then
         set filename = `\ls ${base}.eam.h0.*.nc | tail -1`
         echo running ncclimo:  $filename
         # consider -a jfd?
         ncclimo -a jfd -m eam -s $Y1 -e $Y2 -c $base -i . -o climo -O regrid \
            -r ~ac.zender/data/maps/map_ne30pg2_to_cmip6_180x360_nco.20200901.nc 
      else
         echo cant find h0 files
         exit -1
      endif
   else
      echo found climo files, running plots
   endif
    
   cd climo   # plot native grid data below
   \ls ${base}_ANN*_0*${Y1M}_0*${Y2M}_climo.nc 
   if ( $status ) then
      echo cant find climo files
      exit -1
   endif
   
   set filename = `\ls ${base}_ANN*_0*${Y1M}_0*${Y2M}_climo.nc  | tail -1`
   echo running OMEGA contour plots:  $filename
   set name = `basename $filename .nc`
   echo $name

   foreach opt ( andes himalaya latlon )
      echo $opt
      contour.py -i $name.nc -c -.15,.15 -r 300x600 -m $opt OMEGA500
      contour.py -i $name.nc -c -.15,.15 -r 300x600 -m $opt -p 500 DYN_OMEGA
      if ( $Y1 == 2 && $Y2 == 6 ) then
         mv $name.DYN_OMEGA.pdf $base.DYN_OMEGA.$opt.pdf
         mv $name.OMEGA500.pdf $base.OMEGA500.$opt.pdf
      else
         mv $name.DYN_OMEGA.pdf $base.$Y1-$Y2.DYN_OMEGA.$opt.pdf
         mv $name.OMEGA500.pdf $base.$Y1-$Y2.OMEGA500.$opt.pdf
      endif
   end
endif
