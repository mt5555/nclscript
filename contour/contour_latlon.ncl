load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;***********************************************************************;
; detect type of lat/lon grid
;
; gauss_grid =-1    routine failed 
; gauss_grid = 1    Gauss grid (use shaeg routines)
; gauss_grid = 2    Equal angle grid with poles  (use shaec routines)
; gauss_grid = 3    Equal angle grid without poles  (must interpolate first)
;
;***********************************************************************;
function get_gridtype(lat)
local d1,d2
begin
  gauss_grid = -1
  if (lat(0) .lt. -89.99999 ) then
    print ("detected equal angle grid with poles")
    gauss_grid = 2    ; use shaec routines
  else
    d1=lat(1)-lat(0)
    d2=lat(2)-lat(1)
    print(d1-d2)
    if (abs(d1-d2) .lt. 1e-8 ) then
      print ("detected equal angle grid without poles.  lat(0)="+lat(0))
      gauss_grid = 3    ; need to interpolate using fo2fsh
    else
      print("detected Gauss grid")
      gauss_grid = 1  
    end if
  end if
  return(gauss_grid)
end

begin
setvalues NhlGetWorkspaceObjectId() 
   "wsMaximumSize" : 50000000
end setvalues 


;***********************************************************************;
; 
; run via the tcsh script in ~/bin:
;    contour.sh  PRECT *.nc
;
; or:
;     ncl 'field="PRECT"' 'files="*.nc"'  contour_latlon.ncl
; or
;     ncl contour_latlon.ncl
;     and edit field and files variables below
; 
;
;***********************************************************************;
;  output_format = "ncgm"
  output_format = "pdf"
;  output_format = "X11"

if (isvar("projection")) then         ; pressure level on command line?
else
  projection = 0  ; standard lat/lon
  ;projection = 1  ; Orthographic over N pole
  ;projection = 2  ; zoom over NH refinement patch- used for APE vorticity plots
end if

if (isvar("field")) then         ; pressure level on command line?
else
  field = "PS"
end if

if (isvar("files")) then         ; pressure level on command line?
else
  files = "*.nc"
end if
print("field="+field+" file pattern="+files)
fnames =  systemfunc("ls "+files);


  f      = addfiles(fnames,"r")       ; Open netCDF files.	
  print("f="+fnames)
  ; change the name of the field to plot here
  ; current field names are Q, ps, ke, u, v, T, zeta - names are case sensitive
  ; prec, accum


  pnew = 0d0  ;  set to > 0 to interpolate to pressure level
  pscale = 0d0;  set to > 0 to scale input before plotting
  scale = 1d0 ;  var scaling
  res                     = True



  plot_spectra=False

; default is to plot level 18. especially for OMEGA
; override below if desired
  firstlev=18
  lastlev=18
  levinc=1

  ; timelevels to plot:
  firsttime=0
  lasttime_max=999999
  lasttime_max=0
;  firsttime=lasttime
  timeinc=1

  longname = ""
  units = ""
  sumprec = False
  compute_vor = False




 if ( field .ne. "OMEGA500" ) then
     firsttime=0
     lasttime_max=18

     firsttime=0
     lasttime_max=36
     timeinc=1
     
 ;    firsttime=9
;     lasttime_max=9

     ;hdeg=10.
     hdeg=3.
     ;res@mpMinLatF = 90-hdeg
     ;res@mpMaxLatF = 90.
     ;res@mpMinLonF = 295.-hdeg
     ;res@mpMaxLonF = 295. + hdeg
;     res@mpProjection      = "Orthographic"
;       res@mpCenterLatF      =  80.
;       res@mpCenterLonF         = 296.

;     res@mpLimitMode="Angles"
;       res@mpLeftAngleF=hdeg
;      res@mpRightAngleF=hdeg
;      res@mpBottomAngleF=hdeg
;      res@mpTopAngleF=hdeg

       
  end if













  res@cnLevelSelectionMode = "AutomaticLevels"

  if (field .eq. "TBOT" ) then
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF =  160.
    res@cnMaxLevelValF = 300.
    res@cnLevelSpacingF =  5.
  end if
  if (field .eq. "PRECT" ) then
    ; convert from m/s to mm/day
    scale=1000d0*(24*3600)
    units="mm/day"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF =  0.
    res@cnMaxLevelValF = 10.
    res@cnMaxLevelValF = 10.
    res@cnLevelSpacingF =  0.25
  end if
  if (field .eq. "PRECC" ) then
    scale=1000d0*(24*3600)
    units="mm/day"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF =  0.
    res@cnMaxLevelValF =  18 / 1000.  / (24*3600.)
    res@cnLevelSpacingF =  res@cnMaxLevelValF/25.
  end if
  if (field .eq. "PSL" ) then
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF =  96000.
    res@cnMaxLevelValF = 105000.
    res@cnLevelSpacingF =  200.
  end if
  if (field .eq. "PS" ) then
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF =  52000.
    res@cnMaxLevelValF = 105000.
    res@cnLevelSpacingF =  2000.
    res@cnMinLevelValF =  99000.
    res@cnMaxLevelValF = 103000.
    res@cnLevelSpacingF =   200.
  end if
  if (field .eq. "div" ) 
    firstlev=23  ; 18 is the default, 23 common choice
    lastlev=firstlev
    levinc=1

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -4.e-5
    res@cnMaxLevelValF = 4.e-5
    res@cnLevelSpacingF = 4.e-5/10

    ;res@cnMinLevelValF = -2.0e-6 ; means
    ;res@cnMaxLevelValF = 2.0e-6
    ;res@cnLevelSpacingF = 2.0e-6/10
    if (firstlev .eq. 18) 
      res@cnMinLevelValF = -2.0e-5
      res@cnMaxLevelValF = 2.0e-5
      res@cnLevelSpacingF = .2e-5
    end if

  end if
  if (field .eq. "OMEGA" ) 
      pnew=850
;; monthly means: Pa/s
    res@cnLevelSelectionMode = "ManualLevels"
     res@cnMinLevelValF = -.2
     res@cnMaxLevelValF = .2
     res@cnLevelSpacingF = .01

     res@cnMinLevelValF = -2.0
     res@cnMaxLevelValF = 2.0
     res@cnLevelSpacingF = .2


  end if
  if (field .eq. "OMEGA500" ) 
    field = "OMEGA500"
    res@cnLevelSelectionMode = "ManualLevels"
;; monthly means: Pa/s
     res@cnMinLevelValF = -.2
     res@cnMaxLevelValF = .2
     res@cnLevelSpacingF = .01
;; monthly means 850mb
;    pnew=850.
;    res@cnMinLevelValF = -.2
;    res@cnMaxLevelValF =  .2
;    res@cnLevelSpacingF = .01
;; snapshots
;    res@cnMinLevelValF = -1.0
;    res@cnMaxLevelValF = 1.0
;    res@cnLevelSpacingF = .05
   end if

  if (field .eq. "PHIS_SM")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -500
    res@cnMaxLevelValF = 55000
    res@cnLevelSpacingF = 1000
    plot_spectra=True
    kmax= 512  ; max wave number for spectra plot
  end if

  if (field .eq. "TMQ")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.
    res@cnMaxLevelValF = 65.
    res@cnMaxLevelValF = 50.
    res@cnLevelSpacingF = .5
  end if

  if (field .eq. "VOR250")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.0001
    res@cnMaxLevelValF =  0.0001
    res@cnLevelSpacingF = 0.0001/20

    ; VOR250 and DIV250 computed via my interpvordiv.job script
    ; will have the long name and units of U250. fix that:
    longname = "VOR250"
    units = "1/s"
  end if
  if (field .eq. "VOR850")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.0004
    res@cnMaxLevelValF =  0.0004
    res@cnLevelSpacingF = 0.0004/20

    ; VOR250 and DIV250 computed via my interpvordiv.job script
    ; will have the long name and units of U250. fix that:
    longname = "VOR850"
    units = "1/s"
  end if

  if (field .eq. "U250")
; this is for u250,v250 data that was overwritten by vor/div using HOMME's 
; interpolation utility vor/div option
; BECAUASE: ncrename command corrupts data so we cant rename these fields
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.0001
    res@cnMaxLevelValF =  0.0001
    res@cnLevelSpacingF = 0.0001/20
    longname = "VOR250"
    units = "1/s"
  end if
  if (field .eq. "U850")
; this is for u250,v250 data that was overwritten by vor/div using HOMME's 
; interpolation utility vor/div option
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.0004
    res@cnMaxLevelValF =  0.0004
    res@cnLevelSpacingF = 0.0004/20
    longname = "VOR850"
    units = "1/s"
  end if


  if (field .eq. "SPEEDBOT")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 40.
    res@cnLevelSpacingF = 40/40.

;    res@cnMinLevelValF = -.005
;    res@cnMaxLevelValF = .005
;    res@cnLevelSpacingF = .005/20.

    longname = "Wind Speed"
    units = "m/s"
  end if

  if (field .eq. "mask")
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = .05
;    longname = "Wind Speed"
;    units = "m/s"
  end if

  if (field .eq. "zeta" ) 
    ; pnew=500
    ;longname = "VOR500"
    firstlev=23  ; 18 is the default, 23 common choice
    lastlev=firstlev
    levinc=1

;; monthly means: Pa/s
    res@cnLevelSelectionMode = "ManualLevels"
    zrange=.0001  ; instantaneous
;    zrange=.00001  ; means (level 23), no topo
    res@cnMinLevelValF = -zrange
    res@cnMaxLevelValF =  zrange
    res@cnLevelSpacingF = zrange/20

  end if





wks = gsn_open_wks(output_format,field)
print ("Writing file "+field+"."+output_format)
fcnt = dimsizes(fnames)
nlat = dimsizes(f[0]->lat)
nlon = dimsizes(f[0]->lon)

nlat_max = nlat
do fn=1,fcnt-1
   if (dimsizes(f[fn]->lat) .gt. nlat_max ) then
     nlat_max = dimsizes(f[fn]->lat)
   end if
end do

spc_all = new ( (/fcnt,nlat_max/), "double", 1d-99 )
spc_all = spc_all@_FillValue
spc_all!0 = "wavenumber"
x = new ( (/nlat_max/), "double", 1d-99 )
x = x@_FillValue

do fn=0,fcnt-1
  f1=f[fn]
  
  nlat = dimsizes(f1->lat)
  nlon = dimsizes(f1->lon)
  
  use_filename=True  
  if (isatt(f1,"case")) then
    if (f1@case .eq. "camrun") then
    else
      f1name = f1@case
      use_filename=False
    end if
  end if
  if (use_filename) then
                                ; remove up to last /:
    tmp1 = stringtochar(fnames(fn))
    l1 = 0
    l2 = dimsizes(tmp1)
    l2 = l2-1
    do while isStrSubset(chartostring(tmp1(l1:l2)),"/")
                                ;print("range: "+l1+" "+l2) 
      ii = indStrSubset(chartostring(tmp1(l1:l2)),"/")
                                ;print("ii = "+ii)
      l1 = l1+ii+1
                                ; print("substr "+tmp1(l1:l2))
    end do
                                ; remove .cam2* 
    if isStrSubset(chartostring(tmp1(l1:l2)),".cam2")  then
      iv = indStrSubset(chartostring(tmp1(l1:l2)),".cam2")
      ii = iv(0)
      l2 = l1+ii-1
    end if
    f1name = chartostring(tmp1(l1:l2))
    delete(tmp1)
  end if
  

  
  if (isfilevar(f1,"time")) then
    ntimes = dimsizes(f1->time)
    time = f1->time
  else
    ntimes = 1
    time = (/0/)
  end if
  if (isfilevar(f1,"lev")) then
    nlev = dimsizes(f1->lev)
  else
    nlev =1
  end if
  if (isfilevar(f1,"lev1")) then  ; ECMWF or NCEP observations from AMWG diagnostics
    pscale = 1./(3600.*24./100.)   ; data in mb/day on pressure levels, scale to Pa/s
    units = "Pa/s"
  end if

  if (ntimes-1 .le. lasttime_max) then
    lasttime=ntimes-1
  else
    lasttime=lasttime_max
  end if


  
;
; note: spherepack:  gauss grid or equally-spaced grid.
;  both grids do not include periodic point in longitidue (lon(0) <> lon(nlon-1))
;  equally spaced grid includes points at both poles
gauss_grid=get_gridtype(f1->lat)
  
  

  field_read=field

  compute_speed=False
  if (field .eq. "SPEEDBOT" ) then
    compute_speed=True
    field_read="UBOT"
    field_read2="VBOT"
  end if

  if (field .eq. "VOR250" ) then
    if (isfilevar(f1,"VOR250")) then
       compute_vor = False
    else
       compute_vor = True
       field_read="U250"
       field_read2="V250"
     end if
  end if

  if (field .eq. "PRECT" ) then
    if (isfilevar(f1,"PRECT")) then
    else
      ; PRECT = PRECL + PRECC
      ; read PRECC, then read PRECL and sum them below:
      sumprec=True        
      field_read = "PRECC"
      longname = "PRECT" 
    end if
  end if

  ndim = dimsizes( getfilevardimsizes(f1,field_read))
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon+"  ndim="+ndim)


  if (ndim .eq. 1) then
    print("error computing ndim.  variable not in file?")
    print("file = "+fnames(fn))
    exit
  end if

  dimnames = getfilevardimnames(f1,field_read)
  if (dimnames(0) .ne. "time")
    firstdim_is_time=False
  else
    firstdim_is_time=True
  end if

; ndim=2  assume:   lat,lon  
; ndim=3  assume:   time,lat,lon   OR   levels,lat,lon  
; ndim=4  assume:   time, levels, lat,lon
  if (ndim .le.  3 .and. firstdim_is_time) then
    firstlev=0
    lastlev=0
  end if
  
  
  print("time index: "+firsttime+":"+timeinc+":"+lasttime)
  print("lev index: "+firstlev+":"+levinc+":"+lastlev)
  do n=firsttime,lasttime,timeinc

     var = new((/nlat,nlon/),"double")
     p = new((/nlat,nlon/),"double")
;  var!0="time"
;  var!0="lev"
;  var!0="lat"
;  var!1="lon"
;  var&lat = lat1d  
;  var&lon = lon1d
;  var&lat@units="degrees_north"
;  var&lon@units="degrees_east"



    if (isfilevar(f1,"hyam")) then
    if ( ndim .eq. 4 ) then
      hyam=(f1->hyam(:))
      ps0 = 1000.0d0
      if (isfilevar(f1,"P0")) then
        ps0 = (f1->P0)
      end if
      if (isfilevar(f1,"PS")) then
        PS=f1->PS(n,:,:)
      end if
      if (isfilevar(f1,"ps")) then
        PS=f1->ps(n,:,:)
      end if
      hybm=(f1->hybm(:))
      print("dimsize PS: "+dimsizes(PS))
      do k=1,nlev-1,2
        p(:,:)=  hyam(k)*ps0 + hybm(k)*PS(:,:)
        print("k="+k+" p(k) min/max " + min(p(:,:)/100)+" "+max(p(:,:))/100)
      end do
      delete(ps0)
    end if
    end if


    
;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;     gsn_define_colormap(wks,"BlWhRe")
    
;  gsn_define_colormap(wks,"tbr_var_0-500")
;  res@gsnSpreadColorStart = 60           ; at color index 24.
    
;  gsn_reverse_colormap(wks)
; Any time you don't want a plot drawn, just set the resources:
;  res@gsnDraw  = False
;  res@gsnFrame = False
;
    
    res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnPaperOrientation = "portrait"
    res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 50           ; at color index 24.
    
    
    res@cnFillOn            = True         ; Turn on contour fill
;   res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
    res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
;    res@cnFillMode          = "CellFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLinesOn           = False     ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    
    res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
    res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
    
    
    
    res@mpOutlineOn          = True
    res@mpFillOn             = False
;    res@mpShapeMode          = "FreeAspect"


    if ( projection .eq.1) then     
       res@mpProjection      = "Orthographic"   
       res@mpCenterLatF      =  90.
       res@mpCenterLonF         = 0.
    end if
    if ( projection .eq. 2) then    
       res@mpMinLatF = -40
       res@mpMaxLatF = 80
       res@mpMinLonF = -180.
       res@mpMaxLonF =  -60.
    end if
    if ( projection .eq. 3) then    
       res@mpMinLatF = -50
       res@mpMaxLatF = 50
       res@mpMinLonF = 90.
       res@mpMaxLonF =  180.
    end if
    if ( projection .eq. 4) then    
       res@mpMinLatF = -40.
       res@mpMaxLatF = 15.
       res@mpMinLonF = -100.
       res@mpMaxLonF =  -40.
    end if
    if ( projection .eq. 5) then    
       res@mpMinLatF = 0.
       res@mpMaxLatF = 60.
       res@mpMinLonF = 50.
       res@mpMaxLonF = 110.
    end if
;   res@mpCenterLonF         = -180
;    res@mpMinLatF = -45
;    res@mpMaxLatF = 60
;    res@mpMinLonF = -180-90
;   res@mpMaxLonF = -180+120
; US 
;    res@mpMinLatF = -30
;    res@mpMaxLatF = 75
;    res@mpMinLonF = -180
;    res@mpMaxLonF =  0

;   ZOOM into region in PolvaniMWR2004:
;  res@gsnAddCyclic = False                  ; This is not global data.
;    res@mpProjection = "CylindricalEquidistant"
;    res@mpLimitMode = "LatLon"
;    res@mpMinLatF = 20
;    res@mpMaxLatF = 80
                                ;res@vpWidthF            = 1           ; change aspect ratio of plot






;    res@vpHeightF           = .3
  if (res@cnLevelSelectionMode .eq. "ManualLevels") then
    print("Manual Levels: [" + res@cnMinLevelValF +","+res@cnMaxLevelValF +"] spacing="+res@cnLevelSpacingF)
  else
    print("Using auto levels")
  end if


    

    
    do k=firstlev,lastlev,levinc
      if (isfilevar(f1,"hyam")) then
         print("ndim="+ndim+" k="+k+" eta="+(f1->hyam(k)+f1->hybm(k)))
      end if

      if(ndim.eq.2) then
        varmin = 1.e9
        varmax = -1.
        var(:,:)= (f1->$field_read$(:,:))
        copyatt( var, f1->$field_read$(:,:) )
      end if
      if(ndim.eq.3) then
        if (firstdim_is_time) 
          fi=n
        else
          fi=k
        end if
        print("NOTE: ndim=3 plotting first index = "+fi)
        varmin = 1.e9
        varmax = -1.
        var(:,:)= (f1->$field_read$(fi,:,:))
        copyatt( var, f1->$field_read$(fi,:,:) )
        if (sumprec) then
          var(:,:) = var(:,:) + f1->PRECL(fi,:,:)
        end if
        if (compute_vor) then
          U=var
          V = f1->$field_read2$(fi,:,:)
        end if
        if (compute_speed) then
          var(:,:) = sqrt( (f1->$field_read$(fi,:,:))^2 + (f1->$field_read2$(fi,:,:))^2 )
        end if
      end if
      if(ndim.eq.4) then
        if isStrSubset(field,"ICLDIWP") then
          var = 0
          do kk=0,nlev-1
;          sum over levels:  
            var(:,:) = var(:,:) +  f1->$field_read$(n,kk,:,:)
          end do
          copyatt( var, f1->$field_read$(n,0,:,:) )
        else
          if (pnew .gt. 0 ) then
            print("INTERPOLATING to pressure level = "+pnew) 
            varlev = f1->$field_read$(n,:,:,:)
            interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
            extrap = True  ; is extrapolation desired if data is outside the range of PS
            P0mb = 1000    ; ps in Pa, but this argument must be in mb
            var=vinth2p(varlev,f1->hyam,f1->hybm,pnew,PS(:,:),interp,P0mb,1,extrap)          
            copyatt( var, f1->$field_read$(n,k,:,:) )
          else
            var(:,:)=(f1->$field_read$(n,k,:,:))
            copyatt( var, f1->$field_read$(n,k,:,:) )
          end if
          if (pscale .gt.0) then
            var=var*pscale
          end if
        end if
      end if

      if (compute_vor) then
        print("computing vor via spherical harmonic transform")
        if (gauss_grid .eq. 1) then
          vor = uv2vrG_Wrap (U,V)
          ;div = uv2dvG_Wrap (U,V)
        else 
          if (gauss_grid .eq. 2) then
            vor = uv2vrF_Wrap (U,V)
            ;div = uv2dvF_Wrap (U,V)
          else 
            if (gauss_grid .eq. 3) then
              fo2fshv_Wrap(U,V,unew,vnew)
              vor = uv2vrF_Wrap(unew,vnew)
              ;div = uv2dvF_Wrap(unew,vnew)   
            else
              print("ERROR: bad gauss_grid="+gauss_grid)
            end if
          end if
        end if
        var=vor
        delete(vor)  
        delete(U)
        delete(V)
      end if

      var=var*scale
      if (longname .ne. "") 
          var@long_name = longname
      end if
      if (units .ne. "") 
          var@units = units
      end if

      print("name = "+f1name+" field="+field+" time="+time(n))
      print ("n,k="+n+" "+k+" min/max = "+min(var(:,:))+", "+max(var(:,:)) )
      print ("range lon: "+res@mpMinLonF+"-"+res@mpMaxLonF)
      print ("range lat: "+res@mpMinLatF+"-"+res@mpMaxLatF)

;    res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
;    res@cnLevelSpacingF = .15
     ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 
      res@tiMainString        = f1name
      if (ndim .eq. 3) then
        res@tiMainString        = f1name + sprintf(" %6.3f", time(n))+" days"    ;
;        res@tiMainString        = f1name + " level="+k
      end if 
      if (ndim .eq. 4) then
        res@tiMainString        = f1name+" level = "+k+" "+time(n)+"days"  
      end if
      if (pnew .gt. 0) then
        res@tiMainString        = f1name+" p="+pnew+"mb "+time(n)+"days"  
      end if
      if (isfilevar(f1,"lev1")) then
        res@tiMainString        = f1name+" p="+f1->lev1(k)+"mb "+time(n)+"days"  
      end if

      ;res@tiMainString=" "
      
      res@gsnDraw        = False
      res@gsnFrame       = False 

      contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
;    contour = gsn_contour(wks,var(:,:),res)     ; Draw a contour plot.  ;
      print("Contour complete: "+res@tiMainString)

     if ( False ) then  ; overlay wind vector plot
       res_vc = res
       res_vc@vcGlyphStyle               = "LineArrow"
       res_vc@vcLineArrowThicknessF      = 5
       res_vc@vcMinDistanceF             = 0.01
       res_vc@vcRefLengthF               = 0.03

       vector  = gsn_csm_vector(wks,f[0]->UBOT(n,:,:),f[0]->VBOT(n,:,:),res_vc)
       overlay(contour,vector)
     end if 
     draw(contour)
     frame(wks)





      if (plot_spectra) then
        print("computing spectra")
        if (gauss_grid .eq. 1) then
          var_hat = shagC (var(:,:))
        else 
          if (gauss_grid .eq. 2) then
            var_hat = shaeC (var(:,:))
          else 
            if (gauss_grid .eq. 3) then
              ; interpolate from fixed-offset to fixed grid
              ; areg = fo2fshv(aoff)
              ; interpolate from fixed-offset to fixed grid
              ; for vectors: fo2fshv_Wrap(uoff,voff,ureg,vreg)
              var_hat = shaeC (fo2fsh(var(:,:)) )
            else
              print("ERROR: bad gauss_grid="+gauss_grid)
              var_hat = 0
            end if
          end if
        end if
        

        spc = new ( (/nlat/), "double", 1d-99 )
        spc = spc@_FillValue
        spc!0 = "wavenumber"

        tmpdim = dimsizes(var_hat)
        nmax = tmpdim(2)
        pwr = var_hat(0,:,:)*var_hat(0,:,:) 
        do n1=0,nmax-1
          x(n1) = n1
          spc(n1) = (/pwr(n1,0) /)
          do m=1,n1
            spc(n1) = spc(n1) + 2.*(pwr(n1,m))
          end do
                                ;print(n1+" power="+spc(n1))
        end do
        
        athick=4
        lthick=2
        
        resxy=True
                                ;res@gsLineThicknessF = 3.0
        resxy@xyXStyle ="Log"
        resxy@xyYStyle ="Log"
        resxy@tmXBMinorPerMajor = 8   
        resxy@tmXTMinorPerMajor = 8   
        resxy@tmYRMinorPerMajor = 8   
        resxy@tmYLMinorPerMajor = 8   
        resxy@tmXMajorGridThicknessF = athick
        resxy@tmYLLabelFont = 21    ;  21 = helvetica
        resxy@tmXBLabelFont = 21    ;  22 = helvetica-bold 
        
        resxy@tmXBMajorThicknessF = athick
        resxy@tmYLMajorThicknessF = athick
        resxy@tmXBMinorThicknessF = athick 
        resxy@tmYLMinorThicknessF = athick 
        resxy@tmBorderThicknessF = athick
                                ;resxy@tiMainString           = ":F21:CLIMATOLOGY"
        resxy@tiXAxisString          = ":F21:spherical wavenumber"
        resxy@tiYAxisString          = ":F21:Energy"
        resxy@gsnFrame = False
        resxy@xyLineThicknessF = lthick
        resxy@xyDashPattern = 0
        
        resxy@trXMinF = 1       
        resxy@trXMaxF = 1e3
        resxy@trYMinF = 1e-1      
        resxy@trYMaxF = 1e8
;        plot  = gsn_xy (wks,x,spc,resxy) ; create plot
;        frame(wks)                                     ; now advance
        spc_all(fn,0:nlat-1)=spc(0:nlat-1)
        delete(spc)
        delete(var_hat)
        delete(pwr)
      end if
    end do  ; vertical levels
    delete(var)
    delete(p)
  end do ; loop over time level
  delete(time)
; loop over files
end do 


if ( plot_spectra) then
  resxy@xyLineColors = (/"black","purple","blue","green","orange","red3","cyan","brown"/)

  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = 1e-1      
  resxy@trYMaxF = 1e8
  plot  = gsn_xy (wks,x(0:kmax-1),spc_all(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance

  if (fcnt.gt.1)  then
  spc_r = spc_all(:,:)
  do n1=0,fcnt-1
    do j=0,nlat_max-1
      spc_r(n1,j)=0
      if (.not. ismissing(spc_all(0,j)) .and. spc_all(0,j) .gt. 0) then
        spc_r(n1,j) = spc_all(n1,j)/spc_all(0,j)
      end if
    end do
    print("kf="+n1+" data="+x(0)+" "+spc_all(n1,0))
  end do
  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = .01
  resxy@trYMaxF = 2e0
  plot  = gsn_xy (wks,x(0:kmax-1),spc_r(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance
  end if 

  do n1=0,nlat_max-1
    spc_all(:,n1) = spc_all(:,n1)*(n1*(n1+1))
  end do
  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = 1e6      
  resxy@trYMaxF = 1e10
  plot  = gsn_xy (wks,x(0:kmax-1),spc_all(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance

end if


delete(wks)
end


