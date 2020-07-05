load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;***********************************************************************;
; 
; run via the tcsh script in ~/bin:
;    contour_latlon.sh  PRECT *.nc
;
; or:
;     ncl 'field="PRECT"' 'files="*.nc"'  contour_latlon.ncl
; or
;     ncl contour_latlon.ncl
;     and edit field and files variables below
;
; removed code to generate triangle mesh from GLL subcell data, doesn't work  
; see test_elementnodes.ncl
;***********************************************************************;
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
print("fnames = "+fnames)



; default is to read field specified.  but some fields, like
; SPEED, we read other variables and compute derived quantities like SPEED:
field_read=field
compute_speed=False     
compute_theta=False      
sumprec = False
compute_vor = False




;output_format = "ncgm"
output_format = "pdf"
;output_format = "X11"
res=True
res@cnLevelSelectionMode = "AutomaticLevels"
scale=1d0

; default is to plot level 18. especially for OMEGA
; override below if desired
firstlev_read=18
lastlev_read=18
levinc=1

; timelevels to plot:
firsttime=0
lasttime_max=999999
timeinc=1


    ; US 
    res@mpMinLatF = -30
    res@mpMaxLatF = 50
    res@mpMinLonF = -150
    res@mpMaxLonF =  -80


;     res@mpMinLatF = 20.
;     res@mpMaxLatF = 50.
;     res@mpMinLonF = -180.
;     res@mpMaxLonF = -130.
     
     ;res@mpProjection      = "Orthographic"
;       res@mpCenterLatF      =  80.
;       res@mpCenterLonF         = 296.
;     res@mpLimitMode="Angles"
;       res@mpLeftAngleF=hdeg
;      res@mpRightAngleF=hdeg
;      res@mpBottomAngleF=hdeg
;      res@mpTopAngleF=hdeg

       


  
if (field .eq. "PRECT" ) then
; convert from m/s to mm/day
  scale=1000d0*(24*3600)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  0.
  res@cnMaxLevelValF = 10.
  res@cnLevelSpacingF =  0.25
  if (.not.isfilevar(f[0],"PRECT")) then
    ; PRECT = PRECL + PRECC
    ; read PRECC, then read PRECL and sum them below:
    sumprec=True        
    field_read = "PRECC"
    field_read2 = "PRECL"
  end if
end if
if (field .eq. "PRECC" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  0.
  res@cnMaxLevelValF =  18 / 1000.  / (24*3600.)
  res@cnLevelSpacingF =  res@cnMaxLevelValF/25.
end if
if (field .eq. "PSL" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  96000.
  res@cnMaxLevelValF = 103000.
  res@cnLevelSpacingF =  200.
end if
if (field .eq. "TS" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  200.
  res@cnMaxLevelValF = 300.
  res@cnLevelSpacingF =  5.
end if
if (field .eq. "TBOT" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  160.
  res@cnMaxLevelValF = 300.
  res@cnLevelSpacingF =  5.
end if
if (field .eq. "POTTBOT" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  160.
  res@cnMaxLevelValF = 300.
  res@cnLevelSpacingF =  5.
  
  compute_theta=True
  field_read="TBOT"
end if
if (field .eq. "UBOT" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -20.
  res@cnMaxLevelValF = 20.
  res@cnLevelSpacingF =  2.
end if
if (field .eq. "VBOT" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -20.
  res@cnMaxLevelValF = 20.
  res@cnLevelSpacingF =  2.
end if
if (field .eq. "SPEEDBOT" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  0.
  res@cnMaxLevelValF = 25.
  res@cnLevelSpacingF =  1.
  
  compute_speed=True
  field_read="UBOT"
  field_read2="VBOT"
end if
if (field .eq. "QRL" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -0.0005
  res@cnMaxLevelValF = 0.0005
  res@cnLevelSpacingF =  .0005/40
end if
if (field .eq. "DYN_THTEND" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -2.
  res@cnMaxLevelValF = 2.
  res@cnLevelSpacingF =  2./40
end if
if (field .eq. "DYN_TTEND" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -0.01
  res@cnMaxLevelValF = 0.01
  res@cnLevelSpacingF =  .01/40
end if
if (field .eq. "PTTEND" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -0.01
  res@cnMaxLevelValF = 0.01
  res@cnLevelSpacingF =  .01/40
end if
if (field .eq. "DTCORE" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  -0.01
  res@cnMaxLevelValF = 0.01
  res@cnLevelSpacingF =  .01/40
end if
if (field .eq. "PS" ) then
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF =  90000.
  res@cnMaxLevelValF = 100000.
  res@cnLevelSpacingF =  200.
end if
if (field .eq. "OMEGA" ) 
;; monthly means: Pa/s
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -.2
  res@cnMaxLevelValF = .2
  res@cnLevelSpacingF = .01
  
  res@cnMinLevelValF = -.2
  res@cnMaxLevelValF = .2
  res@cnLevelSpacingF = .01
  
  res@cnMinLevelValF = -5.
  res@cnMaxLevelValF = 5.
  res@cnLevelSpacingF = .25
end if
if (field .eq. "OMEGA500" ) 
  field = "OMEGA500"
  res@cnLevelSelectionMode = "ManualLevels"
;; monthly means: Pa/s
;     res@cnMinLevelValF = -.2
;     res@cnMaxLevelValF = .2
;     res@cnLevelSpacingF = .01
;; monthly means 850mb
;    pnew=850.
;    res@cnMinLevelValF = -.2
;    res@cnMaxLevelValF =  .2
;    res@cnLevelSpacingF = .01
;; snapshots
;    res@cnMinLevelValF = -1.0
;    res@cnMaxLevelValF = 1.0
;    res@cnLevelSpacingF = .05
  res@cnMinLevelValF = -5.0
  res@cnMaxLevelValF = 5.0
  res@cnLevelSpacingF = .2
end if
if (field .eq. "DIV")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -.004
  res@cnMaxLevelValF = .004
  res@cnLevelSpacingF = .004/40
end if
if (field .eq. "VOR")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -.001
  res@cnMaxLevelValF = .001
  res@cnLevelSpacingF = .001/40
end if
if (field .eq. "VOR250" ) then
  compute_vor = True
  field_read="U250"
  field_read2="V250"
end if
if (field .eq. "SHFLX")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -500.
  res@cnMaxLevelValF = 500.
  res@cnLevelSpacingF = 20.
end if
if (field .eq. "LHFLX")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -500.
  res@cnMaxLevelValF = 500.
  res@cnLevelSpacingF = 20.
end if
if (field .eq. "FLNS")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -200.
  res@cnMaxLevelValF = 200.
  res@cnLevelSpacingF = 2.
end if
if (field .eq. "Q")
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = 0.0
  res@cnMaxLevelValF = .005
  res@cnLevelSpacingF = .005/50.
  
  res@cnMaxLevelValF = .001
  res@cnLevelSpacingF = .001/50.
end if
if (field .eq. "PHIS")
;    res@cnLevelSelectionMode = "ManualLevels"
;    res@cnMinLevelValF = -500
;    res@cnMaxLevelValF = 55000
;    res@cnLevelSpacingF = 1000
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

latname = "lat"
lonname = "lon"

if (isfilevar(f[0],field)) then
  dimnames = getvardims(f[0]->$field$)
  print("variable dimension names: "+dimnames)
  if (dimsizes(dimnames).eq.3) then
    if ( dimnames(2).eq."ncol_d" ) then
      print("detected ncol_d dimension. using lat_d, lon_d coordiantes.")
      latname = "lat_d"
      lonname = "lon_d"
    end if
  end if
end if   

lat1d=doubletofloat(f[0]->$latname$);*180.0/3.14
lon1d=doubletofloat(f[0]->$lonname$);*180.0/3.14
print("lat1d min/max = "+min(lat1d)+" "+max(lat1d))
print("lon1d min/max = "+min(lon1d)+" "+max(lon1d))
ncols = dimsizes(lat1d)

fcnt = dimsizes(fnames)
if (isfilevar(f[0],"time")) then
   ntimes = dimsizes(f[0]->time)
  time = f[0]->time
else
   ntimes = 1
   time = 0
end if
if (isfilevar(f[0],"lev")) then
   nlev = dimsizes(f[0]->lev)
else
   nlev =1
end if
nxy = dimsizes(lat1d)


if (ntimes-1 .le. lasttime_max) then
  lasttime=ntimes-1
else
  lasttime=lasttime_max
end if


wks = gsn_open_wks(output_format,field+"native")
print ("Writing file "+field+"."+output_format)
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)


res@sfXArray            = lon1d        ; Required to tell NCL where to
res@sfYArray            = lat1d        ; overlay data on globe.
res@cnFillMode          = "RasterFill"   ; Style of fill. You can also

; PG2 data: use cell boundaries from NE30pg2.scrip file:
;    res.cnFillMode = 'CellFill'
;    res.sfXCellBounds = grid_corner_lon
;    res.sfYCellBounds = grid_corner_lat

res@gsnMaximize         = True         ; Maximize size of plot in frame
res@gsnPaperOrientation = "portrait"
res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 24           ; at color index 24.
res@cnFillOn            = True         ; Turn on contour fill
res@cnLinesOn           = False        ; Turn off contour lines
res@cnLineLabelsOn      = False        ; Turn off contour line labels
res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
res@mpOutlineOn          = True
res@mpFillOn             = False
res@mpShapeMode          = "FreeAspect"

if (isfilevar(f[0],"P0")) then
  ps0 = (f[0]->P0)
else
  ps0 = 1000
end if

var = new((/nlev,nxy/),"double")
p = new((/nlev,nxy/),"double")
print ("ntimes="+ntimes+" nlev="+nlev+" nxy="+nxy)

ndim = dimsizes(getfilevardimsizes(f[0],field_read))
; ndim = 1   assume ncols
; ndim = 2   assume time,ncols
; ndim = 3   assume time,nlev,ncols

;  print("sum area = "+sum(f[0]->area(:))) 
;  print("sum FLNT = "+sum(f[0]->area(:)*f[0]->FLNT(0,:))/sum(f[0]->area(:))    ) 

print("field_read = "+field_read+" dimsizes:"+dimsizes(f[0]->$field_read$))
do n=firsttime,lasttime,timeinc
;  print(n+" date="+f[0]->date(n)+" "+f[0]->datesec(n)/3600.+" time="+time(n)+" days.  ndim="+ndim)

  if(ndim .eq. 1) then
    varmin = 1.e9
    varmax = -1.
    var(0,:)= (f[0]->$field_read$(:))
    firstlev=0
    lastlev=0
  end if
  if(ndim .eq. 2) then
    varmin = 1.e9
    varmax = -1.
    var(0,:)= (f[0]->$field_read$(n,:))
    if (sumprec) then
      var(0,:) = var(0,:) + f[0]->$field_read2$(n,:)
      var@longname = "PRECT"
    end if
    if (compute_speed) then
      var(0,:) = sqrt( (var(0,:))^2 + (f[0]->$field_read2$(n,:))^2 )
      var@longname = "SPEEDBOT"
    end if	      
    if (compute_theta) then
       ;     pm = pres_hybrid_ccm_se(f[0]->PS,ps0,f[0]->hyam,f[0]->hybm)  ; (time,lev,ncol)
       pm1=  f[0]->hyam(lastlev_read)*ps0 + f[0]->hybm(lastlev_read)*f[0]->PS(n,:)
       print("pressure layer 1 min/max = "+min(pm1)/100+" "+max(pm1)/100)
       exner = (pm1/ps0)^.2856
       print("exner layer 1 min/max = "+min(exner)+" "+max(exner))
       ; T = theta*exner
       var(0,:) = var(0,:)/exner(:)
    end if					

    firstlev=0
    lastlev=0
  end if
  if(ndim .eq. 3) then
    if isStrSubset(field,"ICLDIWP") then
      var = 0
      do k=0,nlev-1
;        sum over levels:  
        firstlev=0
        lastlev=0
        var(0,:) = var(0,:) +  f[0]->$field_read$(n,k,:)
      end do
    else
      ; for ne256, running out of memory unless we do this:
      if (firstlev_read.eq.lastlev_read) then
         ; so we dont run out of memory, read level into k=0
         firstlev=0
         lastlev=0
         var(0,:) = f[0]->$field_read$(n,firstlev_read,:)
      else
        do k=firstlev_read,lastlev_read,levinc
          print("read k="+k+" n="+n)
          var(k,:)=(f[0]->$field_read$(n,k,:))
        end do
      end if 
    end if
  end if
  print("var min/max " + min(var(0,:)) + " " + max(var(0,:)) )
  imin=minind(var(0,:))
  print("minind="+minind(var(0,:))+" field(minind)="+var(0,imin))
  print("minind lat/lon = "+lat1d(imin)+" " + lon1d(imin))
  var = var * scale

  if (ndim .ge. 2) then
    hyam=(f[0]->hyam(:))
    hybm=(f[0]->hybm(:))
;  do k=0,nlev-1
;    p(k,:)=  hyam(k)*ps0 + hybm(k)*f[0]->PS(n,:)
;    print("k="+k+" min/max " + min(p(k,:)/100)+" "+max(p(k,:))/100)
;  end do
  end if

  if (res@cnLevelSelectionMode .eq. "ManualLevels") then
    print("Manual Levels: [" + res@cnMinLevelValF +","+res@cnMaxLevelValF +"] spacing="+res@cnLevelSpacingF)
  else
    print("Using auto levels")
    res@cnMaxLevelCount = 50
  end if

  do k= firstlev,lastlev,levinc
    print("plotting min/max " + min(var(k,:)) + " " + max(var(k,:)) ) 

    if (k.eq.0) then
       res@tiMainString        = field+" time="+sprintf("%6.3f", time(n))+" days"    ;
    else
       res@tiMainString        = field+" at level = "+k+" time="+sprintf("%6.3f", time(n))+" days"    ;
     end if

    res@gsnFrame = False       ; or advance the frame
    contour = gsn_csm_contour_map(wks,var(k,:),res)     ; Draw a contour plot.
    print("Contour complete: "+res@tiMainString)
    print(" ")
    
    resp=True
    do i=0,ntris-1
    j0=trigrid(i,0)-1 ; convert to zero based for NCL arrays
    j1=trigrid(i,1)-1 ; convert to zero based for NCL arrays
    j2=trigrid(i,2)-1 ; convert to zero based for NCL arrays
    
    xlon=(/lon1d(j0),lon1d(j1),lon1d(j2),lon1d(j0)/)
    xlat=(/lat1d(j0),lat1d(j1),lat1d(j2),lat1d(j0)/)
    gsn_polyline(wks,contour,xlon,xlat,resp)
  end do ;  loop over levels
end do ; loop over time slices
end
