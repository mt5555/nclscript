load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  Mark Taylor 2011
;
;  ke.ncl:
;  read u,v files, compute spectra, plot at snapshots
;  save average spectra to spectra.nc
;  run this on one file, or multiple files from the same simulation
;  as it will average across files
;
;  ke2.ncl:  
;     plots the average spectra as computed by ke.ncl
;     publication quality
;
;
;
;  set filenames below, then run with:
;   ncl ke.ncl 'plvl=250'  't1=31.0'  't2=58.0'
;  plvl = pressure level in mb.  default is 200)
;  [t1,t2] = compute spectra for all times in this range. set correct timeinc below.  
;
;  filenames should all be from the same simulation, but different times
;  run this code once for each different simulation/case
;



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
    ;print(d1-d2)
    if (abs(d1-d2) .lt. 1e-8 ) then
      print ("detected equal angle grid without poles")
      gauss_grid = 3    ; need to interpolate using fo2fsh
    else
      print("detected Gauss grid")
      gauss_grid = 1  
    end if
  end if
  return(gauss_grid)
end



begin
plot_each_time = False  ; set to True to plot spectra from every snapshot
                        ; in addition to time mean 


gname="ke"+systemfunc("echo $$")
output_format = "ncgm"
wks = gsn_open_wks(output_format,gname)
print ("Writing file "+gname+"."+output_format)




; NETCDF file containing U,V  or U250, V250
fnames =  systemfunc( \
 "ls ./*.cam.h2.0001-??-??-10800.latlon.nc; " \
)
if (any(ismissing(fnames))) then
   exit
end if

fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
print("files = "+fnames)
fcnt = dimsizes(fnames)
print ("number of uv input files = " + fcnt)



; check command line argument for pnew:
pnew=250
if (isvar("plvl")) then         ; pressure level on command line?
  pnew = plvl                   ; command line 
end if

; change the levels and or times to plot here
firsttime=32d0
lasttime=32d0
timeinc=0.125d0

if (isvar("t1")) then      
  firsttime=t1
end if
if (isvar("t2")) then      
  lasttime=t2
end if
if (isvar("tinc")) then      
  timeinc=tinc
end if

firsttime_found=999999d0
lasttime_found=-1d0
count = 0

print("Computing spectra at P = "+pnew)
if (isfilevar(fuv_all[0],"U"+pnew)) then
  useVORDIV=False
  hasUplev = True
  usePSmean = True  ; actually PS is not used
  Uplev = "U"+pnew
  Vplev = "V"+pnew
  print("found "+Uplev+" on input file. no vertical interpolation needed")
else
   hasUplev = False
   if (isfilevar(fuv_all[0],"VOR")) then
      useVORDIV=True
      U3d="VOR"
      V3d="DIV"
   else
      useVORDIV=False
      U3d="U"
      V3d="V"
   end if
   usePSmean = True
   if (isfilevar(fuv_all[0],"PS")) then
      sd=dimsizes(fuv_all[0]->PS)
      if (sizeof(sd).eq.3) then
         ; needs to be time, lat, lon
         usePSmean = False
      end if
   end if
end if

if (.not. hasUplev) then
   if (usePSmean) then
      print("WARNING: USING PS=1000mb to interpolate to plvl")
   else
      print("using PS for interpolation to 250mb")
   end if
end if

f1=fuv_all[0]
gauss_grid=get_gridtype(f1->lat)




fuv=fuv_all[0];    file with u, v

ntimes = dimsizes(fuv->time)
nlev = dimsizes(fuv->lev)
nlat = dimsizes(fuv->lat)
nlon = dimsizes(fuv->lon)
print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)


;
; for output to spectra.nc
;
fnames_list = new ( (/fcnt,800/), "character")
fnames_list!0 = "run"
spc = new ( (/nlat/), "double", 1d-99 )
spc = spc@_FillValue
spc!0 = "wavenumber"
x=spc
spc_vr = spc
spc_dv = spc

spc_avg = new ( (/nlat/), "double", 1d-99 )
spc_avg!0 = "wavenumber"
spc_avg=0
spc_avg_vr=spc_avg
spc_avg_dv=spc_avg


x_avg=spc_avg


; single snapshot spectra defaults:
res=True
                                ;res@gsLineThicknessF = 3.0
res@xyXStyle ="Log"
res@xyYStyle ="Log"
res@tmXBMinorPerMajor = 8   
res@tmXTMinorPerMajor = 8   
res@tmYRMinorPerMajor = 8   
res@tmYLMinorPerMajor = 8   
res@tmYLLabelFont = 21    ;  21 = helvetica
res@tmXBLabelFont = 21    ;  22 = helvetica-bold 

res@gsnFrame = False
res@xyLineThicknessF = 1.0
 
colors = (/"white","black","OrangeRed4","Red3",\
"red1","Red3","Red4","Blue",\
"OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
"red1","red2","red3","red4",\
"DarkTurquoise","Turquoise","Cyan","blue",\
"gold","greenyellow",\
"GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
/)
gsn_define_colormap(wks, colors)
 
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
res@xyDashPattern = 0

polyres                  = True
polyres@gsLineThicknessF = 1.0
polyres@gsLineColor = "black"

   



; loop over all times
do time=firsttime,lasttime,timeinc
  ;print("looking for t="+time)
  ; for each time, loop over files looking for the one with this time
  do fk=0,fcnt-1        
    fuv=fuv_all[fk];
    nuv = ind( (fuv->time).eq.time)

    if ( ismissing(nuv)) then
      print("skipping = "+fnames(fk)+" record=(missing)" )
    else
      count=count+1
      print(" file = "+fnames(fk))
      print("time="+time+" record="+nuv+" count="+count )

      if (time .lt. firsttime_found) then
        firsttime_found=time
      end if                      
      if (time .gt. lasttime_found) then
        lasttime_found=time
      end if                      
      
      ; lat/lon data
      deg = 360.0 / nlon
      NE = 1
      if (isatt(fuv,"ne")) then
        if (fuv@ne .gt. 0) then
          NE=fuv@ne
          if (isatt(fuv,"np")) then		
             NPTS = fuv@np
	  else
             NPTS = 4
	  end if
          deg=360.0 / (4*NE*(NPTS-1))
          hdeg=sprintf("%4.2f", deg)
          tmax = doubletoint(180/deg)
          print("NE="+NE+" NPTS="+NPTS+" tmax="+tmax+" degree="+hdeg)
        end if
      end if
      
      

      if ( hasUplev ) then
        var = fuv->$Uplev$(nuv,:,:)
        var2 = fuv->$Vplev$(nuv,:,:)
      else
        if (usePSmean ) then
	  ps = fuv->$U3d$(nuv,0,:,:)*0 + 100000
	  ps@units="Pa"
	else
          ps = fuv->PS(nuv,:,:)
          copyatt( ps, fuv->PS(nuv,:,:))
        end if	    
        if ( "mb" .eq. ps@units) then
          print("converting ps to Pascals. orig units ="+ps@units)
          ps=ps*100
          ps@units="Pa"
        end if
        
        hyam = fuv->hyam
        hybm = fuv->hybm
        
        interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
        extrap = True  ; is extrapolation desired if data is outside the range of PS
        P0mb = 1000    ; ps in Pa, but this argument must be in mb
        
        temp=vinth2p(fuv->$U3d$(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
	var=temp(0,:,:)
        temp=vinth2p(fuv->$V3d$(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
	var2=temp(0,:,:)
      end if
      print ("var1 min/max = "+min(var(:,:))+", "+max(var(:,:)) )
      print ("var2 min/max = "+min(var2(:,:))+", "+max(var2(:,:)) )

      vor_grid=gauss_grid
      if (  useVORDIV  ) then
         ; var = VOR.  var2= DIV
	 vr=var
	 dv=var2
      else
        print("converting U,V to vor,div")
        if (gauss_grid .eq. 1 ) then
          vr = uv2vrG_Wrap (var,var2)
          dv = uv2dvG_Wrap (var,var2)   
        else 
          if (gauss_grid .eq. 2 ) then
            vr = uv2vrF_Wrap (var,var2)
            dv = uv2dvF_Wrap (var,var2)   
          else
            if (gauss_grid .eq. 3 ) then
              ; interpolate from fixed-offset to fixed grid
              fo2fshv_Wrap(var,var2,unew,vnew)
              vr = uv2vrF_Wrap (unew,vnew)
              dv = uv2dvF_Wrap (unew,vnew)
              vor_grid=2   
            else
              print("Error: unsupported grid for spherical harmonic transform")
            end if
          end if
        end if
      end if
      
      ; compute spectra from vor/div
      if (vor_grid .eq. 1 ) then
        ab_vr = shagC (vr(:,:))       ; GAUSS grid
        ab_dv = shagC (dv(:,:))       ; GAUSS grid
      else 
        if (vor_grid .eq. 2 ) then
          ab_vr = shaeC (vr(:,:))       ; equally spaced
          ab_dv = shaeC (dv(:,:))       ; equally spaced
        else
          print("Error: unsupported grid for spherical harmonic transform")
        end if
      end if




      tmp=dimsizes(ab_dv)
      nmax = tmp(1)
      mmax = tmp(2) 
; compute energy in vort and div components, no factor 1/2
;
      cr = ab_vr(0,:,:)                ; real coef  (nlat,nlat) vort
      ci = ab_vr(1,:,:)                ; imaginary  (nlat,nlat) vort
      pwr_vr = (cr^2 + ci^2)        ; (nlat,nlat)  array
      cr = ab_dv(0,:,:)                ; real coef  (nlat,nlat) div
      ci = ab_dv(1,:,:)                ; imaginary  (nlat,nlat) div
      pwr_dv = (cr^2 + ci^2)       ; (nlat,nlat)  array
      
      do n1=0,nlat-1
	x(n1) = n1
        spc_vr(n1) = (/pwr_vr(n1,0) /)
        spc_dv(n1) = (/pwr_dv(n1,0) /)
        factor = 4. * n1*(n1+1)
        do m=1,min( (/n1,mmax-1/) )
          spc_vr(n1) = spc_vr(n1) + 2.*(pwr_vr(n1,m))
          spc_dv(n1) = spc_dv(n1) + 2.*(pwr_dv(n1,m))
        end do
        if (n1 .ne. 0) then
          spc_vr(n1) = (6371229.^2) / factor * spc_vr(n1)  ; scaling from NCAR/TN-388+STR
          spc_dv(n1) = (6371229.^2) / factor * spc_dv(n1)  ; scaling from NCAR/TN-388+STR
        end if
      end do
      ; normalize so it agrees with Vector Spherical Harmonics code
      spc_vr = spc_vr/2
      spc_dv = spc_dv/2


      ; since we now compute vor & div, no reason to also compute
      ; vector harmonic transform.  
      spc = spc_vr + spc_dv 


      ; COMPUTE SPECTRA FROM U,V with vector spherical harmonics
      ; no longer needed, but save this code to verify it gives the
      ; same answer as computing 'spc' from vor/div spectra
      if (False) then	
        print ("computing vector spherical harmonic transform...t="+time)
        if (gauss_grid .eq. 1 ) then
          ab = vhagC(var(:,:),var2(:,:));   ; GAuss grid
        else if (gauss_grid .eq. 2) then
          ab = vhaeC(var(:,:),var2(:,:));  ; for equally spaced grid, including both poles
        else
           print("Error: unsupported grid for spherical harmonic transform")
        end if
        end if
      
                                ; compute energy in u,v components
        cr = ab(0,:,:)                ; real coef  (nlat,nlat)
        ci = ab(1,:,:)                ; imaginary  (nlat,nlat)
        pwr = (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
        cr = ab(2,:,:)                ; real coef  (nlat,nlat)
        ci = ab(3,:,:)                ; imaginary  (nlat,nlat)
        pwr = pwr + (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
        ab=0
      
                                ; for clarity use do loops
        do n1=0,nlat-1
          spc(n1) = (/pwr(n1,0)/)  ; so we preserve attributes in spc
          do m=1,min( (/n1,mmax-1/) )
            spc(n1) = spc(n1) + 2.*pwr(n1,m)
          end do
          spc(n1) = 0.25*spc(n1)
        end do

;        print("sum of spc = "+sum(spc))
;        fgw = addfile("~/scratch1/cam/eul2/camrun.cam2.h0.0000-09.nc","r")
;        gw=fgw->gw
;        varl = ( gw # (var^2 + var2^2) )/sum(gw)
;        s2 =  dimsizes(var) 
;        output = sum(varl)/s2(1)
;        output = output/2
;        print(" .5u^2 = "+output)
        
      end if



      spc_avg(:)= ((count-1)*spc_avg(:) + spc(:) )/(count)
      spc_avg_vr(:)= ((count-1)*spc_avg_vr(:) + spc_vr(:) )/(count)
      spc_avg_dv(:)= ((count-1)*spc_avg_dv(:) + spc_dv(:) )/(count)
      x_avg(:) = x(:)        

      
      
      
      if ( plot_each_time .eq. True ) then

        print("Plotting time="+time)

        res@tiMainString           = ":F21:DAY "+time
        res@tiXAxisString          = ":F21:spherical wavenumber"
        res@tiYAxisString          = ":F21:Kinetic Energy"
        plot  = gsn_xy (wks,x,spc,res) ; create plot

        xref = (/10.,300./)
        spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
        gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
        frame(wks)                                     ; now advance
        
        spc2=spc
        power=5./3.
        do n1=0,nlat-1
          spc2(n1) = spc(n1)*(n1^power)
        end do
        res@tiYAxisString          = ":F21:Compensated Kinetic Energy (5/3)"
        res@trYMinF = 1e0
        res@trYMaxF = 1e4
        plot  = gsn_csm_xy (wks,x,spc2,res) ; create plot
        xref = (/10.,300./)
        spc_ref = (xref^power)*100*(xref/10)^(-power) ; make it 6 @ k=15
        gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
        frame(wks)                                     ; now advance 


        power=3.
        do n1=0,nlat-1
          spc2(n1) = spc(n1)*(n1^power)
        end do
        res@tiYAxisString          = ":F21:Compensated Kinetic Energy (3)"
        res@trYMinF = 1e2
        res@trYMaxF = 1e6
        plot  = gsn_csm_xy (wks,x,spc2,res) ; create plot
        xref = (/10.,300./)
        spc_ref = (xref^power)*100*(xref/10)^(-power) ; make it 6 @ k=15
        gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
        frame(wks)                                     ; now advance 


      end if
      
      break    
    end if  ; compute specta, plot
    
    
  end do ; loop over files, looking for time
end do  ; loop over time level

if (firsttime_found .eq. 999999d0)
   print("no times found in ["+firsttime+","+lasttime+"]")
   exit
end if


res@trYMinF = 1e-7      
res@trYMaxF = 1e3   
res@tiMainString     = "Climatology"
res@xyLineColors = (/"blue","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg_dv,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg_vr,res) ; create plot

xref = (/10.,300./)
spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
frame(wks)


;save spectrum data
   fk=0
   if (  useVORDIV  ) then
     outname = "spectra-vordiv"+pnew+"tmax"+(nlon/2)+"t"+firsttime_found+"-"+lasttime_found+".nc"
   else
     outname = "spectra-uv"+pnew+"tmax"+(nlon/2)+"t"+firsttime_found+"-"+lasttime_found+".nc"
   end if
   system("/bin/rm "+outname)           ; remove any pre-existing file
   ncdf     = addfile(outname ,"c")  ; open output netCDF file
   ncdf->x_avg = x_avg(:)
   ncdf->spc_avg = spc_avg(:)
   ncdf->spc_avg_vr = spc_avg_vr(:)
   ncdf->spc_avg_dv = spc_avg_dv(:)
   ncdf->ne = NE
   ncdf->deg = deg
   ncdf->count = count
do fk=0,fcnt-1
   c=sizeof(stringtochar(fnames(fk)));
   fnames_list(fk,0:c-1)=stringtochar(fnames(fk));
end do
   ncdf->fnames = fnames_list


end

