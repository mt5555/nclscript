load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  Mark Taylor 2011
;
;  ke-demo.ncl:
;  read u,v files, interpolate to 250mb, compute spectra, plot at snapshots
;


;***********************************************************************;
; detect type of lat/lon grid
;
; latlon_gridtype =-1    routine failed 
; latlon_gridtype = 1    Gauss grid (use shaeg routines)
; latlon_gridtype = 2    Equal angle grid with poles  (use shaec routines)
; latlon_gridtype = 3    Equal angle grid without poles  (must interpolate first)
;
;***********************************************************************;
function get_gridtype(lat)
local d1,d2
begin
  latlon_gridtype = -1
  if (lat(0) .lt. -89.99999 ) then
    print ("detected equal angle grid with poles")
    latlon_gridtype = 2    ; use shaec routines
  else
    d1=lat(1)-lat(0)
    d2=lat(2)-lat(1)
    ;print(d1-d2)
    if (abs(d1-d2) .lt. 1e-8 ) then
      print ("detected equal angle grid without poles")
      latlon_gridtype = 3    ; need to interpolate using fo2fsh
    else
      print("detected Gauss grid")
      latlon_gridtype = 1  
    end if
  end if
  return(latlon_gridtype)
end



begin




; NETCDF file containing U,V  or U250, V250
fnames =  systemfunc( \
 "ls *test.nc" \
)
if (any(ismissing(fnames))) then
   exit
end if

fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
fuv=fuv_all[0]


pnew=250  ; interpoalte to 250mb
latlon_gridtype=get_gridtype(fuv->lat)
ntimes = dimsizes(fuv->time)
nlev = dimsizes(fuv->lev)
nlat = dimsizes(fuv->lat)
nlon = dimsizes(fuv->lon)
print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

spc = new ( (/nlat/), "double", 1d-99 )
spc = spc@_FillValue
spc!0 = "wavenumber"
x=spc
spc_vr = spc
spc_dv = spc



wks = gsn_open_wks("pdf","ke-demo")
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

ntime=0  ; time index to use

;
;
; interpolate VOR/DIV to 250mb
;
;
ps=fuv->PS(ntime,:,:)      
hyam = fuv->hyam
hybm = fuv->hybm
interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
extrap = True  ; is extrapolation desired if data is outside the range of PS
P0mb = 1000    ; ps in Pa, but this argument must be in mb
        
temp=vinth2p(fuv->VOR(ntime,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
var=temp(0,:,:)
temp=vinth2p(fuv->DIV(ntime,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
var2=temp(0,:,:)
print ("var1 min/max = "+min(var(:,:))+", "+max(var(:,:)) )
print ("var2 min/max = "+min(var2(:,:))+", "+max(var2(:,:)) )




;
; if  you read in U_DYN and V_DYN instead of VOR,DIV, convert to VOR,DIV:
; using vector spherical harmonics  (as expensive as two transforms)
;
; one can compute spectra directly from the vector transform (see code below), but we
; split into vorticity and divergence in order to plot each those components
; seperately.  
;
convert_uv_to_vordiv=False
if (convert_uv_to_vordiv) then
   print("converting U,V to vor,div")
   if (latlon_gridtype .eq. 1 ) then
      vr = uv2vrG_Wrap (var,var2)
      dv = uv2dvG_Wrap (var,var2)   
   else 
      if (latlon_gridtype .eq. 2 ) then
         vr = uv2vrF_Wrap (var,var2)
         dv = uv2dvF_Wrap (var,var2)   
      else
         if (latlon_gridtype .eq. 3 ) then
            ; interpolate from fixed-offset to fixed grid
            fo2fshv_Wrap(var,var2,unew,vnew)
            vr = uv2vrF_Wrap (unew,vnew)
            dv = uv2dvF_Wrap (unew,vnew)
            latlon_gridtype=2
         else
            print("Error: unsupported grid for spherical harmonic transform")
	    exit
         end if
      end if
   end if
else
   vr=var
   dv=var2
end if


; compute spectra from vor/div
if (latlon_gridtype .eq. 1 ) then
   ab_vr = shagC (vr(:,:))       ; GAUSS grid
   ab_dv = shagC (dv(:,:))       ; GAUSS grid
else 
   if (latlon_gridtype .eq. 2 ) then
      ab_vr = shaeC (vr(:,:))       ; equally spaced
      ab_dv = shaeC (dv(:,:))       ; equally spaced
   else
      print("Error: unsupported grid for spherical harmonic transform")
      exit
   end if
end if



;
; compute energy in vort and div components
;
tmp=dimsizes(ab_dv)
nmax = tmp(1)
mmax = tmp(2) 

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
   print ("computing vector spherical harmonic transform...")
   if (latlon_gridtype .eq. 1 ) then
      ab = vhagC(var(:,:),var2(:,:));   ; GAuss grid
   else
      if (latlon_gridtype .eq. 2) then
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
end if

      

res@trYMinF = 1e-7      
res@trYMaxF = 1e3   
res@tiMainString     = "Climatology"
res@xyLineColors = (/"blue","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x,spc,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x,spc_dv,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x,spc_vr,res) ; create plot

xref = (/10.,300./)
spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
frame(wks)

end

