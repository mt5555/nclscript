;----------------------------------------------------------------------
; sfElementNodes test code
;
; at one point I could specify the triangularization by chopping the
; GLL subcells into triangles.
; this used to work, for example:
;    https://trac.mcs.anl.gov/projects/parvis/wiki/NCLExamples
; 2020/7 this doesn't work, cant figure out why.
;
; In this script if you dont set res@sfElementNodes and NCL will generate its
; own triangularization using the Shewchuk "Triangle" package, which
; probably does a better job
; http://www.cs.cmu.edu/~quake/triangle.html
; With the internal triangularization, the plot works fine
;
; NCL tests such as trimesh_5.ncl also work fine:
;   https://www.ncl.ucar.edu/Applications/trimesh.shtml
;
; The plot quality is not very good - with RasterFill, it is using a
; solid fill of the dual grid to the triangle grid, creating star like
; shapes.
;
; much better option for GLL data: a scheme which can plot vertex data
; using linear shading within each subcell.  NCL doesn't support this,
; use mapplotlib instead.
; 

begin
    
; read a "latlon" file
fnames =  systemfunc("ls ne4np4_latlon*.nc");
template_latlon=addfile(fnames,"r")       ; Open netCDF files.
if (.not. ismissing(template_latlon)) then
  nd=dimsizes(template_latlon->element_corners)
  print("nd="+nd)
  ncells=nd(1)
  ntris = ncells*2
  print("computing triangularization ncells="+ncells)
  trigrid = new((/ntris,3/),"integer")
  trigrid(::2,0)=template_latlon->element_corners(0,:)
  trigrid(::2,1)=template_latlon->element_corners(1,:)
  trigrid(::2,2)=template_latlon->element_corners(2,:)
  trigrid(1::2,0)=template_latlon->element_corners(0,:)
  trigrid(1::2,1)=template_latlon->element_corners(2,:)
  trigrid(1::2,2)=template_latlon->element_corners(3,:)
  print("trigrid min/max="+min(trigrid)+" "+max(trigrid))
end if


  
  lat   = template_latlon->lat      ; ditto
  lon   = template_latlon->lon      ; ditto
;  ele   = f->ele      ; 13044 x 3 (triangles that connect the cells)
  depth = lat
  depth@lat1d = lat
  depth@lon1d = lon
  
  printVarSummary(depth)
  printVarSummary(lat)
  printVarSummary(trigrid)

  wks = gsn_open_wks("pdf","tritest")

  res                     = True
  res@gsnMaximize         = True

  res@sfYArray            = lat
  res@sfXArray            = lon
  res@sfElementNodes      = trigrid
  res@sfFirstNodeIndex    = 1          ; indexing starts at 1, not 0
  res@cnFillMode = "RasterFill"

  res@cnFillOn            = True        ; Turn on contour fill.
  res@cnLinesOn           = False       ; Turn off contour lines
  res@mpOutlineOn = False
  res@mpFillOn = False

  ; so we can add triangles
  res@gsnFrame = False     
  map = gsn_csm_contour_map(wks,depth,res)

  resp=True

  ; draw the triangles
  do i=0,ntris-1
    j0=trigrid(i,0)-1 ; convert to zero based for NCL arrays
    j1=trigrid(i,1)-1 ; convert to zero based for NCL arrays
    j2=trigrid(i,2)-1 ; convert to zero based for NCL arrays
    
    xlon=(/lon(j0),lon(j1),lon(j2),lon(j0)/)
    xlat=(/lat(j0),lat(j1),lat(j2),lat(j0)/)
    gsn_polyline(wks,map,xlon,xlat,resp)
  end do

  frame(wks)

end

