import numpy, os, sys, getopt
#
#  interface to various NGL routines
#  collected here so MPL code can run in a conda env with no NGL support
#
#
import importlib.util
if importlib.util.find_spec("Ngl") is not None: 
    import Ngl
from netCDF4 import Dataset


def ngl_open(wks_type,outname):
    wks = Ngl.open_wks(wks_type,outname)
    return wks

def ngl_read_colormap(cmap):
    cmap2=Ngl.read_colormap_file(cmap)
    return cmap2

def ngl_end():
    Ngl.end()

def ngl_vinth2p(dataf2,hyam,hybm,plev,PS2,v_interp,P0mb,i_opt,extrap):
    x=Ngl.vinth2p(dataf2,hyam,hybm,plev,PS2,v_interp,P0mb,i_opt,extrap)
    return x



def extract_level(dataf,klev,plev,PS,hyam,hybm,kidx=0):

    if kidx==0:  # intepolate first index:
        if plev == None:
            data2d=dataf[klev,...]
        else:
            print("Interpolating 1st dimension")
            # vertical interpolation
            v_interp = 2   # type of interpolation: 1 = linear, 2 = log, 3 = loglog
            extrap = True  # is extrapolation desired if data is outside the range of PS
            P0mb = 1000    # ps in Pa, but this argument must be in mb
            
            if len(dataf.shape)==2:   # lev,ncol
                dataf2=numpy.expand_dims(dataf,axis=2)  # lev,ncol,1
                PS2=numpy.expand_dims(PS,axis=1)        # ncol,1
                data2d=numpy.squeeze(ngl_vinth2p(dataf2,hyam,hybm,plev,PS2,v_interp,P0mb,1,extrap))
            elif len(dataf.shape)==3:  # lev,lat,lon
                data2d=numpy.squeeze(ngl_vinth2p(dataf,hyam,hybm,plev,PS,v_interp,P0mb,1,extrap))
            else:
                print("ERROR: extract_level: dataf() needs to be 2 or 3 dimensiosn")
                
        return data2d

    if kidx==len(dataf.shape)-1:
        print("need to interpolating last dimension!")
        if plev == None:
            data2d=dataf[...,klev]
        else:
            # vertical interpolation
            print("vertical interpolation last dimension, not coded")
            sys.exit(1)
        return data2d

    print("Error: level dimension was not first or last")
    sys.exit(1)


def ngl_plot(wks,data2d,lon,lat,title,longname,units,
             projection,clev,cmap,scrip_file,se_file,contour_opt,coutlines,data2d_2=numpy.array([])):

    se_num=0
    if os.path.isfile(se_file):
        print("adding to plot:  SE grid from:",se_file)
        infile = Dataset(se_file,"r")
        se_coord  = infile.variables["coord"][:,:]
        se_connect  = infile.variables["connect1"][:,:]
        se_num=se_connect.shape[0]
        print("number of elements =",se_num)
        #print("se_coord shape: ", se_coord.shape)
        #print("se_connectshape: ", se_connect.shape)
        # j=se_connect(i,0:3) is the index of the 4 corners of cell i
        # se_coord(0:2,j) are the x,y,z coords of vertex j
        se_lat = numpy.arcsin(se_coord[2,:])*180/numpy.pi
        se_lon = numpy.arctan2(se_coord[1,:],se_coord[0,:])*180/numpy.pi
        

    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Dataset(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        if clon.shape[0] == len(lat):
            cellbounds=True
    
    res = Ngl.Resources()

    #if cmap!=None:
    res.cnFillPalette   = cmap

    nlevels=-1  # not specified
    if len(clev)==1:
        res.cnLevelSelectionMode  = "AutomaticLevels"
        res.cnMaxLevelCount = clev[0]
        nlevels=clev[0]
    elif len(clev)==4:
        # log spacing
        cbase=clev[3]
        c0=numpy.log(clev[0])/numpy.log(cbase)
        c1=numpy.log(clev[1])/numpy.log(cbase)
        cinc=clev[2]

        #nlevels=(clev[1]-clev[0])/clev[2]
        nlevels=(c1-c0)/cinc
        clevs=[ c0+i*cinc  for i in range(1+round(nlevels))]
        clevs=numpy.power(cbase,clevs)
        print("nlevels=",nlevels," lev spacing ratio= ",cbase,"^",cinc)
        print("log clevs=",clevs)

        res.cnLevelSelectionMode = "ExplicitLevels" 
        res.cnLevels=clevs
    else:
        res.cnLevelSelectionMode  = "ManualLevels"
        res.cnMinLevelValF=clev[0]
        res.cnMaxLevelValF=clev[1]
        res.cnLevelSpacingF=clev[2]
        nlevels=(clev[1]-clev[0])/clev[2]




# defaults. some projection options might change:
    res.cnFillOn              = True           # Turn on contour fill.
    if contour_opt=='lo':
        res.cnFillOn              = False      # Turn on contour fill off.

    res.cnLinesOn             = False          # Turn off contour lines
    if contour_opt=='lo' or contour_opt=='la':
        res.cnLinesOn             = True       # Turn on contour lines

    res.cnLineLabelsOn        = False          # Turn off line labels.
    res.cnInfoLabelOn         = False          # Turn off info label.

    res.mpOutlineOn          = True
    if coutlines==0:     res.mpOutlineOn          = False

    res.mpFillOn             = False
    res.mpGridAndLimbOn      = False    # dont draw grid lines
    #res.mpShapeMode          = "FreeAspect"
    



    if projection == "latlon" :
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
        #    res.mpMinLatF = -90.
        #    res.mpMaxLatF = 90.
        #    res.mpMinLonF = -180.
        #    res.mpMaxLonF = 180.
    elif projection == "latlon-nc1":
        # no continenents, turn on contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
    elif projection == "latlon-nc2":
        # no continenents, turn off contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
    elif projection == "US1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -30
        res.mpMaxLatF = 75
        res.mpMinLonF = -180
        res.mpMaxLonF =  0
    elif projection == "europe":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 20.
        res.mpMaxLatF = 75.
        res.mpMinLonF = -40.
        res.mpMaxLonF =  40.
    elif projection == "andes":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -40.
        res.mpMaxLatF = 15.
        res.mpMinLonF = -100.
        res.mpMaxLonF =  -40.
    elif projection == "pacific1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -20.
        res.mpMaxLatF = 60.
        res.mpMinLonF = -180.
        res.mpMaxLonF =  -80.
    elif projection == "andes2":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -17.
        res.mpMaxLatF =  3.
        res.mpMinLonF = -85.
        res.mpMaxLonF =  -65.
    elif projection == "himalaya":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 0.
        res.mpMaxLatF = 60.
        res.mpMinLonF = 50.
        res.mpMaxLonF = 110.
    elif projection == "oro":
        res.mpProjection      = "Orthographic"   
        res.mpCenterLatF      =  45.
        res.mpCenterLonF      = -45.
    elif projection == "oro-debug1":
        res.mpProjection      = "Orthographic"   
        res.mpLimitMode="Angles"
        hdeg=10
        res.mpLeftAngleF=hdeg
        res.mpRightAngleF=hdeg
        res.mpBottomAngleF=hdeg
        res.mpTopAngleF=hdeg
    elif projection == "debug1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 30.
        res.mpMaxLatF = 50.
        res.mpMinLonF = 85.
        res.mpMaxLonF = 105.
    elif projection == "debug2":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -10.
        res.mpMaxLatF = 75.
        res.mpMinLonF = 45.
        res.mpMaxLonF = 175.
    elif projection == "debug3":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 10.
        res.mpMaxLatF = 50.
        res.mpMinLonF = 60.
        res.mpMaxLonF = 100.
    elif projection == "baroclinic":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpCenterLonF         = 100.
        res.mpMinLatF = 25.
        res.mpMaxLatF = 75.
        res.mpMinLonF = 25.
        res.mpMaxLonF = 175.
    elif projection == "barotopo":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        #res.mpCenterLonF         = -90.
        res.mpMinLatF = 15. 
        res.mpMaxLatF = 75. 
        res.mpMinLonF = -150.
        res.mpMaxLonF = 20.
    elif projection == "barotopo-nc2":
        # no continenents, turn off contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        #res.mpCenterLonF         = -90.
        res.mpMinLatF = 15. 
        res.mpMaxLatF = 75. 
        res.mpMinLonF = -150.
        res.mpMaxLonF = 20.
    else:
        print("Bad projection argument: ",projection)
        sys.exit(3)

        
    res.nglFrame = False # Don't advance frame.
    #res.tiXAxisString = "~F25~longitude"
    #res.tiYAxisString = "~F25~latitude"
    res.nglPointTickmarksOutward = True
    

    if cellbounds:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = clon
        res.sfYCellBounds = clat
    else:
        if contour_opt=='area' or contour_opt=='la':
            res.cnFillMode            = "AreaFill"
        elif contour_opt=='raster':
            res.cnFillMode            = "RasterFill"
        else:
            res.cnFillMode            = "AreaFill"
            if (nlevels>25):
                print("more than 25 contours, switching to raserfill")
                res.cnFillMode            = "RasterFill"
        res.cnRasterSmoothingOn = True
        res.sfXArray = lon[:]
        res.sfYArray = lat[:]



    #res.sfCopyData = False
 
#    "not a valid resource in contour at this time...       
#    if wks_type == "pdf":
#        res.gsnMaximize           = True        
#        res.gsnPaperOrientation   = "portrait"

    #res.lbLabelAutoStride   = True         # Clean up labelbar labels.
    #res.lbAutoManage = True
    #res.lbLabelStride       = 10
    res.lbBoxLinesOn        = False        # Turn of labelbar box lines.
    res.lbOrientation       = "horizontal"
    


    print("Title: ",title)
    print("Longname: ",longname)
    
    res.tiMainString = title

    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))        
    if res.cnLevelSelectionMode == "ManualLevels":
        print("contour levels: manual [",res.cnMinLevelValF,",",\
              res.cnMaxLevelValF,"] spacing=",res.cnLevelSpacingF)
        print("number of contour levels:",nlevels)
    elif res.cnLevelSelectionMode=="AutomaticLevels":
        print("contour levels: auto. number of levels:",res.cnMaxLevelCount)


    if nlevels>20:
        res.lbLabelStride       = nlevels/8
        

    # for lat/lon plots, add cyclic point:
    res2=res
    if hasattr(res,"sfXArray"):
        if len(res.sfXArray)*len(res.sfYArray) == numpy.prod(data2d.shape):
            print("NGL structured data plot.  Adding cyclic point")
            data2d,lon2= Ngl.add_cyclic(data2d[:,:],res.sfXArray[:])
            res2.sfXArray=lon2
        else:
            print("NGL unstructered plot with internal triangulation")
    elif hasattr(res,"sfXCellBounds"):
        print("NGL unstructered plot with cell bounds")
    else:
        print("Error with resource coordinate data")

    if data2d_2.size == 0:
        # contour data2d
        map1 = Ngl.contour_map(wks,data2d,res2)
    else:
        # plot variable, then add contour lines from data2d_2:
        res2.nglDraw  = False
        map1 = Ngl.map(wks,res2)
        map2 = Ngl.contour(wks,data2d,res2)
        print("Adding contour line plot from data2d_2...")
        print("data2d_2 min/max=",numpy.amin(data2d_2),numpy.amax(data2d_2))

        res3=res
        res3.mpOutlineOn          = False
        res3.cnFillOn             = False         # Turn off contour fill.
        res3.cnLinesOn            = True          # Turn on contour lines
        #res3.cnLineColor          = "White"

        #res3.cnLevelSelectionMode  = "AutomaticLevels"
        #res3.cnMaxLevelCount = 5
        #res3.cnLevelSpacingF=1e5   # ignored, but set to prevent warnings

        res3.cnMinLevelValF=200
        res3.cnMaxLevelValF=300
        res3.cnMaxLevelCount = 100
        res3.cnLevelSpacingF=5  
        res3.cnLineThicknessF = 1.5

        map3 = Ngl.contour(wks,data2d_2,res3)
        Ngl.overlay(map1,map2)
        Ngl.overlay(map1,map3)
        Ngl.maximize_plot(wks,map1)    # Maximize size of plot in frame.
        Ngl.draw(map1)
        del res3

    del res2
    print("Contour done.")
        
    #-- write variable long_name and units to the plot
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.016
    Ngl.text_ndc(wks,longname,0.20,0.95,txres)
    Ngl.text_ndc(wks,units, 0.85,0.95,txres)
    del txres

    if se_num>0:
        gsres = Ngl.Resources()
        gsres.gsLineColor            = "black"
        gsres.gsLineThicknessF       = 0.05
        for x in se_connect:
            x1=se_lon[x-1]
            x2=x1[0:1]
            plon=numpy.concatenate( [x1,x2])
            x1=se_lat[x-1]
            x2=x1[0:1]
            plat=numpy.concatenate( [x1,x2])

            # se_coord[0:2,x-1] = for element corners in cartesian coordinate.
            # compute diagonal distance:
            #xcoord=se_coord[:,x-1]
            #d1=numpy.sum((xcoord[:,0]-xcoord[:,2])**2)
            #d2=numpy.sum((xcoord[:,1]-xcoord[:,3])**2)
            #dist2 = ((numpy.sqrt(d1/2)+numpy.sqrt(d2/2))/2)* 6.4e3  # dist in km
            #if (dist2> 80):    # NE30 elements are about 325km
            Ngl.polyline(wks,map1,plon,plat,gsres)
        del gsres
            
    
    Ngl.frame(wks)       # advance frame
    return map1



def ngl_vertprofile(wks,coldata,lev,xbnds,ybnds,title,units,time):

    res = Ngl.Resources()
    res.tiMainString           = title+"~C~"+"time="+str(time)
    res.tiXAxisString          = units
    res.tiYAxisString          = "hPa"

    print("vert profile title=",title)
    if len(xbnds)==2:
        res.trXMinF=xbnds[0]     
        res.trXMaxF=xbnds[1]
    if len(ybnds)==2:
        res.trYMinF=ybnds[0]     
        res.trYMaxF=ybnds[1]
     
    res.trYReverse        = True      

    #res.vpWidthF               =  0.9                   #-- viewport width
    #res.vpHeightF              =  0.6                   #-- viewport height

    #res.caXMissingV            =  -999.                 #-- indicate missing value
    #res.caYMissingV            =  -999.                 #-- indicate missing value

    #-- marker and line settings
    res.xyLineColors           =  ["blue","green","red"] #-- set line colors
    res.xyLineThicknessF       =  3.0                    #-- define line thickness
    #res.xyDashPatterns         =  [0,0,2]                #-- ( none, solid, cross )
    #res.xyMarkLineModes        =  ["Markers","Lines","Markers"] #-- marker mode for each line
    #res.xyMarkers              =  [16,0,2]               #-- marker type of each line
    #res.xyMarkerSizeF          =  0.01                   #-- default is 0.01
    #res.xyMarkerColors         =  ["blue","green","red"] #-- set marker colors
    
    #-- legend settings
    #res.xyExplicitLegendLabels = [" data"," linear"," square"]  #-- set explicit legend labels
    #res.pmLegendDisplayMode    = "Always"               #-- turn on the drawing
    #res.pmLegendOrthogonalPosF = -1.13                  #-- move the legend upwards
    #res.pmLegendParallelPosF   =  0.15                  #-- move the legend to the right
    #res.pmLegendWidthF         =  0.2                   #-- change width
    #res.pmLegendHeightF        =  0.10                  #-- change height
    #res.lgBoxMinorExtentF      =  0.16                  #-- legend lines shorter
    
    #-- draw the plot
    plot = Ngl.xy(wks,coldata[:,:].T,lev[:,:].T/100,res)
    #Ngl.frame(wks)       # advance frame


