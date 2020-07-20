import numpy, os, sys, getopt
import Ngl
import Nio

#
# setup resources for contour plots
#
# can handle 3 types of plots:
#
# 1. latlon data:  data(nlat,nlon) with lat(nlat) and lon(nlon)
# 2. unstructured:  data(ncol), lat(ncol), lon(ncol)
# 3. unstructured with cell fill:  data(ncol), lat(ncol), lon(ncol)
#        and scrip_file is a file containing the cell bounds:
#        grid_corner_lat(ncol,4), grid_corner_lon(ncol,4)
#
#

def myargs(argv):
    inputfile = ''
    scripfile = ''
    varname = ''
    name = argv[0]
    try:
        opts, args = getopt.getopt(argv[1:],"i:s:")
    except getopt.GetoptError:
        print (name,' -i <inputfile> -s <scriptfile> varname')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-s", "--sfile"):
            scripfile = arg
                
    return inputfile,scripfile,args



def map_setup(wks_type,projection,lat,lon,scrip_file):
    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Nio.open_file(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        if clon.shape[0] == len(lat):
            cellbounds=True
            
    
    res = Ngl.Resources()
    if projection == "latlon" :
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
        #    res.mpMinLatF = -90.
        #    res.mpMaxLatF = 90.
        #    res.mpMinLonF = -180.
        #    res.mpMaxLonF = 180.
    elif projection == "US1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -30
        res.mpMaxLatF = 75
        res.mpMinLonF = -180
        res.mpMaxLonF =  0
    elif projection == "andes":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -40.
        res.mpMaxLatF = 15.
        res.mpMinLonF = -100.
        res.mpMaxLonF =  -40.
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

        
    res.nglFrame = False # Don't advance frame.
    #res.tiXAxisString = "~F25~longitude"
    #res.tiYAxisString = "~F25~latitude"
    res.nglPointTickmarksOutward = True
    
    res.cnFillOn              = True           # Turn on contour fill.
    #res.cnFillPalette         = "BlueYellowRed"  # good for symmetric data
    #res.cnFillPalette         = 'MPL_viridis'
    #res.cnFillPalette         = "wgne15"
    #res.cnFillPalette         = "StepSeq25"
    res.cnFillPalette         = "WhiteBlueGreenYellowRed"
    #res.cnFillPalette         = "BlAqGrYeOrReVi200"
    res.cnLinesOn             = False          # Turn off contour lines
    res.cnLineLabelsOn        = False          # Turn off line labels.
    res.cnInfoLabelOn         = False          # Turn off info label.
    res.cnLevelSelectionMode  = "AutomaticLevels"

    if cellbounds:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = clon
        res.sfYCellBounds = clat
    else:
        res.cnFillMode            = "RasterFill"
        res.sfXArray = lon[:]
        res.sfYArray = lat[:]
    
    if wks_type == "pdf":
        res.gsnMaximize           = True        
        res.gsnPaperOrientation   = "portrait"

    res.lbLabelAutoStride   = True         # Clean up labelbar labels.
    res.lbLabelStride       = 10 
    res.lbBoxLinesOn        = False        # Turn of labelbar box lines.
    res.lbOrientation       = "horizontal"
    
    res.mpOutlineOn          = True
    res.mpFillOn             = False
    res.mpGridAndLimbOn      = False    # dont draw grid lines
    #res.mpShapeMode          = "FreeAspect"
    
    return res



def map_plot(wks,res,data2d,title):
    res.tiMainString = "Vortex test field"
    print("Title: ",res.tiMainString)

    longname=""
    units=""
    if hasattr(data2d,"long_name"):
        longname=data2d.long_name
    if hasattr(data2d,"units"):
        units=data2d.units

    if res.cnLevelSelectionMode == "ManualLevels":
        print("contour levels: manual [",res.cnMinLevelValF,",",\
              res.cnMaxLevelValF,"] spacing=",res@cnLevelSpacingF)
    else:
        print("contour levels: auto")
        res.cnMaxLevelCount = 50
        


    # for lat/lon plots, add cyclic point:
    res2=res
    if hasattr(res,"sfXArray"):
        if len(res.sfXArray)*len(res.sfYArray) == numpy.prod(data2d.shape):
            print("Detected structured data.  Adding cyclic point")
            data2d,lon2= Ngl.add_cyclic(data2d[:,:],res.sfXArray[:])
            res2.sfXArray=lon2
        else:
            print("Unstructered plot with internal triangulation")
    elif hasattr(res,"sfXCellBounds"):
        print("Unstructered plot with cell bounds")
    else:
        print("Error with resource coordinate data")
            
    # plot:
    map = Ngl.contour_map(wks,data2d,res2)
    del res2
        
    #-- write variable long_name and units to the plot
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.012
    Ngl.text_ndc(wks,longname,0.14,0.82,txres)
    Ngl.text_ndc(wks,units, 0.95,0.82,txres)
    del txres
    
    Ngl.frame(wks)       # advance frame
    return map


    

