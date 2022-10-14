import numpy, os, sys, getopt

# needed for ngl_plot
import Ngl

def ngl_vertprofile(wks,coldata,lev,vrange,title,units):

    res = Ngl.Resources()
    res.tiMainString           = title   
    res.tiXAxisString          = units
    res.tiYAxisString          = "hPa"

    if len(vrange)==2:
        res.trXMinF=vrange[0]     
        res.trXMaxF=vrange[1]
     
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


