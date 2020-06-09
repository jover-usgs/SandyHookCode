# -----------------
# Bootstrap Python
# -----------------
# First step in uncertainty calculations for 3D Analysis.
# Programatically removes 1 transect, generates a DEM, and compares interpolated
# values with the original merged raster. 
#
#
# See '3D_Analysis_Full _Documentation' in 3D_Analysis Folder or '3D_Analysis_Bootstrap_Documentation'(newer) or "Workflow_Resilience_revised.docx" 
# in 17_Software_Manuals for further information on required directories and set-up.
#
import os
import sys
try:
    import archook #The module which locates arcgis
    archook.get_arcpy()
    import arcpy
except ImportError:
    print "can't find arcpy"
from arcpy.sa import *
from arcpy import env
import datetime
import re
import shutil
import time
import exceptions, traceback
import glob
import random
import numpy

from numpy import matrix
from numpy import linalg
from numpy import genfromtxt

print (time.strftime("%H:%M:%S"))

arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")
env.overwriteOutput = True
arcpy.overwriteOutput = True

#-----------------------set up file locations and user edits------------------------------
backslash = '\\'
sr = str("PROJCS['NAD_1983_UTM_Zone_18N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-75.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

wkspc = 'C:/Users/geo/Desktop/Over_JINSI/KM_2019/KM_2012-2019/'

Datenum=str(20191125)	#Only one year can be run at a time	
Site= "KM"				#Change site as needed
baseras = "ngvd29"  	#change if needed
tran_num = int(103)  	#number of transects

locshp = wkspc + '18_BS/03_Shp/'
loctin = wkspc + '18_BS/05_TINS/'+ Datenum + "/"
locras = wkspc + '18_BS/09_Rasters/'+ Datenum + "/"
locfras = wkspc + '09_Rasters/m_' + Datenum
locdif = wkspc + '18_BS/11_Diff/' + Datenum + "/" 
loccomp = wkspc + '13_Compartments/'
locdifpoints = wkspc + '18_BS/14_DiffPoints/'
locelvpoints = wkspc + '18_BS/15_ElPoints/'

in_dem = wkspc + '12_Infill/' + baseras
compname = "Compartments_KM_v3.shp"	#Changeable
compartments = loccomp + compname
integ=str(random.randint(1,12000000))
print "Directory:"
print wkspc
#print str(str(os.getcwd()).replace(backslash,"/") + "/")
print" "
#---------------------------------------------------------------------------

#get shapefiles from shapefiles folder
n=1
data = []
gpsdata = glob.glob(locshp + "*.shp")
for n in xrange (0,len(gpsdata)):
	data.append(gpsdata[n][-12:-4])
print "Shapefiles:"
print data
print " "

#set environment
arcpy.MakeRasterLayer_management(in_dem, "extlyr")
env.extent = "extlyr"

# Check if required directories exist, and make them if they don't
reqDir = [loctin,locras,locfras,locdif]
for dir in reqDir:
	if not os.path.exists(dir):
		print "Creating missing directory: " + dir
		os.makedirs(dir)
		
#-------------------Bootstrap----------------------------------
#pull one transect out of data and create a TIN&DEM
for num in range(2,tran_num):
	Transectnum = str(num)
	y=num
	tid = '\"TransectLa\" <> \'' + str(y) + '\''			#Make sure the attribute label is Transect, the newer surveys may have it as TransectLa
	#tid  = str('"TransectLa" <> ' + '\''+ str(y) + '\'') # Need to put single quotes around str(y) in SQL query
	print tid
	lyr = locshp + Site +"_3D_" + data[0]+ ".shp"
	arcpy.MakeFeatureLayer_management(locshp + Site +"_3D_" + data[0]+ ".shp", "lyr")
	arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", tid)
	
	TIN = loctin +"no_" + Transectnum #create TIN without one transect
	RAS = locras +"no_" + Transectnum #create RASTER without one transect
	arcpy.CreateTin_3d(TIN, sr, [["lyr", "z", "masspoints"]],"CONSTRAINED_DELAUNAY")
	arcpy.TinRaster_3d(TIN, RAS, "FLOAT", "LINEAR", "CELLSIZE 1.0") #1 m cell size
	
	mRAS = locras + 'm_' + "no_" + Transectnum #merge transectless DEM with base
	temp = "in_memory/t1" + integ
	tempn = "t1" + integ
	mRASname = 'm_' + "_no_" + Transectnum
	arcpy.MosaicToNewRaster_management([RAS,'extlyr'], "in_memory", tempn, sr, "32_BIT_FLOAT", "1.0", "1", "MAXIMUM","FIRST") #1 m cell size
	arcpy.gp.ExtractByMask_sa(temp,compartments,mRAS)
	arcpy.Delete_management("in_memory")
	#print "Raster without transect " + Transectnum + "created"
	
	#subtract from full dataset DEM
	RASdif = locdif + "diff_no_" + Transectnum
	outMinus = Minus(mRAS,locfras)
	outMinus.save(RASdif)
	print "Subtracted from full DEM"
	
	#pull elevation along transect at survey points from interpolated DEM
	tid2 = '\"TransectLa\" = \'' + str(y) + '\''
	#tid2  = str('"TransectLa" <> ' + '\''+ str(y) + '\'') # Need to put single quotes around str(y) in SQL query
	arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", tid2)
	transectpoints = locelvpoints + "Diff_Points_" + "no_" + Transectnum
	ExtractValuesToPoints("lyr", mRAS, transectpoints)
	print "Spot Elevation Extracted"

difdata = glob.glob(locelvpoints + "*.shp")
#print difdata
spotdif=locelvpoints + "Diff_Elev_Points_all.shp"
arcpy.Merge_management(difdata, spotdif)


