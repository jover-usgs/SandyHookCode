###############################################
###############################################

#3D change analysis program v2.0 edited 9/12/2018

#Sandy Hook Cooperative Research Programs
#Psuty, Schmelz, Greenberg, Liu, Ames, Habeck

###############################################
###############################################

# ################### USER EDITS <START> ######################################################### #
# Declare Site and path to working directory
#site 				= "JF"
#path 				= r'M:/11_Papers_and_Meetings/2018/ESCO_FortTilden/01_Riis_Tilden/'
site 				= "KM"
path 				= r'M:/01_NPS/01_Gateway/01_SHU/05_Projects/KingmanMills/KM_2012-2019/'


# Name of compartments and baseline to be used in analysis
compname 			= "Compartments_KM_v3.shp"
baselname			= "KM_Baseline.shp"

# Master CSV containing values for all site-specific variables
csvname 			= "3D_Analysis_Variables.csv"

# Multiplier of volumetric output for creation of vectors
vectormultiplier = .1

# Analysis Type (Must run total before running Dune or Tidal)
AnalysisType = "Tidal"			# Type of anlaysis (Total, Dune, Tidal) user determined

# Option to Skip TIN and Merged Raster Creation (saves part of the processing time)
# If these files already exist and no changes were made to them, set this parameter to 1 to skip their creation step
TinRASCreation		= 0			# 1 = SKIP; 0 = DO NOT SKIP
# ################### USER EDITS  <END>  ##########################################################
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
from decimal import Decimal
import numpy
import csv
import pandas

from numpy import matrix
from numpy import linalg
from numpy import genfromtxt
############################################################################################ #
# Declare Functions
def ElapsedTime(start,end,statement):
	elapsed 		= end - start
	if elapsed >= 60:
		elapsed 	= round(Decimal((elapsed)/60),2)
		arcpy.AddMessage(statement + " in " + str(elapsed) + " minute(s).")
	else:
		elapsed 	= round(Decimal(elapsed),2)
		arcpy.AddMessage(statement + " in " + str(elapsed) + " seconds.")
	return elapsed
# The following global variables are needed for this function to run:
# loctab, locdif, data, compcount, difs
def makePivotTable(type,AnalysisType):
	print("Making " + type + " Pivot Table...")
	# Set table names and paths
	voltab = loctab + voltabname
	pivtab = loctab + type + "_" + AnalysisType + '_Pivot'
	pivtabname = type + "_" + AnalysisType + '_Pivot'

	# Create table and add fields
	arcpy.CreateTable_management(loctab, pivtabname)
	arcpy.AddField_management(pivtab, "COMP", "SHORT")
	arcpy.DeleteField_management (pivtab, "FIELD1")
	arcpy.DeleteField_management (pivtab, "OBJECTID")
	
	# Insert Compartment Numbers into the Pivot Table
	iCur = arcpy.da.InsertCursor(pivtab,"COMP")
	for y in xrange(1,compcount + 1):
		iCur.insertRow([y])
	del iCur
	
	# Make duration fields
	for n in xrange (0, len(data)):
		difname = str(difs[n])
		dif_F_name = 'D' + difname[-13:]
		arcpy.AddField_management(pivtab, dif_F_name, "Double")
		n = n + 1
	
	for y in xrange (1, compcount + 1):
		compnum = str(y)
		# Insert Volumes for each comparison duration by Compartment Number
		for n in xrange (0, len(data)):
			# Set duration strings
			dif = locdif + difs[n]
			difname = str(difs[n])
			dif_F_name = 'D' + difname[-13:]
			
			# Set SQL clauses for compartment number and duration
			compid = str('"COMP" = ' + str(y))
			difcompid = str('"COMP" = ' + str(y)) + " AND " + str('"DURATION" = ' + "'" + difname + "'")
			
			# Get volume for compartment/duration
			with arcpy.da.SearchCursor(voltab,type,difcompid) as sCur:
				for row in sCur:
					comptimvol = row[0]
			
			#comptimvol = [row[0] for row in arcpy.da.SearchCursor(voltab,type,difcompid)][0]
		
			# Set corresponding volume entry in pivot table
			with arcpy.da.UpdateCursor(pivtab,dif_F_name,compid) as uCur:
				for row in uCur:
					row[0] = comptimvol
					uCur.updateRow(row)
			
		sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
		
	arcpy.TableToTable_conversion (pivtab, loctab, pivtabname)
	arcpy.TableToExcel_conversion (pivtab, pivtab + '.xls') 	
	print '.'
	return pivtab
###############################################

	elapsed = end - start
	if elapsed >= 60:
		elapsed = round(Decimal((elapsed)/60),2)
		arcpy.AddMessage(statement + " in " + str(elapsed) + " minute(s).")
	else:
		elapsed = round(Decimal(elapsed),2)
		arcpy.AddMessage(statement + " in " + str(elapsed) + " seconds.")
	return elapsed
# The following variables are global variables that must exist for this function to work:
# loctab, locdif, data, difs, startx, starty, baselinedistance, baselineradian
# vectormultiplier, vectorangle1, vectorangle2, vectorradian
# baseline variable indicates whether the baseline is "CURVED" or "STRAIGHT"
def DrawVolumeVectors(type,baseline, AnalysisType):
	print("Now drawing " + type + " vectors")
	voltab = loctab + voltabname
	deptabp = loctab + type + "_Pivot"
	
	# Read in CSV if baseline is curved
	if baseline == "CURVED":
		vecstab = loctab + "Vec_ori" + ".csv"
		with open(vecstab, 'rb') as csvfile:
			reader2 = csv.reader(csvfile)
			vecdistsets = list(reader2)
	else:
		baselineradian = math.atan2(deltay,deltax)
	
	# Loop for each comparison 
	for n in xrange (0, len(data)):
		
		# Set duration names
		dif = locdif + difs[n]
		difname = str(difs[n])
		dif_F_name = 'D' + difname[-13:]
		
		# Set output vector name
		vectors = "Vectors_" + type + "_"  + AnalysisType + "_" + difname
		dunedifout = voltab
		
		# Set initial Start XY Point for vector drawing
		feature_info = []
		
		# Get StartXY, length, and baselineradian from total vectors if AnalysisType is Dune
		if AnalysisType == "Dune":
			
			# Read in existing total vectors and calculate End XY
			baselname = vec + 'Vectors_' + type + "_Total_" + difname
			
			# Make starting coordinates
			vecstartx=[]
			vecstarty=[]
			vecendx=[]
			vecendy=[]
			length=[]
			dunebaselineradian=[]
			
			# Calculate Start and End XY fields in baseline, and get values
			with arcpy.da.SearchCursor(baselname,["SHAPE@","SHAPE@LENGTH","Comp"]) as sCur:
				for row in sCur:
					if row[1] == None:
						vecstartx.append(0)
						vecstarty.append(0)
						vecendx.append(0)
						vecendy.append(0)
						dunebaselineradian.append(0)
						length.append(0)
					else:
						StartX = row[0].firstPoint.X
						StartY = row[0].firstPoint.Y
						EndX = row[0].lastPoint.X
						EndY = row[0].lastPoint.Y
						DuneRadian = math.atan2(EndY - StartY, EndX - StartX)
					
						# Append variables 
						vecstartx.append(StartX)
						vecstarty.append(StartY)
						vecendx.append(EndX)
						vecendy.append(EndY)
						dunebaselineradian.append(DuneRadian)
						length.append(row[1])
		
		# Set starting XY point based on baseline type
		if baseline == "CURVED":
			q = float(vecdistsets[0][0][0:12])
			w = float(vecdistsets[0][1][0:13])
		else:
			q = startx + (baselinedistance * 0.5) * (math.cos(baselineradian))
			w = starty + (baselinedistance * 0.5) * (math.sin(baselineradian))
		
		# Get vector coordinates for each compartment based on volumetric change
		for y in range(1, compcount + 1):
			
			# Get baseline radian if vector is curved
			if baseline == "CURVED" and AnalysisType != "Dune":
				baselineradian = float(vecdistsets[y-1][2][0:12])
				
			# Calculate vector drawing angle (separate method for Dune Analysis)
			if AnalysisType == "Dune":
				if vectordirection == "left":
					vectorangle1 = dunebaselineradian[y-1]
					vectorangle2 = dunebaselineradian[y-1]
				elif vectordirection == "right":
					vectorangle1 = dunebaselineradian[y-1]
					vectorangle2 = dunebaselineradian[y-1]
			else:
				if vectordirection == "left":
					vectorangle1 = baselineradian + math.pi/2
					vectorangle2 = baselineradian - math.pi/2
				elif vectordirection == "right":
					vectorangle1 = baselineradian - math.pi/2
					vectorangle2 = baselineradian + math.pi/2
				
			# Set Starting XY Coordinates based on baseline and analysis type
			if AnalysisType == "Dune":
				q = vecstartx[y-1] + length[y-1] * math.cos(vectorangle1) #northing start
				w = vecstarty[y-1] + length[y-1] * math.sin(vectorangle1) #easting start
			else:
				if baseline == "CURVED":
					q = float(vecdistsets[y-1][0][0:12])
					w = float(vecdistsets[y-1][1][0:13])
			
			# Set SQL statements to retrieve volume value
			compid = str('"COMP" = ' + str(y))
			difcompid = str('"COMP" = ' + str(y)) + " AND " + str('"DURATION" = ' + "'" + difname + "'")
			
			# Use cursor to retrieve volume value
			with arcpy.da.SearchCursor(voltab,type,difcompid) as sCur:
				for row in sCur:
					volumedif = float(row[0])
				
					# Set volumedif to negative for erosion vectors
					if type == "Erosion":
						volumedif = volumedif * -1
				
			# Calculate Vector Distance and Radian to determine drawing coordinates
			if volumedif >= 0:
				vectordistance = vectormultiplier * volumedif * -1
				vectorradian = vectorangle2
			else:
				vectordistance = vectormultiplier * volumedif
				vectorradian = vectorangle1
		
			# Get vector endpoint values for drawing
			e = q + vectordistance * math.cos(vectorradian)
			r = w + vectordistance * math.sin(vectorradian)
			
			# Record and append feature coordinates to list
			addFeature = [[q, w], [e, r]]
			feature_info.append(addFeature)
			
			# Add set distance to prepare for drawing next vector
			if baseline == "CURVED":
				if y >= compcount:
					sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
					continue
				else:
					q = float(vecdistsets[y][0][0:12])
					w = float(vecdistsets[y][1][0:13])
			else:
				q = q + baselinedistance * (math.cos(baselineradian))
				w = w + baselinedistance * (math.sin(baselineradian))
			
			# Output processing messages and add to counter
			sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
		print "."
		
		# Create and add vectors
		features = []
		for feature in feature_info:
			features.append(
				arcpy.Polyline(
					arcpy.Array([arcpy.Point(*coords) for coords in feature]),sr))
		
		# Copy features to output vector and add fields
		arcpy.CopyFeatures_management(features, vectors)
		arcpy.AddField_management(vectors, "Comp", "Integer")
		arcpy.AddField_management(vectors, type, "Double")
		
		# Fill in values
		n = 1
		for y in range(1, compcount + 1):
		
			# Set SQL clauses to retrieve volume value
			zone = str(y)
			compid = str('"COMP" = ' + zone)
			objid = str("OBJECTID = " + zone)
			difcompid = str('"COMP" = ' + zone + "AND" + str('"DURATION" = ' + "'" + difname + "'"))
			
			# Get volume value
			with arcpy.da.SearchCursor(voltab,type,difcompid) as sCur:
				for row in sCur:
					volumedif = float(row[0])
				# Multiply volumetric change by -1 
				if type == "Erosion":
					volumedif = volumedif * -1
		
			# Make make new featurelayer and write in value
			arcpy.MakeFeatureLayer_management(vectors,"lyr1",objid)
			arcpy.CalculateField_management("lyr1", "Comp", zone, "PYTHON_9.3")
			arcpy.CalculateField_management("lyr1", type, volumedif, "PYTHON_9.3")
			sys.stdout.write("\rCalculated %i of %i vector attributes for duration %s" % (y, compcount, difname))
			y = y + 1
		print "."
	return vectors
############################################################################################ #
print("Now running 3D Analysis for " + site + "...\n")
print (time.strftime("%H:%M:%S"))
start = time.time()
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")
env.overwriteOutput = True
arcpy.overwriteOutput = True
sr = str("""PROJCS['NAD_1983_UTM_Zone_18N',
	GEOGCS['GCS_North_American_1983',
	DATUM['D_North_American_1983',
	SPHEROID['GRS_1980',6378137.0,298.257222101]],
	PRIMEM['Greenwich',0.0],
	UNIT['Degree',0.0174532925199433]],
	PROJECTION['Transverse_Mercator'],
	PARAMETER['False_Easting',500000.0],
	PARAMETER['False_Northing',0.0],
	PARAMETER['Central_Meridian',-75.0],
	PARAMETER['Scale_Factor',0.9996],
	PARAMETER['Latitude_Of_Origin',0.0],
	UNIT['Meter',1.0]]""")
	
# Read in csv
mcsv 				= r'M:/17_Software_Manuals/3D_Analysis/3D_Analysis_Scripts/' + csvname
df 					= pandas.read_csv(mcsv,header=0,index_col=1)
#print df

# Get variables from dataframe
sitename 			= df.loc[site,"site_name"]				# Full name of site, no spaces
sitenum				= df.loc[site,"site_number"]			# Number of site used for directory naming
Total_baseras 		= df.loc[site,"rasname_Total"]			# Site base raster name
Dune_baseras		= df.loc[site,"rasname_Dune"]			# Dune threshold raster
Tidal_baseras		= df.loc[site,"rasname_Tidal"]			# Tidal flat threshold raster
tin_length			= int(df.loc[site,"tin_length"])		# For creating tins (must be greater than largest transect spacing)
startdate			= str(df.loc[site,"date_start"])		# Start date of analysis
enddate				= str(df.loc[site,"date_end"])			# End date of analysis
BaselineType		= str(df.loc[site,"baselinetype"])		# Type of baseline
vectordirection		= str(df.loc[site,"vectordirection"])	# Direction of vectors
baselinedistance 	= int(df.loc[site,"baselinedistance"])
# Set necessary strings for GDB and working directory
gdb					= site + "_wkspc.gdb/"
if sitenum < 9:
	sitenum = "0" + str(sitenum)
else:
	sitenum = str(sitenum)
print("Setting workspace to the following directory: " + path) #+ sitenum + "_" + sitename)
wkspc				= path # + sitenum + "_" + sitename + "/"

# Set paths and names
locvol 				= wkspc + '20_Volumes/'				# Total volume tables
locgdb 				= wkspc + gdb						# Geodatabase location
locshp				= wkspc + '03_Shapefiles/'			# Shapefiles and survey data
loctin 				= wkspc + '04_TINs/'				# TINs output location
locras 				= wkspc + '09_Rasters/'				# Rasters output location
locdif 				= wkspc + '10_Differences/'			# Difference rasters for comparisons
loctab 				= wkspc + '11_Tables/'				# Volume tables (differences)
loccomp 			= wkspc + '13_Compartments/'		# Compartment shape files
locvec 				= wkspc + '14_Vectors/'				# Change vector shape files
loclyr 				= wkspc + '15_Layers/'				# Layer files for symbology

# Additional paths based on analysis type
if AnalysisType =="Dune":
	locdif = wkspc + '/10_Differences/Dune/'
	locvec = wkspc + '/14_Vectors/Dune/'
	rastype = "d_" # prefix for rasters ('m_' for merged, 'd_' for dune, 't_' for tidal)
	vec = wkspc + '14_Vectors/Vectors.gdb/'	# Location of total vectors, whose endpoints will be used to draw dune vectors
	baseras = Dune_baseras
elif AnalysisType =="Tidal":
	locdif = wkspc + '/10_Differences/Tidal/'
	locvec = wkspc + '/14_Vectors/Tidal/'
	rastype = "t_" # prefix for rasters ('m_' for merged, 'd_' for dune, 't_' for tidal)
	baseras = Tidal_baseras
else:
	rastype = "m_" # prefix for rasters ('m_' for merged, 'd_' for dune, 't_' for tidal)
	baseras = Total_baseras
voltabname = "Volumes_" + AnalysisType	# name of volume table (match to AnalysisType)					

integ=str(random.randint(1,12000000))

###############################################
# Add and number compartment fields
in_dem = wkspc + '/12_Infill/' + baseras
baseline1 = locvec + baselname
compartments = loccomp + compname
#arcpy.AddField_management(compartments, 'COMP', 'DOUBLE')
#arcpy.CalculateField_management(compartments, 'COMP', '!FID! +1', "PYTHON_9.3")

print "Directory:"
backslash = '\\'
print str(str(os.getcwd()).replace(backslash,"/") + "/")
print" "
# Get the total number of compartments
compcount = int(arcpy.GetCount_management(loccomp + compname).getOutput(0))
# Get name of all shape files to be used in comparisons
gpsdata = glob.glob(locshp + "*.shp")

# Get all survey dates
n=1
data = []
for n in xrange (0,len(gpsdata)):
	data.append(gpsdata[n][-12:-4])
print "Shapefiles:"
print data
print " "

# Get all Shapefile names
data2 = []
for n in xrange (0,len(gpsdata)):
	data2.append(gpsdata[n][-18:-4])
print "Shapefiles:"
print data2
print " "

# Get titles for all the analyses
months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
titles = []
for n in xrange (0, len(data)):
	# If on the last date, do net comparison. Otherwise, compare consecutive dates.
	if n == len(data)-1:
		shoname = startdate
		shoname2 = enddate
	else:
		shoname = data[n]
		shoname2 = data[n+1]
	shomm = shoname[-4:-2]
	shomm2 = shoname2[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	title = str(months[int(shomm)-1] + " " + shoyyyy + " to " + months[int(shomm2)-1] + " " + shoyyyy2)
	titles.append(title)

# Get labels (Month and Year)
labels = []
for n in xrange (0, len(data)):
	shoname = data[n]
	shomm = shoname[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	label = str(months[int(shomm)-1] + " " + shoyyyy)
	labels.append(label)
print "Surveys:"
print list(labels)
print " "
print "Analyses:"
print list(titles)
print " "
	
#set environment
arcpy.MakeRasterLayer_management(in_dem, "extlyr")
env.extent = "extlyr"
print "Running " + AnalysisType + " Analysis"

# Create TINs and merged rasters if the TinRASCreation is set to "0" 
if TinRASCreation != 1:
	if AnalysisType == "Total":
		#Create TINs
		for n in xrange (1, len(data) + 1):
			shp = locshp + data2[n-1] + ".shp"
			#breakl = locshp + data2[n-1] + "_Cliff.shp"
			TIN = loctin + data[n-1]
			RAS = locras + data[n-1]
			arcpy.CreateTin_3d(TIN, sr, [[shp, "Z", "masspoints"]],"CONSTRAINED_DELAUNAY")
			arcpy.DelineateTinDataArea_3d (TIN, tin_length)
			print data[n-1] + ' TIN created'

		#Tins to rasters
		for n in xrange (1, len(data) + 1):
			TIN = loctin + data[n-1]
			RAS = locras + data[n-1]
			arcpy.TinRaster_3d(TIN, RAS, "FLOAT", "LINEAR", "CELLSIZE 1.0") #1 m cell size 
			print data[n-1] + ' raster created'
	else:
		print "Using Rasters created from previous 'Total Analysis' Run"

	#Merge rasters
	for n in xrange (1, len(data) + 1):
		if AnalysisType == "Tidal":
			RAS = locras + 'm_' + data[n-1]
		else:
			RAS = locras + data[n-1]
		mRAS = locras + rastype + data[n-1]
		temp = "in_memory/t1" + integ
		tempn = "t1" + integ
		mRASname = rastype + data[n-1]
		if AnalysisType == "Tidal":
			arcpy.MosaicToNewRaster_management([RAS,'extlyr'], "in_memory", tempn, sr, "32_BIT_FLOAT", "1.0", "1", "MINIMUM","FIRST") #1 m cell size
			print "Taking minimum merged elevation"
		else:
			arcpy.MosaicToNewRaster_management([RAS,'extlyr'], "in_memory", tempn, sr, "32_BIT_FLOAT", "1.0", "1", "MAXIMUM","FIRST") #1 m cell size
			print "Taking maximum merged elevation"
		arcpy.gp.ExtractByMask_sa(temp,compartments,mRAS)
		arcpy.Delete_management("in_memory")
		print data[n-1] + ' data raster merged with base'

# Delete existing difference and vector folders
if os.path.exists(locdif):
	shutil.rmtree(locdif)
os.mkdir(locdif)
if os.path.exists(locvec):
	shutil.rmtree(locvec)
shutil.rmtree(locdif)
os.makedirs(locdif)
os.mkdir(locvec)

arcpy.CopyFeatures_management(locgdb + baselname[0:11], locvec + baselname)
#Subtract rasters
difs = []												# Make difs list for looping
for n in xrange (0, len(data)):
	##If on the last date, do net comparison. Otherwise, compare consecutive dates.
	if n == len(data)-1:
		startRAS = startdate							# Get YYYYMMDD
		startRASm = locras + rastype + startRAS			# Set location of merged raster
		endRAS = enddate
		endRASm = locras + rastype + endRAS
	else:
		startRAS = data[n]
		startRASm = locras + rastype + startRAS
		endRAS = data[n+1]
		endRASm = locras + rastype + endRAS
	
	RASdifname = endRAS[:6] + "_" + startRAS[:6]		# Set string for difference raster
	RASdif = locdif + RASdifname						# Set location for difference raster
	difs.append(RASdifname)
	arcpy.Minus_3d(endRASm,startRASm,RASdif)			# Output difference
	print RASdifname + ' difference raster created'

#Create data table
print difs
voltab = loctab + voltabname
arcpy.CreateTable_management(loctab, voltabname)
arcpy.AddField_management(voltab, "COMP", "SHORT")
arcpy.AddField_management(voltab, "DURATION", "TEXT")
arcpy.AddField_management(voltab, "Deposition", "Double")
arcpy.AddField_management(voltab, "Erosion", "Double")
arcpy.AddField_management(voltab, "Volume", "Double")
arcpy.DeleteField_management (voltab, "FIELD1")
arcpy.DeleteField_management (voltab, "OBJECTID")

#Calculate differences
for n in xrange (0, len(data)):
	dif = locdif + difs[n]
	difname = str(difs[n])
	for y in xrange (1, compcount + 1):
		rowid = str(((n-2)*compcount)+y)								# Get COMP ID and set clause
		comp = loccomp + 'Compartments_KM_v3.shp'
		compnum = str(y)
		compid  = str('"COMP" = ' + str(y))
		temp1 = "in_memory/t2" + integ
		arcpy.MakeFeatureLayer_management(comp, "lyr",compid)			# Create Mask for isolating compartment
		temp1 = ExtractByMask(dif,"lyr")								# Get Raster for comp				
		result = arcpy.SurfaceVolume_3d(temp1,"#","ABOVE","0","1","0")	# Get Deposition Volume for compartment

		depvolp = str(arcpy.GetMessages())
		volpos = depvolp.find('Volume=')
		depvol1 = str(depvolp[volpos+7:volpos+17])						# Get volume value from arcpy message
		try:
			depvol = float(depvol1)										# Convert to float
		except:
			depvol = 0
		result = arcpy.SurfaceVolume_3d(temp1,"#","BELOW","0","1","0")	# Get Erosion Volume for compartment
		arcpy.GetMessages()
		erovolp = str(arcpy.GetMessages())
		volpos = erovolp.find('Volume=')
		erovol1 = str(erovolp[volpos+7:volpos+17])
		try:
			erovol = float(erovol1)
		except:
			erovol = 0
		totvol = depvol - erovol										# Get Net Volumetric Change for compartment
		fields = ["Rowid", "COMP","DURATION", "Deposition", "Erosion", "Volume"]
		cursor = arcpy.da.InsertCursor(voltab, fields)					# Use InsertCursor to insert data into Volume Table
		cursor.insertRow([rowid, compnum, difname, depvol, erovol, totvol])
		del(cursor)
		sys.stdout.write("\rProcessed %i of %i for duration %s" % (y, compcount, difname))
		arcpy.Delete_management("in_memory")
	
	print '.'	

#Make Pivot Tables
analysisComponents = ["Volume","Deposition","Erosion"]
for component in analysisComponents:
	makePivotTable(component,AnalysisType)
print (time.strftime("%H:%M:%S"))

#Create Vectors
arcpy.CreateFileGDB_management(locvec, 'Vectors')
env.workspace = locvec + 'Vectors.gdb'

# Copy baseline features over and add fields
baseline = 'Baseline'
arcpy.CopyFeatures_management(baseline1, baseline)
env.extent = "MAXOF"
arcpy.AddField_management(baseline, 'StartX', 'DOUBLE')
arcpy.AddField_management(baseline, 'StartY', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndX', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndY', 'DOUBLE')

# Calculate Start and End XY fields in baseline, and get values
with arcpy.da.UpdateCursor(baseline,["SHAPE@","StartX","StartY","EndX","EndY"]) as uCur:
	for row in uCur:
		startx = row[1] = row[0].firstPoint.X
		starty = row[2] = row[0].firstPoint.Y
		endx = row[3] = row[0].lastPoint.X
		endy = row[4] = row[0].lastPoint.Y
		uCur.updateRow(row)

# Calculate the angle of the baseline
deltax = endx - startx
deltay = endy - starty

DrawVolumeVectors("Deposition",BaselineType,AnalysisType) 
print "Deposition Vectors Completed \n"
DrawVolumeVectors("Erosion",BaselineType,AnalysisType)
print "Erosion Vectors Completed \n"
print "Vectors drawn."

end = time.time()
statement = AnalysisType + " Analysis has run to completion "
ElapsedTime(start,end,statement)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #