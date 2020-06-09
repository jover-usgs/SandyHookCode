#This is a script that will update the symbology of the vector layers all at once. I don't understandy why it doesn't work. It works if you put the individual parts into the python window in ArcMap and setting mxd to "current", but no erros pop up, but the layers aren't updated. See last section for direct to ArcMap
#Jin-Si Over 1/20/2020

#Edited to work in python 3.x

import os
import sys
try:
    import archook #The module which locates arcgis
    archook.get_arcpy()
    import arcpy
except ImportError:
    print ("can't find arcpy")
	
import arcpy.mapping as mapping

#---------------------------------------------------------------------------------------------
path 				= r'M:/01_NPS/01_Gateway/01_SHU/05_Projects/KingmanMills/KM_2012-2019/'

mxd = mapping.MapDocument(path + "06_ArcMap\KM_VolumetricChange.mxd") 	#Mapping document
print("Accesing " + mxd.title)

vectorType = "Dune"  #Type of anlaysis (Total, Dune, Tidal) user determined
type = "Erosion"     #Type of vector (Erosion or Deposition) user determined
layerLoc = (path + "15_Layers/" + vectorType + "_" + type + ".lyr")
source = mapping.Layer(layerLoc)
#target = mapping.ListLayers(mxd,'Vectors_' + type + '*',df)
print (source)
#------------------------------------------------------------------------------------------

for df in mapping.ListDataFrames(mxd):
	if (df.name == 'ElevationChange'):         #This is a user controlled variable - its whatever the dataframe you put the vectors into is called
		targets = mapping.ListLayers(mxd,'Vectors_' + type + '*',df)
		for layer in targets:
			print ("Updating " + layer.name)
			test = mapping.UpdateLayer(df,layer,source,True) #make this a variable and call it in and see if its working
			print (test)
print ("Complete")		

#---------------------------------------------------------------------------------------------------------
#For use in the python window in ArcMap - this saves time if you have more than 2 comparisons, otherwise manually changing is faster. You have to change the .lyr name and the targets name when switching to tidal and total. If you've already separated te vetors into individual data frames - then do not include line 46 and make the target layer more specific ex: 'Vectors_Erosion_Dune*', this will do a search through all the dataframes you have.

#>>> import arcpy.mapping as mapping
#>>> mxd = mapping.MapDocument("CURRENT")
#>>> source = mapping.Layer(r"M:\01_NPS\01_Gateway\01_SHU\05_Projects\KingmanMills\KM_2012-2019\15_Layers\Dune_erosion.lyr")
#>>> for df in mapping.ListDataFrames(mxd):
#...     if (df.name == 'ElevationChange'):
#...         targets = mapping.ListLayers(mxd,'Vectors_Erosion*',df)
#...         for layer in targets:
#...             mapping.UpdateLayer(df,layer,source,True)
