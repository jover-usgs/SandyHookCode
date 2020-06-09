# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Make sure to remove existing xlsx spreadsheet 
import pandas
import math
import os
import shutil
import openpyxl

# XXX - USER EDITED VARIABLES <START> - XXX # --------------------------------

# Target CSV
wkdir = 'C:/Users/geo/pytest/2D_Test/'
csvname   = 'GATE_JBU_BP_2D_20181003_combined.csv'
outdir = 'PointsByTransect/'
site = 'BP'
survey = 'F2018' # Survey season that will be pasted to end of xlsx name
delim = '_' # Delimiter b/w Transect and Point Label
ptsplit = 0 # set to 1 if TransectLabel nd PointLabel are already split
pointcol = 'PointLabel' # Column Label of points
thres = 0 # elevation threshold seaward points will be interpolated to
# Transect DB
bmdir = 'M:/17_Software_Manuals/2D_3DDataProcessing/'
bmcsvname = 'GATE_Transect_DB.csv'

# Options
overwritedir = 1 # Set to 1 if you want existing xlsx files to be overwritten


# Clean characters [TEMP FIX](need to figure out encoding)
# (characters to remove from the point label)
# if any unicode characters are throwing errors, add them to this list
repchar = [' ', '\x90']

# XXX - USER EDITED VARIABLES <END>   - XXX # --------------------------------

# Set directory locations
csv = wkdir + csvname
bmcsv = bmdir + bmcsvname 
loccout = wkdir + outdir

# Delete existing directory if condition is set
if overwritedir == 1:
    shutil.rmtree(loccout)
    os.makedirs(loccout)
else:
    if not os.path.exists(loccout):
        print 'Creating ' + outdir + ' in working directory.'
        os.makedirs(loccout)

# Read in csv
df = pandas.read_csv(csv, header=0)
bdf = pandas.read_csv(bmcsv, header=0)

# Get list of columns
cols = df.columns.tolist()

# Replace odd characters
for char in repchar:
    df.PointLabel = df.PointLabel.apply(lambda x: x.replace(char,''))


# Split PointLabel based on setting
if ptsplit == 1:
    df['Transect'] = df['TransectLabel']
    df['Point'] = df['PointLabel']
else:
    # Split point label into transect and point columns
    new = df[pointcol].str.split(delim,n=1,expand=True)

    # Reassign to original dataframe
    df['Transect'] = new[0]
    df['Point'] = new[1]
    
# Create 'Source' Column to indicate origin
df['Source'] = survey
        
# Get list of unique transects (filter out any ref stations or QAQC)
trans = df.Transect.unique()
transid = trans.tolist()
sitetrans = [t for t in transid if site in t]

# Loop for each transect
for t in sitetrans:
    print 'Processing ' + t
    # Select all points in that transect
    pts = df[df['Transect'] == t]
    
    # Try to find the 'BM'
    bm = pts[pts['Point'] == 'BM']
    print 'length of bm is ' + str(len(bm))
    if len(bm) < 1:
        print 'Could not find point labeled [BM]...',
        print 'Drawing from Benchmark Database'
        bm = bdf[bdf['Transect'] == t]      
    elif len(bm) > 1:
        print 'WARNING: Found more than 1 [BM]...',
        print 'Taking the first BM'
        
        # Get data from first row
        bm = bm.drop_duplicates(subset=['Point'], keep='first')
        #bm = bm.iloc[1] # turns 1st row into a series.
        #starty = bm.loc['Northing_meters'].tolist()
        #startx = bm.loc['Easting_meters'].tolist()
        #startz = bm.loc['Elevation_meters'].tolist()
    
    #if len(bm) <= 1:   
        # Get start x,y,z
    startx = bm['Easting_meters'].tolist()[0]
    starty = bm['Northing_meters'].tolist()[0]
    startz = bm['Elevation_meters'].tolist()[0]
        
        
    # Remove benchmarks and Nonetypes from point list to allow for sorting
    pts = pts.dropna(subset=['Point'])
    pts = pts[~pts.Point.str.contains('BM')]
    
    # Remove any spaces/foreign characters present in the points column
    for char in repchar:
        pts.Point = pts.Point.apply(lambda x: x.replace(char,''))
    
    # Get inland points (if they exist)
    loop = ['seaward']
    inlandpts = pts[pts['Point'].str.contains('-')]  # USING '-'
    if len(inlandpts) < 1:
        pts = pts.astype({'Point':'float64'}, errors='ignore')
        inlandpts = pts.query('Point >999')              # For Leica labeling
        seaward = pts.query('Point >0')
    else:
        seaward = pts[~pts['Point'].str.contains('-')]
    
    if len(inlandpts) > 0:
        loop.append('inland')
    else:
        print 'No inland extent detected!'
    # print loop
    
    # Loop for seaward and inland extents
    for section in loop:
        print '\tWorking on ' + section 
        # Detect which section
        if section == 'seaward':
            p = seaward
            sheetname = 'Seaward_Extent'
            ascending = 1 # Sort by ascending
        elif section == 'inland':
            p = inlandpts
            sheetname = 'Inland_Extent'
            ascending = 0 # Sort by descending
            
        # Convert point column to float and sort
        p = p.astype({'Point':'float64'}, errors='ignore')
        ps = p.sort_values('Point', ascending=ascending)
    
        # Read out xy coordinates as lists
        northing = [starty] + ps['Northing_meters'].tolist()
        easting = [startx] + ps['Easting_meters'].tolist()
        elevations = [startz] + ps['Elevation_meters'].tolist()
    
        # Loop for each point
        dists = []
        for i in xrange(0,len(northing)):
            
            # Set first point as 0
            if i == 0:
                d1 = 0
                i +=1
                dists.append(d1)
                continue
            
            # Set point values
            x0 = easting[i-1]
            x1 = easting[i]
            y0 = northing[i-1]
            y1 = northing[i]
            d0 = dists[i-1]
            
            # Compute distance and append to list
            dx = x1-x0
            dy = y1-y0
            d1 = d0 + math.sqrt(math.pow(dx,2) + math.pow(dy,2))
            dists.append(d1)
            
#            # For last point in seaward, make sure elevation goes below 0
#            
#            if i == len(northing)-1:
#                if section == 'seaward':
#                    dz = elevations[i-1] - elevations[i-2] # height diff b/w last two points (ltp)
#                    dd = dists[i-1] - dists[i-2] # change in distance between ltp
#                    slope = dz/dd # slope (d_height, d_distance)
#                    rz = abs(thres) + elevations[i-2] - dz # remaining height till threshold
#                    interdist = rz * abs(1/slope) # (abs value needed?)
#                    dn = northing[i-1] - northing[i-2] # change in northing 
#                    azi = 90 - math.asin(dn/dd) # get azimuth
#                    internor = northing[i-1] + interdist*sin(azi) # new northing
#                    intereas = easting[i-1] + interdist*cos(azi) # new easting
#                    
#                    # Make separate df
#                    infillpt = pandas.DataFrame(columns = cols)
#                    infillpt['PointLabel'] = 'Interp_NAVD88'
#                    infillpt['Northing_meters'] = internor
#                    infillpt['Easting_meters'] = intereas
#                    infillpt['Elevation_meters'] = thres
#                    Infillpt['Distance']
                    
                    
            # add to counter
            i += 1
            
            
        # Readd benchmark row so distances will be properly placed.
        expps = bm.append(ps)
        
        # Reindex so that the distances will join correctly
        expps = expps.reset_index(drop=True)
        
        # Merge distance column into original dataframe
        expps['Distance'] = pandas.Series(dists)
        
        # Export frame into an excel spreadsheet
        outcsv = loccout + t + '_' + survey + '.xlsx'
        if os.path.exists(outcsv): 
            print '\tExporting to existing sheet!'
            writer = pandas.ExcelWriter(outcsv, engine='openpyxl')
            writer.book = openpyxl.load_workbook(outcsv)
        else:
            print '\tCreating new sheet!'
            writer = pandas.ExcelWriter(outcsv)
        expps.to_excel(writer, sheet_name=sheetname, columns=['PointLabel','Northing_meters','Easting_meters','Elevation_meters','Distance','Source'], encoding='utf-8')
        writer.save()
        print '\tExported csv for ' + section + ' extent of transect ' + t

        