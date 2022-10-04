# Python code used in:
# Towards reproducible analysis of benthos structural complexity: 
# A case study on Antarctic polychaete reefs using action cameras and remotely operated vehicles

# Authors: J.C. Montes-Herrera, G. Johnstone, J. Stark, N. Hill, V. Cummings, V. Lucieer

# Contact: juancarlos.montesherrera@utas.edu.au

# Data: https://doi.org/10.5281/zenodo.7115132

##### Load raster files ######

#Spatial resolution export:
## Orthomosaic: 0.5 mm/pix
## DEM: 1 cm/pix
import os
from osgeo import gdal

os.chdir(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Photogrammetry\AAD2019_EF_PolychaeteReef\processed_data\EF-G-AggressiveFilter\EF-G-PolygonExport\EF-G-04")
out_fp = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Photogrammetry\AAD2019_EF_PolychaeteReef\processed_data\EF-G-AggressiveFilter\EF-G-AgressiveFilter-Tiles\EF-G-04/"

dem = gdal.Open(r"C:\Users\jcmontes\Desktop\test_tiles\DEM_reef_slope.tif")

gt = dem.GetGeoTransform()

xmin = gt[0]
ymax = gt[3]
res = gt[1]

# raster lengths
xlen = res * dem.RasterXSize
ylen = res * dem.RasterYSize

# number of tiles in x and y direction
xdiv = int(xlen)
ydiv = int(ylen)

# size of single tile
xsize = xlen/xdiv
ysize = ylen/ydiv

# create lists of x and y coordinates
xsteps = [xmin + xsize * i for i in range(xdiv+1)]
ysteps = [ymax - ysize * i for i in range(ydiv+1)]

for i in range(xdiv):
    for j in range(ydiv):
        xmin = xsteps[i]
        xmax = xsteps[i+1]
        ymax = ysteps[j]
        ymin = ysteps[j+1]
        
        print("xmin: "+str(xmin))
        print("xmax: "+str(xmax))
        print("ymin: "+str(ymin))
        print("ymax: "+str(ymax))
        print("\n")
        
        # use gdal warp
        #gdal.Warp("dem"+str(i)+str(j)+".tif", dem, outputBounds = (xmin, ymin, xmax, ymax), dstNodata = -9999)
        
        # or gdal translate to subset the input raster
        gdal.Translate(out_fp + "dem_tile"+str(i)+str(j)+".tif", dem, projWin = (xmin, ymax, xmax, ymin), xRes = res, yRes = -res)
        
#close the open dataset to flush memory
dem = None

# Now the orthomosaic

ortho = gdal.Open(r"C:\Users\jcmontes\Desktop\EF-07_DEM_Ortho\ortho_1.tif")

gt = ortho.GetGeoTransform()

xmin = gt[0]
ymax = gt[3]
res = gt[1]

# raster lengths
xlen = res * ortho.RasterXSize
ylen = res * ortho.RasterYSize

# number of tiles in x and y direction
xdiv = int(xlen)
ydiv = int(ylen)

# size of single tile
xsize = xlen/xdiv
ysize = ylen/ydiv

# create lists of x and y coordinates
xsteps = [xmin + xsize * i for i in range(xdiv+1)]
ysteps = [ymax - ysize * i for i in range(ydiv+1)]

for i in range(xdiv):
    for j in range(ydiv):
        xmin = xsteps[i]
        xmax = xsteps[i+1]
        ymax = ysteps[j]
        ymin = ysteps[j+1]
        
        print("xmin: "+str(xmin))
        print("xmax: "+str(xmax))
        print("ymin: "+str(ymin))
        print("ymax: "+str(ymax))
        print("\n")
        
        # use gdal warp
        #gdal.Warp("dem"+str(i)+str(j)+".tif", dem, outputBounds = (xmin, ymin, xmax, ymax), dstNodata = -9999)
        
        # or gdal translate to subset the input raster
        gdal.Translate(out_fp + "ortho_tile"+str(i)+str(j)+".tif", ortho, projWin = (xmin, ymax, xmax, ymin), xRes = res, yRes = -res)
       
#close the open dataset to flush memory
ortho = None