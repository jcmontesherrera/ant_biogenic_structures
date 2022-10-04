# Metashape Python API code used in:
# Towards reproducible analysis of benthos structural complexity: 
# A case study on Antarctic polychaete reefs using action cameras and remotely operated vehicles

# Authors: J.C. Montes-Herrera, G. Johnstone, J. Stark, N. Hill, V. Cummings, V. Lucieer

# Contact: juancarlos.montesherrera@utas.edu.au

import Metashape

chunk = Metashape.app.document.chunk

#specify directory to save polygons
directory = r"C:\Users\jcmontes\Desktop/"

##Spatial resolution export:
### Orthomosaic: 0.5 mm/pix
### DEM: 1 cm/pix
### DEM_lowres: To match dimensions of ortho raster for automated classification tests

## Export jpeg from metashape

compression = Metashape.ImageCompression()
compression.jpeg_quality = 90
compression.tiff_big = False
jpeg_format = Metashape.ImageFormat.ImageFormatJPEG

for shape in chunk.shapes:
    shape.boundary_type = Metashape.Shape.BoundaryType.OuterBoundary
    path_ortho = directory + "ortho_" + shape.label + ".tif" #if there are unique shape.labels, otherwise you can use str(shape.key)
    path_ortho_lowres = directory + "ortho_" + shape.label + "_lowres.tif"
    path_dem = directory + "DEM_" + shape.label + ".tif"
    path_dem_hires = directory + "DEM_" + shape.label + "_hires.tif"
    path_jpg = directory + "ortho_" + shape.label + "_jpg.jpeg"
    ## Changing the orthomosaic resolution - Testing classification algorithms require raster to be same dimensions
    ## Full resolution export - required for human annotation  
    chunk.exportRaster(path_ortho, resolution_x=0.001, resolution_y=0.001, clip_to_boundary = True, source_data = Metashape.OrthomosaicData)
    chunk.exportRaster(path_dem, resolution_x=0.001, resolution_y=0.001, clip_to_boundary = True, source_data = Metashape.ElevationData)
    ## Export for different resolution between orthomosaic and DEMs/Structural metrics
    #chunk.exportRaster(path_ortho_lowres, resolution_x=0.001, resolution_y=0.001, clip_to_boundary = True, source_data = Metashape.OrthomosaicData)
    #chunk.exportRaster(path_dem_hires, resolution_x=0.0005, resolution_y=0.0005, clip_to_boundary = True, source_data = Metashape.ElevationData)
    ## JPEG
    #chunk.exportRaster(path_jpg, resolution_x=0.0025, resolution_y=0.0025, image_format=jpeg_format, image_compression=compression, source_data = Metashape.OrthomosaicData)
    shape.boundary_type = Metashape.Shape.BoundaryType.NoBoundary

print('Done')

