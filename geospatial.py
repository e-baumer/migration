import numpy as np
from osgeo import gdal, osr, ogr


def raster_2array(raster_file, band=None, replace_nodata_val=None):
    '''
    Convert raster file to numpy array
    '''

    # Open raster file
    raster = gdal.Open(raster_file)
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    no_data = raster.GetRasterBand(1).GetNoDataValue()

    if band is None:
        rast_array = np.zeros((raster.RasterCount, rows, cols))
        for b in range(raster.RasterCount):
            b += 1
            rast_array[b-1, :, :] = raster.GetRasterBand(b).ReadAsArray()

    else:
        band = raster.GetRasterBand(band)
        rast_array = band.ReadAsArray()

    # Close raster file
    raster = None

    if not (replace_nodata_val is None):
        rast_array[rast_array == no_data] = replace_nodata_val

    return rast_array


def get_raster_proj(raster_file):
    '''
        Get projection and spatial system of reference raster file:
        -- Spatial reference system
        -- Origin X,Y
        -- Pixel width, height
        -- Number rows, cols

        This data is used in creating the new raster file with the same spatial characteristics
        as the reference raster file.
    '''
    projdata = {}

    # ----------------------------------------------------------
    # Opening the raster file can not be a separate Luigi task,
    # because the return object of gdal read is of Swig type
    # which cannot be pickled. Also, maybe it doesn't make sense
    # passing a large raster file as a parameter
    # ----------------------------------------------------------
    raster = gdal.Open(raster_file)

    # ----------------------------
    # Get spatial reference system
    # ----------------------------
    srs_wkt = raster.GetProjectionRef()

    # ----------------------------
    # Get grid (pixel) coordinates
    # ----------------------------
    geotransform = raster.GetGeoTransform()
    originx = geotransform[0]
    originy = geotransform[3]
    pixelwidth = geotransform[1]
    pixelheight = geotransform[5]

    # ------------------------------
    # Find number of rows and pixels
    # ------------------------------
    ncols = raster.RasterXSize
    nrows = raster.RasterYSize
    raster = None

    projdata['srs'] = srs_wkt
    projdata['pixel'] = (originx, pixelwidth, 0, originy, 0, pixelheight)
    projdata['ncolsrows'] = (ncols, nrows)

    return projdata


def array_2raster(ref_proj, ref_array, no_data_val=-9999, proj_type='srs', fname_out='out.tif'):
    '''Converts a numpy array to a raster file (geo-tiff)'''

    (ncols, nrows) = ref_proj['ncolsrows']

    driver = gdal.GetDriverByName('GTiff')

    out_raster = driver.Create(fname_out, ncols, nrows, 1, gdal.GDT_Float32)
    out_raster.SetGeoTransform(ref_proj['pixel'])
    outband = out_raster.GetRasterBand(1)
    outband.SetNoDataValue(no_data_val)
    outband.WriteArray(ref_array, 0, 0)

    proj = osr.SpatialReference()

    if proj_type.lower() == 'wkp':
        proj.SetWellKnownGeogCS(ref_proj['wkp'])

    elif proj_type.lower() == 'srs':
        proj.ImportFromWkt(ref_proj['srs'])

    else:
        raise Exception("Unrecongnized projection type for creating output raster file, "
                        "must be wkp or srs")

    out_raster.SetProjection(proj.ExportToWkt())
    outband.FlushCache()

    driver = None
    out_raster = None
    outband = None


def rasterize_vectorfile(ref_proj, vector_file, attr_field=None, no_data_val=-9999,
                         fname_out='out.tif', proj_type='srs'):
    """
    Base class for rasterizing a vector file.
    -- Vector file should be an ERSI Shapefile
    -- Raster file will be a geotiff

    TO DO: Add support for other raster and vector file types
    """
    # prj_wkt = '+proj=longlat +datum=WGS84 +no_defs'
    (ncols, nrows) = ref_proj['ncolsrows']

    driver = ogr.GetDriverByName("ESRI Shapefile")

    data_source = driver.Open(vector_file, 0)
    layer = data_source.GetLayer()

    target_ds = gdal.GetDriverByName('GTiff').Create(
        fname_out, ncols, nrows, 1, gdal.GDT_Float32
    )

    target_ds.SetGeoTransform(ref_proj['pixel'])

    proj = osr.SpatialReference()
    if proj_type.lower() == 'wkp':
        proj.SetWellKnownGeogCS(ref_proj['wkp'])

    elif proj_type.lower() == 'srs':
        proj.ImportFromWkt(ref_proj['srs'])

    else:
        raise Exception("Unrecongnized projection type for creating output raster file, "
                        "must be wkp or srs")

    target_ds.SetProjection(proj.ExportToWkt())

    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(no_data_val)
    band.Fill(no_data_val)
    band.FlushCache()

    if attr_field is None:
        gdal.RasterizeLayer(target_ds, [1], layer)
    else:
        gdal.RasterizeLayer(target_ds, [1], layer,
                            options=["ATTRIBUTE=%s" % attr_field])

    data_source = None
    layer = None
    target_ds = None
    band = None
