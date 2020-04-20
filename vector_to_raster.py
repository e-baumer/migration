from geospatial import get_raster_proj, rasterize_vectorfile


# shp_fname = "/data/shapefiles/yemen/yemen_admin1_1_buffer.shp"
# shp_fname = "/data/vp/pop_migration/yemen_admin2_severity.shp"
# ref_rast = '/data/vp/pop_migration/yemen_slope.tif'
# out_fname = '/data/vp/pop_migration/yemen_admin2_severity.tif'
shp_fname = "/data/shapefiles/yemen_source.shp"
ref_rast = '/data/vp/pop_migration/yemen_admin2_clipped.tif'
out_fname = '/data/vp/pop_migration/yemen_source.tif'


ref_proj = get_raster_proj(ref_rast)

# rasterize_vectorfile(ref_proj, shp_fname, attr_field='inter_clus', fname_out=out_fname)
rasterize_vectorfile(ref_proj, shp_fname, fname_out=out_fname)
