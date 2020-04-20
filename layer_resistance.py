import json
import numpy as np
import os
from scipy.ndimage.filters import gaussian_filter

from geospatial import get_raster_proj, raster_2array, array_2raster


# ----------------
# Layer file names
# ----------------
slope_rast = '/data/vp/pop_migration/yemen_slope_clipped_resized.tif'
admin2_rast = '/data/vp/pop_migration/yemen_admin2_clipped_resized.tif'
severity_rast = '/data/vp/pop_migration/yemen_severity_clipped_resized.tif'
roads_rast = '/data/vp/pop_migration/yemen_roads_clipped_resized.tif'
water_rast = '/data/vp/pop_migration/yemen_waterways_clipped_resized.tif'

# slope_rast = '/data/vp/pop_migration/yemen_slope_clipped.tif'
# admin2_rast = '/data/vp/pop_migration/yemen_admin2_clipped.tif'
# severity_rast = '/data/vp/pop_migration/yemen_severity_clipped.tif'
# roads_rast = '/data/vp/pop_migration/yemen_roads_clipped.tif'
# water_rast = '/data/vp/pop_migration/yemen_waterways_clipped.tif'

out_dir = '/data/vp/pop_migration/comb_rasters_red'
ver_fname = 'yemen_raster_version_history.txt'

new_scale = [1, 1000]
sigma = [1.0, 1.0]

# -----------------------------------------
# Get reference project, use slope for this
# -----------------------------------------
ref_proj = get_raster_proj(slope_rast)

# -------------------------
# Read in rasters as arrays
# -------------------------
slope_array = raster_2array(slope_rast, replace_nodata_val=np.nan)
admin2_array = raster_2array(admin2_rast, replace_nodata_val=np.nan)
severity_array = raster_2array(severity_rast, replace_nodata_val=np.nan)
roads_array = raster_2array(roads_rast, replace_nodata_val=np.nan)
water_array = raster_2array(water_rast, replace_nodata_val=np.nan)

slope_array = np.squeeze(slope_array)
admin2_array = np.squeeze(admin2_array)
severity_array = np.squeeze(severity_array)
roads_array = np.squeeze(roads_array)
water_array = np.squeeze(water_array)

# -------------------------------------------------
# Merge raster files. Slope is the base layer.
# Roads, waterways, getsliced in. Admin2 boundaries
# and severity is added to the base slope layer.
# -------------------------------------------------
admin2_inds = np.where(admin2_array == 255)
roads_inds = np.where(roads_array == 255)
water_inds = np.where(water_array == 255)

# --------------------------------------------
# Loop through various scenerios of resistance
# --------------------------------------------
version = 0
ver_dic = {}
iend = 100
for i in [1, 5, iend]:
    for j in [0, 70, 200]:
        for scl_fct in [0, 3, 6]:

            # Load in slope data
            comb_array = slope_array.copy()

            # Convert slope from degrees to ratio
            comb_array = np.tan(comb_array * np.pi / 180)

            # Calculate resistance from slope
            comb_array = np.exp(3 * comb_array) + 5

            # Rescale slope data between 1 and 100
            # comb_array = (
            #     (new_scale[1] - new_scale[0]) * (comb_array - comb_array.min())
            # ) / (comb_array.max() - comb_array.min()) + new_scale[0]

            # Water ways are always not desirable to pass through
            comb_array[water_inds] = 400

            # Set resistance values for roads
            if i == 1:
                comb_array[roads_inds] = i
            elif i == iend:
                comb_array[roads_inds] = comb_array[roads_inds] + i

            # Add resistance from admin2 boundaries
            comb_array[admin2_inds] = comb_array[admin2_inds] + j

            # Add resistance from needs severity
            if scl_fct == 0:
                comb_array = comb_array + severity_array * scl_fct
            else:
                temp_sev = severity_array - 3
                comb_array = comb_array + temp_sev * scl_fct

            comb_array[comb_array > 1000] = 1000
            comb_array[comb_array < 1] = 1

            # Smooth array
            smoothed_array = gaussian_filter(comb_array, sigma, mode='reflect')

            # Write combinded rasters to file
            array_2raster(ref_proj, comb_array,
                          fname_out=os.path.join(out_dir, f'yemen_resistance_{version}.tif'))

            # Write combinded rasters to file
            array_2raster(
                ref_proj, smoothed_array,
                fname_out=os.path.join(out_dir, f'yemen_resistance_smoothed_{version}.tif')
            )

            # Print combo
            ver_dic[version] = {'roads_factor': i, 'admin2_factor': j, 'severity_scl': scl_fct}

            version += 1

# Write version history to text file
with open(os.path.join(out_dir, ver_fname), 'w') as fid:
    fid.write(json.dumps(ver_dic))
