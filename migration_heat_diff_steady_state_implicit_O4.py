import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import sparse
from scipy.sparse.linalg import spsolve

import geopandas as gpd
import rasterio
from rasterio.mask import raster_geometry_mask

from geospatial import raster_2array, get_raster_proj, array_2raster


def partial_k(k):

    k_x = (-k[4:, 2:-2] + 8 * k[3:-1, 2:-2] - 8 * k[1:-3, 2:-2] + k[:-4, 2:-2]) / 12
    k_y = (-k[2:-2, 4:] + 8 * k[2:-2, 3:-1] - 8 * k[2:-2, 1:-3] + k[2:-2, :-4]) / 12

    # pad derivatives
    k_x = np.pad(k_x, (2, 2), 'constant')
    k_y = np.pad(k_y, (2, 2), 'constant')

    return k_x, k_y


def get_raster_inds(shp_file, fname):

    # Grab vector file of region of interest
    df_shp = gpd.read_file(shp_file)
    extract_shp = df_shp[df_shp[shp_dic['colname']] == shp_dic['geoname']]

    # Convert extracted feature geometry to geojson for mask
    extract_shp = extract_shp['geometry'].values[0]

    # Use Rasterio to open raster file again
    with rasterio.open(os.path.join(base_dir, fname)) as src:
        crop_mask, crop_transform, window = raster_geometry_mask(
            src, [extract_shp], invert=True, all_touched=False
        )

    inds = np.where(crop_mask)
    inds_tpl = [(i, j) for i, j in zip(inds[0], inds[1])]

    return inds, inds_tpl


if __name__ == "__main__":

    # Files and directories
    ver = 0
    base_dir = '/data/vp/pop_migration/comb_rasters'
    source_shp = '/data/shapefiles/yemen_source.shp'
    # sink_shp1 = '/data/shapefiles/yemen_sink1.shp'
    # sink_shp2 = '/data/shapefiles/yemen_sink2.shp'
    shp_dic = {'colname': 'name',
               'geoname': 'initial_condition'}
    fname = f'yemen_resistance_smoothed_{ver}.tif'
    # fname = f'yemen_resistance_{ver}.tif'
    save_dir = f'/data/vp/pop_migration/figs_steady_state_{ver}'

    scl_fct = .001
    # new_scale = [0.01, 0.015]

    # Boundary temperatures
    t_cool = 0
    t_hot = 200

    # --------------------
    # Set source condition
    # --------------------
    source_inds, source_inds_tpl = get_raster_inds(source_shp, fname)

    # ------------------
    # Set sink condition
    # ------------------
    sink_inds1, sink_inds_tpl1 = get_raster_inds(sink_shp1, fname)
    sink_inds2, sink_inds_tpl2 = get_raster_inds(sink_shp2, fname)

    all_inds = source_inds_tpl  # + sink_inds_tpl1 + sink_inds_tpl2

    # ------------------------
    # Read in friction surface
    # ------------------------
    # Get Friction surface from raster
    k = raster_2array(os.path.join(base_dir, fname), band=1, replace_nodata_val=None)
    k = k.astype(float)

    # Replace missing values
    missing_inds = np.where(k < 0)

    # Some small fraction of values are 0
    # k[k == 0] = 0.01
    k[k < 1] = 1

    # Inverse friction to represent thermal conductivity
    k = 1 / k * scl_fct
    # k = ((new_scale[1] - new_scale[0]) * (k - k.min())) / (k.max() - k.min()) + new_scale[0]

    # Find partials of thermal conductivity
    k_x, k_y = partial_k(k)
    k_x[np.isnan(k_x)] = 0
    k_y[np.isnan(k_y)] = 0

    # Number of rows and columns
    n_rows = np.shape(k)[0]
    n_cols = np.shape(k)[1]

    # Construct the coefficient matrix
    A = sparse.eye(n_rows * n_cols)
    A = A.tolil()

    # find range ignoring top and bottom boundaries
    l_bound = 2 * n_cols
    u_bound = n_rows * n_cols - 2 * n_cols

    # -------------------------
    # Fill coefficient matrix
    # -------------------------
    for l in range(l_bound, u_bound):

        ind_mat = np.unravel_index(l, (n_rows, n_cols))
        i = ind_mat[0]
        j = ind_mat[1]

        # Ignore if part of source from shapefile
        if not ((i, j) in all_inds):

            # Ignore left right boundaries
            if (j > 2) and (j < n_cols - 2):
                # U_i+2,j term
                A[l, l+2*n_cols] = -1 / 12 * (k[i, j] + k_x[i, j])

                # U_i+1,j term
                A[l, l+n_cols] = 2 / 3 * (2 * k[i, j] + k_x[i, j])

                # U_i,j term
                A[l, l] = -5 * k[i, j]

                # U_i-1,j term
                A[l, l-n_cols] = 2 / 3 * (2 * k[i, j] - k_x[i, j])

                # U_i-2, j term
                A[l, l-2*n_cols] = - 1 / 12 * (k[i, j] - k_x[i, j])

                # U_i,j+2
                A[l, l+2] = -1 / 12 * (k[i, j] + k_y[i, j])

                # U_i,j+1
                A[l, l+1] = 2 / 3 * (2 * k[i, j] + k_y[i, j])

                # U_i,j-1
                A[l, l-1] = 2 / 3 * (2 * k[i, j] - k_y[i, j])

                # U_i,j-2
                A[l, l-2] = -1 / 12 * (k[i, j] - k_y[i, j])

    # --------------------------------------
    # Specify temperatures at boundaries (b)
    # --------------------------------------
    b = np.ones(n_rows * n_cols) * t_cool

    # -----------------------------
    # Specify temperature at source
    # -----------------------------
    source_inds_flat = np.ravel_multi_index(source_inds, (n_rows, n_cols))
    b[source_inds_flat] = t_hot

    # -------------------------
    # Solve system of equations
    # -------------------------
    A = A.tocsr()
    x = spsolve(A, b)

    x = np.reshape(x, (n_rows, n_cols))

    # --------------------
    # Save array as raster
    # --------------------
    ref_proj = get_raster_proj(os.path.join(base_dir, fname))
    array_2raster(ref_proj, x, fname_out=os.path.join(save_dir,
                  f'yemen_heat_steadystate_ver_{ver}.tif'))

    # -------------------------
    # Plot and save numpy array
    # -------------------------
    fig, ax = plt.subplots()
    im = ax.imshow(x, cmap=plt.get_cmap('viridis'), vmin=t_cool, vmax=t_hot)
    ax.set_axis_off()
    ax.set_title(f'Steady State Solution')
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    cbar_ax.set_xlabel('$T$ / K', labelpad=20)
    fig.colorbar(im, cax=cbar_ax)
    fig.savefig(os.path.join(save_dir, f"heat_steady_state_solution_ver_{ver}.svg"))
    plt.close(fig)
    np.save(os.path.join(save_dir, f"heat_array_steady_state_ver{ver}.npy"), x)
