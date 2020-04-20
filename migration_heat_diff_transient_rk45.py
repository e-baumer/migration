import click
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import re

import geopandas as gpd
from osgeo import gdal
import rasterio
from rasterio.mask import raster_geometry_mask
from osgeo import gdal

from geospatial import raster_2array, get_raster_proj, array_2raster


def partial_k(k):
    k_x = k.copy()
    k_y = k.copy()

    k_x = 0.5 * (k[2:, 1:-1] - k[:-2, 1:-1])
    k_y = 0.5 * (k[1:-1, 2:] - k[1:-1, :-2])

    return k_x, k_y


def f_star(u, k, k_x, k_y):
    f = u.copy()
    f[1:-1, 1:-1] = (
        (k[1:-1, 1:-1] + 0.5 * k_x) * u[2:, 1:-1] -
        4 * k[1:-1, 1:-1] * u[1:-1, 1:-1] +
        (k[1:-1, 1:-1] - 0.5 * k_x) * u[:-2, 1:-1] +
        (k[1:-1, 1:-1] + 0.5 * k_y) * u[1:-1, 2:] +
        (k[1:-1, 1:-1] - 0.5 * k_y) * u[1:-1, :-2]
    )

    return f


def create_output(t_cool, t_hot, save_dir, ver, i, u_plot):

    fig, ax = plt.subplots()
    im = ax.imshow(u_plot, cmap=plt.get_cmap('viridis'), vmin=t_cool, vmax=t_hot)
    ax.set_axis_off()
    ax.set_title(f'Iteration {i}')
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    fig.savefig(os.path.join(save_dir, f"heat_fig_ver_{ver}_iter_{i:.0f}.svg"))
    plt.close(fig)
    np.save(os.path.join(save_dir, f"heat_array_ver_{ver}_iter_{i:.0f}.npy"), u_plot)


@click.group()
def cli():
    pass


@click.command()
@click.argument('i_ver', type=int)
@click.argument('f_ver', type=int)
@click.argument('init_flg', default=False, type=bool)
def rk45(i_ver, f_ver, init_flg):
    # ---------------------
    # Files and directories
    # ---------------------
    for ver in range(i_ver, f_ver):
        base_dir = '/data/vp/pop_migration/comb_rasters_red'
        shp_file = '/data/shapefiles/yemen_ic.shp'
        shp_dic = {'colname': 'name',
                   'geoname': 'initial_condition'}
        fname = f'yemen_resistance_smoothed_{ver}.tif'
        # fname = f'yemen_resistance_{ver}.tif'
        save_dir = f'/data/vp/pop_migration/figs_rk45_{ver}_red'

        # ---------
        # Constants
        # ---------
        # Scale factor to apply to resistance surface
        scl_fct = 0.01

        # Acceptable tolerance for error
        e_tol = 4e-3

        # Max allowable time step
        # del_t_max = 1.0
        del_t_max = 60.0
        # Min allowable time step
        del_t_min = 0.0001
        # Initial time step to try
        del_t = 0.01

        # Time interval for output
        t_out_int = 25000
        # Iteration interval for output
        i_int = 10

        # Upper bound for del_t scaling factor
        s_upper = 4.0
        # Lower bound for del_t scaling factor
        s_lower = 0.125

        # Max number of iterations
        max_iters = 100000

        # Set the boundary conditions
        t_cool = 0
        t_hot = 200

        # Test if save directory exists, if not create it
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        # Get projection from reference raster
        ref_proj = get_raster_proj(os.path.join(base_dir, fname))

        # ----------------------------
        # Find region to set condition
        # ----------------------------
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

        source_inds = np.where(crop_mask)

        # ------------------------
        # Read in friction surface
        # ------------------------
        # Get Friction surface from raster
        k = raster_2array(os.path.join(base_dir, fname), band=1, replace_nodata_val=None)
        k = k.astype(float)

        # Inverse friction to represent thermal conductivity
        k = 1 / k * scl_fct

        # Find partials of thermal conductivity
        k_x, k_y = partial_k(k)
        k_x[np.isnan(k_x)] = 0
        k_y[np.isnan(k_y)] = 0

        # Number of rows and columns
        n_rows = np.shape(k)[0]
        n_cols = np.shape(k)[1]

        # Initialize temperature matrix
        if not init_flg:
            # Load latest temperature field and set as ic
            flist = [os.path.basename(f) for f in glob.glob(os.path.join(save_dir, "*.tif"))]
            regex = re.compile(r'\d+')
            flist.sort(key=lambda x: int(regex.findall(x)[1]))
            last_f = os.path.join(save_dir, flist[-1])

            t_total = int(regex.findall(flist[-1])[1])
            t_out_ref = t_total

            u_0 = raster_2array(last_f, band=1)

        else:
            u_0 = t_cool * np.ones((n_rows, n_cols))

            # Set initial conditions
            u_0[source_inds] = t_hot

            # -----------
            # Calculation
            # -----------
            t_total = 0
            t_out_ref = 0

        for i in range(max_iters):

            # Calculate first slope
            k_1 = del_t * f_star(u_0, k, k_x, k_y)

            # Calculate intermediate estimate of function (1/4)
            u_1 = u_0 + 0.25 * k_1

            # Calculate second slope
            k_2 = del_t * f_star(u_1, k, k_x, k_y)

            # Calculate intermediate estimate of function (3/8)
            u_2 = u_0 + (3 / 32) * k_1 + (9 / 32) * k_2

            # Calculate third slope
            k_3 = del_t * f_star(u_2, k, k_x, k_y)

            # Calculate intermediate estimate of function (12/13)
            u_3 = u_0 + (1932 / 2197) * k_1 - (7200 / 2197) * k_2 + (7296 / 2197) * k_3

            # Calculate fourth slope
            k_4 = del_t * f_star(u_3, k, k_x, k_y)

            # Calculate intermediate estimate of function (12/13)
            u_4 = u_0 + (439 / 216) * k_1 - 8 * k_2 + (3680 / 513) * k_3 - (845 / 4104) * k_4

            # Calculate fifth slope
            k_5 = del_t * f_star(u_4, k, k_x, k_y)

            # Calculate intermediate estimate of function (12/13)
            u_5 = (
                u_0
                - (8 / 27) * k_1
                + 2 * k_2
                - (3544 / 2565) * k_3
                + (1859 / 4104) * k_4
                - (11 / 40) * k_5
            )

            # Calculate sixth slope
            k_6 = del_t * f_star(u_5, k, k_x, k_y)

            # 4th order approximation
            u4_update = u_0 + (25 / 216) * k_1 + (1408 / 2565) * k_3 + (2197 / 4101) * k_4 - 0.2 * k_5

            # 5th order approximation
            u5_update = (
                u_0
                + (16 / 135) * k_1
                + (6656 / 12825) * k_3
                + (28561 / 56430) * k_4
                - (9 / 50) * k_5
                + (2 / 55) * k_6
            )

            # ----------------------
            # Find optimal step size
            # ----------------------
            # Find error between 4th and 5th order RK
            u_diff = np.abs(u4_update - u5_update)
            err_cal = u_diff.max()

            # Calculate scaling factor
            s = np.power(((e_tol * del_t) / (2 * err_cal)), 0.25)

            # Check value of s
            s = np.min((s, s_upper))
            s = np.max((s, s_lower))

            # Apply scaling factor to time-step
            del_t_new = s * del_t

            # Check if calculated time-step is in bounds
            del_t_new = np.min((del_t_new, del_t_max))
            del_t_new = np.max((del_t_new, del_t_min))

            # -------------------------------------------
            # Check calculated error compared to min
            # error tolerance to decide if to accept step
            # -------------------------------------------
            # import ipdb; ipdb.set_trace()
            if ((err_cal / del_t) < e_tol) or (del_t == del_t_min):
                t_total += del_t
                # u_0 = u4_update.copy()
                u_0 = u5_update.copy()

            del_t_old = del_t
            del_t = del_t_new.copy()

            if i % i_int == 0:
                print(f'Completed Iteration = {i}')
                print(f'Total time = {t_total}')
                print(f'Calculated error per second = {err_cal / del_t_old}')
                print(f'Current time step = {del_t}')
                print(f'Minimum value of array = {u_0.min()}')
                print(f'Maximum value of array = {u_0.max()}')
                print(f'---------------------------------------------------')

            # --------------------------
            # Output array, svg, and tif
            # --------------------------
            if (t_total >= (t_out_ref + t_out_int)) or (i == 0):
                t_out_ref = t_total
                t_out = np.round(t_total, decimals=0)

                print(f'Saving output for time {t_total}')
                print(f'---------------------------------------------------')
                # Save array as raster
                u_plot = np.reshape(u4_update.copy(), (n_rows, n_cols))
                array_2raster(
                    ref_proj, u_plot,
                    fname_out=os.path.join(save_dir,
                                           f'yemen_heat_transient_ver_{ver}_time_{t_out:.0f}.tif')
                )

                create_output(t_cool, t_hot, save_dir, ver, t_out, u_plot)


if __name__ == "__main__":
    rk45()
