import click
import glob
import numpy as np
import matplotlib.pyplot as plt
import os

from geospatial import raster_2array, get_raster_proj, array_2raster


def partials(u):

    u_x1 = 0.5 * (u[2:, 1:-1] - u[1:-1, 1:-1])
    u_x2 = 0.5 * (u[1:-1, 1:-1] - u[:-2, 1:-1])
    u_y1 = 0.5 * (u[1:-1, 2:] - u[1:-1, 1:-1])
    u_y2 = 0.5 * (u[1:-1, 1:-1] - u[1:-1, :-2])

    # pad derivatives
    u_x1 = np.pad(u_x1, (1, 1), 'constant')
    u_x2 = np.pad(u_x2, (1, 1), 'constant')
    u_y1 = np.pad(u_y1, (1, 1), 'constant')
    u_y2 = np.pad(u_y2, (1, 1), 'constant')

    return u_x1, u_x2, u_y1, u_y2


def create_output(save_dir, fname, u_plot):

    fig, ax = plt.subplots()
    im = ax.imshow(u_plot, cmap=plt.get_cmap('viridis'))
    ax.set_axis_off()
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    fig.savefig(os.path.join(save_dir, 'heat_flux_'+os.path.splitext(fname)[0]+".svg"))
    plt.close(fig)
    np.save(os.path.join(save_dir, 'heat_flux'+os.path.splitext(fname)[0]+".npy"), u_plot)


@click.group()
def cli():
    pass


@click.command()
@click.argument('friction', type=str)
@click.argument('save_dir', type=str)
@click.argument('scl_fct', default=0.01, type=float)
@click.option('--tempdir', default='', type=str)
@click.option('--fname', default='', type=str)
def calc_heat_flux(friction, save_dir, scl_fct, tempdir, fname):

    if tempdir:
        temp_files = glob.glob(os.path.join(tempdir, "*.tif"))

    if fname:
        temp_files = [fname]

    for f in temp_files:

        # ------------------------------------
        # Get projection from reference raster
        # ------------------------------------
        ref_proj = get_raster_proj(f)

        # ------------------------------------
        # Determine which friction file to use
        # ------------------------------------
        if os.path.isdir(friction):
            base_name = os.path.basename(f)
            split_name = os.path.splitext(base_name)[0]
            ver = split_name.split('_')[-1]

            friction_file = os.path.join(friction, 'yemen_resistance_smoothed_'+ver+'.tif')

        else:
            friction_file = friction

        # ---------------------------------
        # Read friction surface from raster
        # ---------------------------------
        k = raster_2array(friction_file, band=1, replace_nodata_val=None)
        k = k.astype(float)

        # Inverse friction to represent thermal conductivity
        k = 1 / k * scl_fct

        # ----------------------------
        # Read temperature from raster
        # ----------------------------
        u = raster_2array(f, band=1, replace_nodata_val=None)

        # Calculate derivatives of temperature
        u_x1, u_x2, u_y1, u_y2 = partials(u)

        # Calculate heat flux
        u_flux = k * (np.abs(u_x1 + u_x2) + np.abs(u_y1 + u_y2))

        # ------------------------------
        # Save heat flux array as raster
        # ------------------------------
        f_base = os.path.basename(f)
        array_2raster(
            ref_proj, u_flux,
            fname_out=os.path.join(save_dir, 'heat_flux_'+os.path.basename(f_base))
        )
        create_output(save_dir, f_base, u_flux)


if __name__ == "__main__":
    calc_heat_flux()
