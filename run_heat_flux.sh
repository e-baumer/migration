#! /bin/bash

# python3
p3=/usr/bin/python3

# heat flux
hf=/home/ebaumer/Code/vp/pop_migration/calculate_heat_flux.py

# Local Dir
DATA_DIR=$1

declare -a ver=("5" "6" "7" "8" "9" "10"
                "11" "12" "13" "14" "15" "16" "17" "18"
                "19" "20" "21" "22" "23" "24" "25" "26")
# declare -a ver=("4")

# Options for gdal command (specify no data values)
# GDALT_OPTS="-of PNG -scale"
GDALT_OPTS="-of JPEG -scale"

for v in "${ver[@]}"
do
    rfile="/data/vp/pop_migration/comb_rasters_red/yemen_resistance_smoothed_${v}.tif"
    odir="/data/vp/pop_migration/heat_flux_${v}"
    tdir="/data/vp/pop_migration/figs_rk45_${v}_red"
    $p3 $hf $rfile $odir 0.01 "--tempdir=${tdir}"
done

# for f in `ls -v *.jpg`; do echo "file '$f'" >> list.txt; done
# ffmpeg -f concat -i list.txt heatflux_ver0.avi
# GDALM_OPTS="-of GTiff -ot Byte -a_nodata 0 -scale -.2 1 1 255"

# ffmpeg -f concat -i list.txt -vcodec libx264 -acodec aac output.mp4
# ffmpeg -f concat -i list.txt -c:v libx264 -pix_fmt yuv420p output3.mp4

# for y in "${years[@]}"
# do
#     INFILE=$DATA_BASE_DIR/"modis_${y}_NDVI.tif"
#     OUTFILE=$DATA_BASE_DIR/"modis_${y}_NDVI_byte.tif"

#     $gdalm $GDALM_OPTS $INFILE $OUTFILE
