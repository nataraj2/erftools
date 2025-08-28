import sys
import os
import argparse
from mpi4py import MPI
import glob

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import Download_ERA5_Data
from erftools.preprocessing import Download_ERA5_ForecastData
from erftools.preprocessing import Download_ERA5_SurfaceData
from erftools.preprocessing import Download_ERA5_ForecastSurfaceData
from erftools.preprocessing import ReadERA5_3DData
from erftools.preprocessing import ReadERA5_SurfaceData

from pyproj import CRS, Transformer
from numpy import *

def CreateLCCMapping(area):

    lat1 = area[2]
    lat2 = area[0]
    lon1 = area[1]
    lon2 = area[3]

    # Build CRS
    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2

    lat_1 = lat1 + delta/6
    lat_2 = lat2 - delta/6

    lambert_conformal = (
        f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
        f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"
    )

    return lambert_conformal
    
if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if len(sys.argv) == 1:
        print("Usage: python3 WriteICFromERA5Data.py <input_filename> [--do_forecast=true]")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Download and process ERA5 data.")
    parser.add_argument("input_filename", help="Input filename, e.g. inputs_for_Laura")
    parser.add_argument("--do_forecast", type=bool, default=False, help="Set to true to download forecast data")

    args = parser.parse_args()

    input_filename = args.input_filename
    do_forecast = args.do_forecast


    os.makedirs("ERA5Data/3D", exist_ok=True)
    os.makedirs("ERA5Data/Surface", exist_ok=True)

    os.makedirs("Output/3D", exist_ok=True)
    os.makedirs("Output/Surface", exist_ok=True)
    os.makedirs("Output/VTK/3D/ERA5Domain", exist_ok=True)
    os.makedirs("Output/VTK/3D/ERFDomain", exist_ok=True)
    os.makedirs("Output/VTK/Surface/ERA5Domain", exist_ok=True)
    os.makedirs("Output/VTK/Surface/ERFDomain", exist_ok=True)

   # Download surface data
    if do_forecast:
        print("Running forecast (pressure-level) download + processing...")
        filenames, area = Download_ERA5_ForecastSurfaceData(input_filename)
    else:
        print("Running surface download + processing...")
        filenames, area = Download_ERA5_SurfaceData(input_filename)
    comm.Barrier();

    lambert_conformal = CreateLCCMapping(area)

    # Each rank scans the local directory for matching files
    filenames = sorted(glob.glob("era5_surf_*.grib"))

    my_files = filenames[rank::size]

    for filename in my_files:
        print(f"[Rank {rank}] Processing file: {filename}")
        ReadERA5_SurfaceData(filename, lambert_conformal)

    comm.Barrier();
    if rank == 0:
        print("All ranks finished successfully. Exiting.", flush=True)
 
    # Download 3d data over pressure levels
    if do_forecast:
        filenames, area = Download_ERA5_ForecastData(input_filename)
        comm.Barrier();

        filenames = sorted(glob.glob("era5_3d_*.grib"))

        my_files = filenames[rank::size]
        # Create the directory if it doesn't exist
        for filename in my_files:
            print(f"Processing file: {filename}")
            ReadERA5_3DData(filename, lambert_conformal)

    else:
        filename, area = Download_ERA5_Data(input_filename)
        print("Filename is ", filename)
        print(f"Processing file: {filename}")
        ReadERA5_3DData(filename, lambert_conformal) 
