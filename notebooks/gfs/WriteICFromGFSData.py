import sys
import os
import argparse
from mpi4py import MPI
import glob

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import Download_GFS_Data
from erftools.preprocessing import Download_GFS_ForecastData
from erftools.preprocessing import ReadGFS_3DData
from erftools.preprocessing import ReadGFS_3DData_UVW

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
        print("Usage: python3 WriteICFromGFSData.py <input_filename> [--do_forecast=true]")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Download and process GFS data.")
    parser.add_argument("input_filename", help="Input filename, e.g. inputs_for_Laura")
    parser.add_argument("--do_forecast", type=bool, default=False, help="Set to true to download forecast data")

    args = parser.parse_args()

    input_filename = args.input_filename
    do_forecast = args.do_forecast

    os.makedirs("Output", exist_ok=True)
    os.makedirs("Output/GFSData_3D", exist_ok=True)
    os.makedirs("Output/VTK/3D/GFSDomain", exist_ok=True)
    os.makedirs("Output/VTK/3D/ERFDomain", exist_ok=True)

    if do_forecast:
        filenames, area = Download_GFS_ForecastData(input_filename)
        lambert_conformal = CreateLCCMapping(area)
        comm.Barrier();

        filenames = sorted(glob.glob("gfs.0p25.*.grib2"))
        print(f"[Rank {rank}] Found {len(filenames)} grib2 files", flush=True)

        my_files = filenames[rank::size]
        # Create the directory if it doesn't exist
        for filename in my_files:
            print(f"[Rank {rank}] Entering ReadGFS_3DData for {filename}", flush=True)

            ReadGFS_3DData(filename, area, lambert_conformal)
        
    else:
        filename, area = Download_GFS_Data(input_filename)
        lambert_conformal = CreateLCCMapping(area)
        print("Filename is ", filename)
        print(f"Processing file: {filename}")
        ReadGFS_3DData(filename, area, lambert_conformal)
        #ReadGFS_3DData_UVW(filename, area, lambert_conformal)
