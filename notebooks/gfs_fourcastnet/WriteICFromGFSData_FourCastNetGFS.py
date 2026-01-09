import sys
import os
from mpi4py import MPI
from numpy import array_split
import argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import ReadGFS_3DData_FourCastNetGFS


def CreateLCCMapping(area):
    lat1 = area[2]
    lat2 = area[0]
    lon1 = area[1]
    lon2 = area[3]

    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2
    lat_1 = lat1 + delta / 6
    lat_2 = lat2 - delta / 6

    lambert_conformal = (
        f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
        f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"
    )

    return lambert_conformal


def get_area(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("area:"):
                area_str = line.split("area:")[1].strip()
                area = [float(x) for x in area_str.split(",")]
                return area
    raise ValueError("Area not found in the file.")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate initial conditions from GFS data (MPI-parallel)"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to input GFS file"
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="Input directory containing GFS data"
    )
    parser.add_argument(
        "--init_hurricane_latlon",
        type=str,
        default=None,
        help="Optional initial hurricane lat,lon (e.g. --init_hurricane_latlon=25,-80)"
    )
    return parser.parse_args()



if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = parse_args()

    input_file = args.input_file
    input_dir  = args.input_dir

    # Parse lat/lon if provided
    if args.init_hurricane_latlon is not None:
        if size > 1:
            print("Error: When hurricane eye tracking is turned only serial run is allowed. Use only 1 MPI rank.")
            sys.exit(1)
        try:
            init_lat, init_lon = map(float, args.init_hurricane_latlon.split(","))
            init_lon = 360.0+init_lon
        except Exception:
            if rank == 0:
                print("Error: --init_hurricane_latlon must be of the form 'lat,lon', e.g. --init_hurricane_latlon=25,-80")
            sys.exit(1)
        
    else:
        init_lat = init_lon = -1e3  # Optional case

    

    if rank == 0:
        print(f"Input file: {input_file}")
        print(f"Input dir:  {input_dir}")
        if init_lat is not None:
            print(f"Initial hurricane location: lat={init_lat}, lon={init_lon}")
        else:
            print("Initial hurricane location not provided.")

    # Rank 0 reads input and prepares file list
    if rank == 0:
        if not os.path.isdir(input_dir):
            print(f"Error: {input_dir} is not a valid directory")
            sys.exit(1)

        area = get_area(input_file)
        lambert_conformal = CreateLCCMapping(area)

        file_list = sorted(
            f for f in os.listdir(input_dir)
            if os.path.isfile(os.path.join(input_dir, f))
        )

        os.makedirs("Output/GFSData_3D", exist_ok=True)
        os.makedirs("Output/VTK/3D/GFSDomain", exist_ok=True)
        os.makedirs("Output/VTK/3D/ERFDomain", exist_ok=True)
        os.makedirs("Output/track_latlon", exist_ok=True)

    else:
        file_list = None
        area = None
        lambert_conformal = None

    # Broadcast info to all ranks
    file_list = comm.bcast(file_list, root=0)
    area = comm.bcast(area, root=0)
    lambert_conformal = comm.bcast(lambert_conformal, root=0)

    # Split files evenly among ranks
    file_chunks = array_split(file_list, size)
    my_files = file_chunks[rank]

    print(f"[Rank {rank}] Processing {len(my_files)} files", flush=True)

    filename = "Output/hurricane_track_latlon.txt"

    if rank == 0:
        # Delete if it exists
        if os.path.exists(filename):
            os.remove(filename)
        # Create an empty file
        open(filename, "w").close()

    # Ensure all ranks wait until root is done
    comm.Barrier()

    track_points = []

    for filename in my_files:
        full_path = os.path.join(input_dir, filename)
        print(f"[Rank {rank}] Processing {filename}", flush=True)
        try:
            ReadGFS_3DData_FourCastNetGFS(full_path, area, lambert_conformal, init_lon, init_lat, track_points)
        except Exception as e:
            print(f"[Rank {rank}] Error processing {filename}: {e}", flush=True)

    comm.Barrier()
    if rank == 0:
        print("All ranks completed successfully.", flush=True)

