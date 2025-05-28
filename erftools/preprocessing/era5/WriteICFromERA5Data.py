import sys
from ReadERA5DataAndWriteERF_IC import ReadERA5_3DData
from Download_ERA5Data import *

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 download_era5.py <input_filename>")
        sys.exit(1)

input_filename = sys.argv[1]
filename  = download_era5_data(input_filename);
print("Filename is ", filename);

file_paths = [filename]

# Loop through each file and process it
count = 0;
is_IC = True
for file_path in file_paths:
    print(f"Processing file: {file_path}")
    #ReadERA5_SurfaceData("Laura_SurfaceData.grib", is_IC)
    ReadERA5_3DData(file_path, is_IC)
    is_IC = False
