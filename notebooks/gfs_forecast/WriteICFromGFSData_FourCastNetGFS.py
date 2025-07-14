import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import Download_GFS_Data
from erftools.preprocessing import ReadGFS_3DData_FourCastNetGFS

def get_area(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("area:"):
                area_str = line.split("area:")[1].strip()
                area = [float(x) for x in area_str.split(",")]
                return area
    raise ValueError("Area not found in the file.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 WriteICFromGFSData.py <input_file> <input_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    input_dir  = sys.argv[2]

    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory")
        sys.exit(1)

    is_IC = True
    area = get_area(input_file)    

    # Sort files in lexicographic order
    file_list = sorted(f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)))

    for filename in file_list:
        full_path = os.path.join(input_dir, filename)
        print(f"Processing file: {filename}")
        ReadGFS_3DData_FourCastNetGFS(full_path, area, is_IC)
