# ERFtools
A collection of scripts for facilitating the usage of ERF.

## Examples

### Converting a WRF namelist into ERF inputs
```python
from erftools.preprocessing import WRFInputDeck
wrf = WRFInputDeck('namelist.input')
wrf.process_initial_conditions('wrfinput_d01',
                               landuse_table_path='/Users/equon/WRF/run/LANDUSE.TBL',
                               write_hgt='terrain_height.txt',
                               write_z0='roughness_height.txt')
wrf.write_inputfile('inputs')
```

### Postprocessing data logs
Data logs are output with the `erf.data_log` param and can include time histories of surface conditions and planar averaged profiles (e.g., for idealized LES simulations)
```python
from erftools.postprocessing import DataLog
log = DataLog(f'{simdir}/surf_hist.dat',
              f'{simdir}/mean_profiles.dat',
              f'{simdir}/flux_profiles.dat',
              f'{simdir}/sfs_profiles.dat')
log.calc_stress()
log.est_abl_height('max_theta_grad')
print(log.ds) # data are stored in an xarray dataset
```

### Reading ERA5 weather data and writing VTK output and initial condition file for ERF
To create an initial condition for ERF from the ERA5 weather data and visualize it

1. Follow the steps here https://cds.climate.copernicus.eu/how-to-api to create a free account
   to download ERA5 data

2. Give the year, month, day, time and the geographical area to download the data in a input text file.
For example, to download the data on August 26, 2020 at 00:00 (24 hour format)
```
year: 2020
month: 08
day: 26
time: 00:00
area: 50,-130,10,-50
```
Note: The geographical area is specified as latitude maximum, longitude minimum, latitude minimum, longitude maximum.

3. Execute the following script from the `notebooks/era5/` folder as 
```
python3 <script_name> <input_file>
```
where the `input_file` is the text file in step 2 above.  
If any packages are missing, install them using `pip install <package>`
```python
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import Download_ERA5_Data
from erftools.preprocessing import ReadERA5_3DData

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 WriteICFromERA5Data.py <input_filename>")
        sys.exit(1)

input_filename = sys.argv[1]
filename = Download_ERA5_Data(input_filename);
print("Filename is ", filename);

is_IC = True
print(f"Processing file: {filename}")
ReadERA5_3DData(filename, is_IC)
```

4. The output VTK files for visualization is written into a directory `Output`. The initial condition binary file (`*bin`) for ERF
is also written into `Output`.

## Examples

Example inputs are given in the input files in `notebooks/era5` -- `input_for_era5_Laura` and `input_for_era5_Henri`.

1. Run `python3 WriteICFromERA5Data.py input_for_era5_Laura`.
2. Visualize the VTK files in the `Output` directory in VisIt or ParaView.


