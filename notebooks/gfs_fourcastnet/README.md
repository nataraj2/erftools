# Read weather data from FourCastNetGFS and write forecast files for ERF

This directory contains the python scripts to create forecast files for time-dependent boundary forcing  
 for hurricane simulations in ERF from the forecast weather data obtained from inference using the 
AI model FourCastNetGFS.

1. Give the year, month, day, time and the geographical area to download the data in a text file.  
For example, to download the data on August 26, 2020 at 00:00 (24 hour format)
```
year: 2020
month: 08
day: 26
time: 00:00
area: 50,-130,10,-50
```
Note: The geographical area is specified as latitude maximum, longitude minimum, latitude minimum, longitude maximum.

2. `python3 WriteICFromGFSData_FourCastNetGFS.py <input_file> <input_dir_with_GFS_forecast_data_grib2>`
The `input_file` is the text file in step 1. The `input_dir_with_GFS_forecast_data_grib2` is the directory which contains the forecast   
data from the FourCastNetGFS inference. If any packages are missing, install them using `pip install <package>`

3. The output VTK files for visualization is written into a directory `Output`. The forecast binary files (`*bin`) for ERF   
is also written into `Output`.

## Examples

Example inputs are given in the input file `input_for_Laura` and `input_for_Henri`. 

1. Run `python3 WriteICFromGFSData_FourCastNetGFS.py input_for_Laura`<input_dir>`
2. Visualize the VTK files in the `Output` directory in VisIt or ParaView.

