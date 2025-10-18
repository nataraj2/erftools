# Read weather data from GFS and write initial condition file for ERF

This directory contains the python scripts to create an initial condition for hurricane simulations in ERF from the GFS weather data.

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

2. `srun -n 32 python3 WriteICFromGFSData.py <input_file> --do_forecast=true`
The `input_file` is the text file in step 1. If any packages are missing, install them using `pip install <package>`. This will download   
forecast data for 72 hours with an interval of 3 hours.

3. The output VTK files for visualization is written into a directory `Output`. The initial condition binary file (`*bin`) for ERF   
is also written into `Output`.

## Examples

Example inputs are given in the input file `input_for_Laura` and `input_for_Henri`. 

1. Run `python3 WriteICFromGFSData.py input_for_Laura --do_forecast=true`.  
2. Visualize the VTK files in the `Output` directory in VisIt or ParaView.

