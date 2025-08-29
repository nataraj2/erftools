# Read weather data from ERA5 and write initial condition file for ERF

This directory contains the python scripts to create an initial condition for hurricane simulations in ERF from the ERA5 weather data.

1. Follow the steps here https://cds.climate.copernicus.eu/how-to-api to create a free account   
   to download ERA5 data

2. Give the year, month, day, time and the geographical area to download the data in a text file.  
For example, to download the data on August 26, 2020 at 00:00 (24 hour format)
```
year: 2020
month: 08
day: 26
time: 00:00
area: 50,-130,10,-50
```
Note: The geographical area is specified as latitude maximum, longitude minimum, latitude minimum, longitude maximum.

3. Run the script.    
`srun -n 32 python3 WriteICFromERA5Data.py input_for_Henri --do_forecast=true --forecast_time_hours=72 --interval_hours=3`  
This uses 32 MPI ranks to download and process the weather data for hurricane Henri for a total of 72 hours with an interval of 3 hours.   

4. The output VTK files for visualization is written into  `Output/VTK/3D` for 3D data and   
`Output/VTK/Surface` for surface data.  

5. The following directories are to be copied into the ERF run directory.
`Output/3D`- The binary files (`*bin`) for lateral forcing.   
`Output/Surface` -  The binary files (`*.bin`) for surface fluxes.

## Examples

Various example inputs for different hurricanes are provided in this folder. For eg. for a 3-day hindcasting  
for Hurricane Henri

1. Run `srun -n 32 python3 WriteICFromERA5Data.py input_for_Henri --do_forecast=true --forecast_time_hours=72 --interval_hours=3`
2. Visualize the VTK files in the `Output/VTK` directory in VisIt or ParaView.

