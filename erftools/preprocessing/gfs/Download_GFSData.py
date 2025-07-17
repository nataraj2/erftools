# --- Read user values from file ---
def read_user_input(filename):
	inputs = {}
	with open(filename, 'r') as f:
		for line in f:
			if ':' not in line:
				continue
			key, value = line.strip().split(':', 1)
			key = key.strip()
			value = value.strip()
			if key == "area":
				inputs[key] = [float(x) for x in value.split(',')]
			else:
				inputs[key] = [value]
	return inputs

from datetime import datetime, timedelta

def generate_urls_for_24_hours(data):
	year = data['year'][0]
	month = data['month'][0].zfill(2)
	day = data['day'][0].zfill(2)
	hour = data['time'][0].split(':')[0].zfill(2)

	yyyymmddhh = f"{year}{month}{day}{hour}"
	yyyymmdd = f"{year}{month}{day}"

	urls = []
	filenames = []
	# 24 hours / 3 = 8 forecast times: 0,3,6,...,21
	for fhour in range(0, 24, 3):
		fhour_str = f"f{fhour:03d}"
		filename = f"gfs.0p25.{yyyymmddhh}.{fhour_str}.grib2"
		url = f"https://data-osdf.rda.ucar.edu/ncar/rda/d084001/{year}/{yyyymmdd}/{filename}"
		urls.append(url)
		filenames.append(filename)
	return urls, filenames 

def construct_url_filename(data):

	year = data['year'][0]		# your read_user_input returns a list for year
	month = data['month'][0].zfill(2)
	day = data['day'][0].zfill(2)
	hour = data['time'][0].split(':')[0].zfill(2)

	yyyymmddhh = f"{year}{month}{day}{hour}"
	yyyymm = f"{year}{month}"
	yyyymmdd = f"{year}{month}{day}"

	# Historical forecast data
	filename = f"gfs.0p25.{yyyymmddhh}.f000.grib2"
	url = f"https://data-osdf.rda.ucar.edu/ncar/rda/d084001/{year}/{yyyymmdd}/{filename}"

	# Final reanalaysis data
	#filename = f"gdas1.fnl0p25.{yyyymmddhh}.f00.grib2"
	#url = f"https://data-osdf.rda.ucar.edu/ncar/rda/d083003/{year}/{yyyymm}/{filename}"
	
	print("URL is ", url)

	return url, filename

def Download_GFS_Data(inputs):
	data = read_user_input(inputs)
	lat_max, lon_min, lat_min, lon_max = data.get('area')

	url, filename = construct_url_filename(data)

	print("Download URL:", url)
	print("Filename:", filename)

	# Optional: actually download
	import urllib.request
	urllib.request.urlretrieve(url, filename)
	print("Download complete.")

	area = [lat_max, lon_min, lat_min, lon_max]
	return filename, area

def Download_GFS_ForecastData(inputs):
	data = read_user_input(inputs)
	lat_max, lon_min, lat_min, lon_max = data.get('area')

	urls, filenames = generate_urls_for_24_hours(data)

	 # Optional: actually download
	import urllib.request

	for url, filename in zip(urls, filenames):
		print(f"Downloading {url} -> {filename}")
		urllib.request.urlretrieve(url, filename)

	area = [lat_max, lon_min, lat_min, lon_max]
	return filenames, area


