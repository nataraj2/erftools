import cartopy.feature as cfeature
import pyvista as pv
import numpy as np

# Define continental US bounds
lon_min, lon_max = -125, -66.5
lat_min, lat_max = 24, 50

# Get US state borders
states_feature = cfeature.STATES.with_scale('50m')
geoms = list(states_feature.geometries())

all_points = []
lines_flat = []
pt_id = 0

for geom in geoms:
    try:
        # MultiLineString
        for line in geom.boundary.geoms:
            coords = [(x, y) for x, y in line.coords if lon_min <= x <= lon_max and lat_min <= y <= lat_max]
            n = len(coords)
            if n < 2:
                continue  # skip very short lines
            all_points.extend([[x, y, 0] for x, y in coords])
            lines_flat.extend([n] + list(range(pt_id, pt_id+n)))
            pt_id += n
    except AttributeError:
        # single LineString
        coords = [(x, y) for x, y in geom.boundary.coords if lon_min <= x <= lon_max and lat_min <= y <= lat_max]
        n = len(coords)
        if n < 2:
            continue
        all_points.extend([[x, y, 0] for x, y in coords])
        lines_flat.extend([n] + list(range(pt_id, pt_id+n)))
        pt_id += n

# Convert to numpy arrays
all_points = np.array(all_points)
lines_flat = np.array(lines_flat, dtype=np.int64)

# Create PolyData
polydata = pv.PolyData()
polydata.points = all_points
polydata.lines = lines_flat

# Save to VTK
polydata.save("US_continental_state_borders.vtk")
print("Saved US_continental_state_borders.vtk")

