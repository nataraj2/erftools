import xarray as xr

def stagger_field(ds,field,surfval=0.0):
    """Interpolate cell-centered field to staggered z locations on the
    interior; extrapolate to domain top; and set the surface value.
    """
    zbot = ds.coords['zstag'].values[0]
    ztop = ds.coords['zstag'].values[-1]

    top = 1.5*ds[field].isel(z=-1) - 0.5*ds[field].isel(z=-2)
    top = top.expand_dims({'z':[ztop]})

    interior = xr.DataArray(
        0.5*(  ds[field].isel(z=slice(0,  -1)).values
             + ds[field].isel(z=slice(1,None)).values),
        coords={'t': ds.coords['t'],
                'z': ds.coords['zstag'].values[1:-1]},
        dims=('t','z')
    )

    surface = surfval * xr.ones_like(top)
    surface = surface.assign_coords(z=[zbot])

    da = xr.concat((surface,interior,top), dim='z')
    da = da.transpose('t','z')
    da = da.rename(z='zstag')

    return da
