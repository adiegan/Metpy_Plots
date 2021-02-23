#Trying to plot some basic 700mb data

#Importing relevant libraries

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
#plt.rcParams({'figure.max_open_warning':0})
from netCDF4 import num2date
import numpy as np
import xarray as xr
from siphon.catalog import TDSCatalog
from datetime import datetime
import datetime as dt
from xarray.backends import NetCDF4DataStore
import metpy.constants as mpconstants

#import metpy
# Any import of metpy will activate the accessors
import metpy.calc as mpcalc
#from metpy.testing import get_test_data
from metpy.units import units
import metpy.constants as mpconstants


def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

#Get data using siphon
best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
best_ds = best_gfs.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric')

data = ncss.get_data(query)



#Parse data using MetPy
ds = xr.open_dataset(NetCDF4DataStore(data))
data = ds.metpy.parse_cf()

filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
coord_names = list(filtered_ds.coords)
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/QG')
    

for i in range(0,41):
    #Get data using siphon
    best_gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/CONUS_20km/catalog.html?dataset=grib/NCEP/GFS/CONUS_20km/Best')
    best_ds = best_gfs.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric')
    
    data = ncss.get_data(query)
    print(list(data.variables))
    ds = xr.open_dataset(xr.backends.NetCDF4DataStore(data)).metpy.parse_cf()
    ds
    
    vtime = ds['Temperature_isobaric'].metpy.time[0]
    
    #lats = ds['lat'].metpy.unit_array
    #lons = ds['lon'].metpy.unit_array
    
    lat_var = data.variables['y']
    lon_var = data.variables['x']
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lon_var, lat_var)
    
    #dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    hght_700 = ds['Geopotential_height_isobaric'].metpy.sel(vertical=700 * units.hPa,time=vtime)
    tmpk_700 = ds['Temperature_isobaric'].metpy.sel(vertical=700 * units.hPa,time=vtime)

    # 700 hPa u-component_of_wind
    uwnd_700 = ds['u-component_of_wind_isobaric'].metpy.sel(vertical=700 * units.hPa,time=vtime)

    # 700 hPa v-component_of_wind
    vwnd_700 = ds['v-component_of_wind_isobaric'].metpy.sel(vertical=700 * units.hPa,time=vtime)
    # Cell content replaced by load magic replacement.
    # Remaining variables needed to compute QG Omega forcing terms
    hght_500 = ds.Geopotential_height_isobaric.metpy.sel(vertical=500 * units.hPa,
                                                        time=vtime)
    uwnd_500 = ds['u-component_of_wind_isobaric'].metpy.sel(vertical=500 * units.hPa,
                                                        time=vtime)
    vwnd_500 = ds['v-component_of_wind_isobaric'].metpy.sel(vertical=500 * units.hPa,
                                                        time=vtime)
    uwnd_900 = ds['u-component_of_wind_isobaric'].metpy.sel(vertical=900 * units.hPa,
                                                        time=vtime)
    vwnd_900 = ds['v-component_of_wind_isobaric'].metpy.sel(vertical=900 * units.hPa,
                                                        time=vtime)
                                                        
    # Set constant values that will be needed in computations

    # Set default static stability value
    sigma = 2.0e-6 * units('m^2 Pa^-2 s^-2')

    # Set f-plane at typical synoptic f0 value
    f0 = 1e-4 * units('s^-1')

    # Use dry gas constant from MetPy constants
    Rd = mpconstants.Rd
    
    # Smooth Heights
    # For calculation purposes we want to smooth our variables
    # a little to get to the "synoptic values" from higher
    # resolution datasets

    # Number of repetitions of smoothing function
    n_reps = 50

    # Apply the 9-point smoother
    hght_700s = mpcalc.smooth_n_point(hght_700, 9, n_reps)#.metpy.unit_array
    hght_500s = mpcalc.smooth_n_point(hght_500, 9, n_reps)#.metpy.unit_array

    tmpk_700s = mpcalc.smooth_n_point(tmpk_700, 9, n_reps)#.metpy.unit_array
    tmpc_700s = tmpk_700s.to('degC')

    uwnd_700s = mpcalc.smooth_n_point(uwnd_700, 9, n_reps)#.metpy.unit_array
    vwnd_700s = mpcalc.smooth_n_point(vwnd_700, 9, n_reps)#.metpy.unit_array

    uwnd_500s = mpcalc.smooth_n_point(uwnd_500, 9, n_reps)#.metpy.unit_array
    vwnd_500s = mpcalc.smooth_n_point(vwnd_500, 9, n_reps)#.metpy.unit_array

    uwnd_900s = mpcalc.smooth_n_point(uwnd_900, 9, n_reps)#.metpy.unit_array
    vwnd_900s = mpcalc.smooth_n_point(vwnd_900, 9, n_reps)#.metpy.unit_array
    
    # Absolute Vorticity Calculation
    avor_900 = mpcalc.absolute_vorticity(uwnd_900s, vwnd_900s, dx, dy, lat_var)
    avor_500 = mpcalc.absolute_vorticity(uwnd_500s, vwnd_500s, dx, dy, lat_var)

    # Advection of Absolute Vorticity
    vortadv_900 = mpcalc.advection(avor_900, (uwnd_900s, vwnd_900s), (dx, dy)).to_base_units()
    vortadv_500 = mpcalc.advection(avor_500, (uwnd_500s, vwnd_500s), (dx, dy)).to_base_units()

    # Differential Vorticity Advection between two levels
    diff_avor = ((vortadv_900 - vortadv_500)/(400 * units.hPa)).to_base_units()

    # Calculation of final differential vorticity advection term
    term_A = (-f0 / sigma * diff_avor).to_base_units()
    print(term_A.units)
    
    # 700-hPa Temperature Advection
    tadv_700 = mpcalc.advection(tmpk_700s, (uwnd_700s, vwnd_700s), (dx, dy)).to_base_units()
    # Laplacian of Temperature Advection
    lap_tadv_700 = mpcalc.laplacian(tadv_700, deltas=(dy, dx))

    # Final term B calculation with constants
    term_B = (-Rd / (sigma * (700 * units.hPa)) * lap_tadv_700).to_base_units()
    print(term_B.units)
    
    # Set some contour intervals for various parameters

    # CINT 500 hPa Heights
    clev_hght_500 = np.arange(0, 7000, 60)
    # CINT 700 hPa Heights
    clev_hght_700 = np.arange(0, 7000, 30)
    # CINT 700 hPa Temps
    clev_tmpc_700 = np.arange(-40, 40, 5)
    # CINT Omega terms
    clev_omega = np.arange(-20, 21, 2)
    
    # Set some projections for our data (Plate Carree)
    # and output maps (Lambert Conformal)

    # Data projection; NARR Data is Earth Relative
    dataproj = ccrs.PlateCarree()

    # Plot projection
    # The look you want for the view, LambertConformal for mid-latitude view
    plotproj = ccrs.LambertConformal(central_longitude=-100.,
                                     central_latitude=40.,
                                     standard_parallels=[30, 60])
    # Cell content replaced by load magic replacement.
    fig=plt.figure(1, figsize=(15.,12.))

    # Upper-Left Panel
    ax=plt.subplot(111,projection=plotproj)
    ax.set_extent([-125.,-73,25.,50.],ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES,linewidth=0.5)

    # Contour #1
    cs = ax.contour(lons, lats, hght_700s, clev_hght_700,colors='k',
                linewidths=1.5, linestyles='solid', transform=dataproj)
    plt.clabel(cs, fontsize=10, inline=1, inline_spacing=3, fmt='%i',
                rightside_up=True, use_clabeltext=True)

    # Contour #2
    cs2 = ax.contour(lons, lats, tmpc_700s, clev_tmpc_700, colors='grey',
                linewidths=1.0, linestyles='dotted', transform=dataproj)
    plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=3, fmt='%d',
           rightside_up=True, use_clabeltext=True)

    # Colorfill
    cf = ax.contourf(lons, lats, (term_A+term_B)*10**12, clev_omega,
                 cmap=plt.cm.RdYlBu_r, extend='both', transform=dataproj)
    plt.colorbar(cf, orientation='horizontal', pad=0.0, aspect=50, extendrect=True)

# Vector
ax.barbs(lons.m, lats.m, uwnd_700s.to('kts').m, vwnd_700s.to('kts').m,
         regrid_shape=15, transform=dataproj)

# Titles
plt.title('700-hPa Geopotential Heights, Temperature (C),\n'
          'Winds (kt), and QG Omega Forcings ($*10^{12}$ kg m$^{-3}$ s$^{-3}$)',loc='left')
plt.title('VALID: ' + vtime_str, loc='right')

plt.show()                            