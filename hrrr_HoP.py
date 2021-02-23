import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from netCDF4 import num2date
import numpy as np
import xarray as xr
from datetime import datetime
import datetime as dt
from xarray.backends import NetCDF4DataStore
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import gaussian_filter
import metpy.calc as mpcalc
import numpy.ma as ma
from metpy.units import units
import scipy.ndimage as ndimage
from metpy.plots import USCOUNTIES
import matplotlib.patches as mpatches
import supplementary_tools as spt

mdate = spt.get_init_time('HRRR')[0]
init_hour = spt.get_init_time('HRRR')[1]
url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+mdate+'/hrrr_sfc.t'+init_hour+'z'
#url='http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr20201231/hrrr_sfc.t00z'
print(url)

# Create new directory
output_dir = mdate+'_'+init_hour+'00'
#output_dir = '20201231_0000'
spt.mkdir_p(output_dir)
spt.mkdir_p(output_dir+'/HRRR_ex_HoP')
#Parse data using MetPy
ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(25,55,0.25)
lons = np.arange(260,310,0.25)
total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402

#Initialize ptype arrays by grabbing first hour of categorical precip
catrain = ds['crainsfc'].squeeze().isel(time=1).squeeze()
catsnow = ds['csnowsfc'].squeeze().isel(time=1).squeeze()
catsleet = ds['cicepsfc'].squeeze().isel(time=1).squeeze()
catice = ds['cfrzrsfc'].squeeze().isel(time=1).squeeze()

total_precip=ds['apcpsfc'].isel(time=1).squeeze()*.0393700787402

acc_rain = np.ma.masked_where(catrain==0,total_precip)
acc_sleet = np.ma.masked_where(catsleet==0,total_precip)
acc_ice = np.ma.masked_where(catice==0,total_precip)
acc_snow = np.ma.masked_where(catsnow==0,total_precip)

acc_rain = acc_rain.filled(0)
acc_sleet = acc_sleet.filled(0)
acc_ice = acc_ice.filled(0)
acc_snow = acc_snow.filled(0)

t2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15
td2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15

u10 = ds['ugrd10m'].isel(time=1).squeeze()*1.94384449
v10 = ds['vgrd10m'].isel(time=1).squeeze()*1.94384449
ws10 = ((u10**2)+(v10**2))**.5
acc_fram = spt.fram(acc_ice,spt.wet_bulb(t2mi,td2mi),ws10)
print("INITIALIZATION SUCCESSFUL")

for i in range(0,49):

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tcdcclm':'tcc',
        'tmpprs': 'temperature',
        'ugrd10m': 'u',
        'vgrd10m': 'v',
        'mslmamsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'capesfc':'cape',
        'gustsfc':'sfcgust',
        'hcdchcll':'high_cloud',
        'mcdcmcll':'mid_cloud',
        'lcdclcll':'low_cloud',
        'vissfc':'sfcvis',
        'hgt263_k':'hgt_m10c',
        'hgt253_k':'hgt_m20c',
        'ltngclm':'lightning',
        'sbt124toa':'simsat',
        'hgt0c':'0chgt',
        'asnowsfc':'snow'
    })

    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    zH5 = data['temperature'].squeeze()
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)

    t2m = data['sfc_temp'].squeeze()
    t2mc = t2m-273.15
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2mc = td2m-273.15
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    wb2mc = spt.wet_bulb(t2mc,td2mc)
    #wb2mc = (wb2m-32.)*(5./9.)

    swg = data['sfcgust'].squeeze()
    cloudcover = data['tcc'].squeeze()
    high_cloud = data['high_cloud'].squeeze()
    mid_cloud = data['mid_cloud'].squeeze()
    low_cloud = data['low_cloud'].squeeze()
    vis = data['sfcvis'].squeeze()*0.000621371
    reflectivity = data['radar'].squeeze()
    cape = data['cape'].squeeze()
    lightning=data['lightning'].squeeze()
    dgz_depth = data['hgt_m20c'].squeeze()-data['hgt_m10c'].squeeze()
    simsat = data['simsat'].squeeze()
    hgt0c = data['0chgt'].squeeze()*3.28084
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    total_precip = total_precip+hrly_precip
    total_snow = data['snow'].squeeze()
   

    rain = np.ma.masked_where(catrain==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(catsnow==0,reflectivity)

    qrain = np.ma.masked_where(catrain==0,hrly_precip)
    qsleet = np.ma.masked_where(catsleet==0,hrly_precip)
    qice = np.ma.masked_where(catice==0,hrly_precip)
    qsnow = np.ma.masked_where(catsnow==0,hrly_precip)

    #Generate running accumulation total arrays for each ptype
    acc_snow = acc_snow+qsnow.filled(0)
    acc_sleet = acc_sleet+qsleet.filled(0)
    acc_ice = acc_ice+qice.filled(0)
    acc_rain = acc_rain+qrain.filled(0)

    mslpc = data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)
    smoothedsnow = ndimage.gaussian_filter(snow.filled(0),sigma=1,order=0)
    wind_slice = slice(36,-36,36)
    wind_slice_ne = slice(18,-18,18)
    wind_slice_me = slice(9,-9,9)
    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()
    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    blank = np.ones(np.shape(swg))

    blizzard = np.ma.masked_where((smoothedsnow>10)&(swg>30)|(swg==30)&(vis<0.25)|(vis==0.25),blank)
    blizzard = blizzard.filled(3)
    blizzard = np.ma.masked_where((smoothedsnow>10)&(swg>20)&(swg<30)&(vis<0.5)*(vis>0.25),blizzard)
    blizzard = blizzard.filled(2)

    print(np.max(blizzard))

    ###COMPUTE ICE ACCRETION###

    fram_accretion=spt.fram(qice,wb2mc,wspd)
    fram_accretion=fram_accretion.filled(0)
    acc_fram = acc_fram+fram_accretion

    ###GRAB LOCAL DATA###
    stations=['HOP','M21','2IO','CEY','PAH','PHT','UCY','MKL']
    coords=[[36.644940, -87.49936],[37.247911, -87.176467],[37.372355, -87.453684],[36.658890, -88.381034],
            [37.099267, -88.786733],[36.377577, -88.414566],[36.371405, -88.989613],[35.619656, -88.924166],
           ]


    station_qpf = []
    station_temp = []
    print(len(stations))
    print(len(coords))
    for i in range(len(stations)):
        station = stations[i]
        lat = coords[i][0]
        lon = coords[i][1]
        qpf = total_precip.interp(lat=lat,lon=lon).values
        temp = t2m.interp(lat=lat,lon=lon).values
        station_qpf.append(np.round(qpf,1))
        station_temp.append(np.round(temp,1))

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(44,15))

    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    gs.update(hspace=0.01,wspace=0.01)
    ax1 = fig.add_subplot(gs[:, 1], projection = zH5_crs)
    ax2 = fig.add_subplot(gs[0, 0], projection = zH5_crs)
    ax3 = fig.add_subplot(gs[1, 0], projection = zH5_crs)
    ax4 = fig.add_subplot(gs[0, 2], projection = zH5_crs)
    ax5 = fig.add_subplot(gs[1, 2], projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))

    ax2.coastlines(resolution='10m')
    ax2.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax2.add_feature(cfeature.STATES.with_scale('10m'))

    ax3.coastlines(resolution='10m')
    ax3.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax3.add_feature(cfeature.STATES.with_scale('10m'))

    ax4.coastlines(resolution='10m')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax4.add_feature(cfeature.STATES.with_scale('10m'))

    ax5.coastlines(resolution='10m')
    ax5.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax5.add_feature(cfeature.STATES.with_scale('10m'))

    #fig.suptitle("HRRR Forecast valid at " + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=36)

    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    print(dtfs)

########## PLOTTING #######################################################
    tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.7)
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=True, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]

    try:
        ra = ax1.contourf(x,y,rain,cmap='Greens',levels=ref_levs,alpha=0.7)
    except:
        print('no rain')
    try:
        sn = ax1.contourf(x,y,snow,cmap='cool',levels=ref_levs,alpha=0.7)
    except:
        print('no snow')
    try:
        ip = ax1.contourf(x,y,sleet,cmap='autumn',levels=ref_levs,alpha=0.7)
    except:
        print('no sleet')
    try:
        zr = ax1.contourf(x,y,ice, cmap='RdPu',levels=ref_levs,alpha=0.7)
    except:
        print('no ice')

    #cax1 = fig.add_axes([.35,.5,.02,.25])
    #cax2 = fig.add_axes([.35,.75,.02,.25])

    refp = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs) #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    capep = ax4.contourf(x, y, cape, levels=[25,50,75,100,150, 250, 500, 750, 1000], extend='max', alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    lgt = ax4.contour(x,y,lightning,levels=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
    cb = fig.colorbar(capep, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000])
    cb.set_label('CAPE (J/kg)', size='large')

    #ax4.barbs(x[ws2],y[ws2],s5u[ws2,ws2],s5v[ws2,ws2], length = 7)
    blue = mpatches.Patch(color='#6ec6fd', label='Snow')
    orange = mpatches.Patch(color='#e3aa52', label='Sleet')
    pink = mpatches.Patch(color='#ff94aa', label='Freezing Rain')
    green = mpatches.Patch(color='#9cc78d', label='Rain')
    leg = ax1.legend(handles=[blue,orange,pink,green],loc=4,title='Precipitation Type',framealpha=1)
    leg.set_zorder(100)

    swgc = ax5.contourf(x,y,swg, levels=range(0,60,5), cmap = 'hot')
    sw3g = ax5.contour(x,y,swg,levels=[30],colors=['lime'],linewidths=1.5)
    visc = ax5.contour(x,y,vis,levels=[0.25],colors=['violet'],linewidths=1.5)
    try:
        sn = ax5.contour(x,y,smoothedsnow,cmap='cool',levels=[15],linewidths=1.5)
    except:
        print('no snow')
    cbar4 = fig.colorbar(swgc, orientation='vertical',pad=0.01,ax=ax5, aspect = 50, extendrect=False, ticks = np.arange(0,70,10))
    blue = mpatches.Patch(color='dodgerblue',label='Snow')
    lime = mpatches.Patch(color='lime', label='35 mph Gusts')
    pink = mpatches.Patch(color='violet', label='1/4mi Visibility')
    leg = ax5.legend(handles=[blue,lime,pink],loc=4,title='Blizzard Ingredients',framealpha=1)
    cbar4.set_label('Surface Wind Gust (kts)')

    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)
    ax1.set_title('Precip Type, 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    cloud_levs = [40,45,50,55,60,65,70,75,80,85,90,95]
    #high_cols = ['#33c7f4f4','#40b6f0f1','#4da5edee','#5990e8ea','#6678ebed','#7359e6e8','#8030e5e8','#8c18d5d8','#9916c2c5','#9914b5b8','#9912a5a5','#99119579','#990e7a7c','#990d7173','#990c686a','#990b6365']
    #mid_cols = ['#4dbcfac1','#59aef9b4','#66a0f8a7','#7391f79a','#8083f68d','#8c75f57f','#9966f472','#a658f365','#b34af258','#b33bf14b','#b32df03d','#b325e935','#b31ce32d','#b31ad52a','#b319c827','#b317ba25']
    #low_cols = ['#66ea99f4','#73e78cf3','#80e47ef1','#8ce170f0','#99db54ed','#a6d846ec','#b3d539ea','#bfd32be9','#ccd01de7','#ccc617de','#ccba16d0','#ccae14c2','#cca113b4','#cc9511a6','#cc881098','#cc7c0f8a']

    high_cols = ['#c7f4f4','#b6f0f1','#a5edee','#90e8ea','#78ebed','#59e6e8','#30e5e8','#18d5d8','#16c2c5','#14b5b8','#12a5a5','#119579']
    mid_cols = ['#bcfac1','#aef9b4','#a0f8a7','#91f79a','#83f68d','#75f57f','#66f472','#58f365','#4af258','#3bf14b','#2df03d','#25e935']
    low_cols = ['#ea99f4','#e78cf3','#e47ef1','#e170f0','#db54ed','#d846ec','#d539ea','#d32be9','#d01de7','#c617de','#ba16d0','#ae14c2']

    blue = mpatches.Patch(color='#119579', label='High Clouds')
    green = mpatches.Patch(color='#25e935', label='Mid-Level Clouds')
    purple = mpatches.Patch(color='#ae14c2',label='Low-Level Clouds')
    leg = ax2.legend(handles=[blue,green,purple],loc=4,framealpha=1)
    leg.set_zorder(100)
    tccp = ax2.contourf(x,y,cloudcover,cmap='Greys',levels=cloud_levs,alpha=0,extend='max')
    lccp = ax2.contourf(x,y,low_cloud, colors=low_cols,levels=cloud_levs,alpha = 0.35,extend='max')
    mccp = ax2.contourf(x,y,mid_cloud, colors=mid_cols,levels=cloud_levs,alpha = 0.25,extend='max')
    hccp = ax2.contourf(x,y,high_cloud, colors=high_cols,levels=cloud_levs,alpha = 0.15,extend='max')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(tccp,orientation='vertical',pad=0.01,ax=ax2,aspect=50,extendrect=False, ticks=np.arange(10,100,10))
    cbar2.set_label('Cloud Cover (%)',fontsize=14)

    tprecip = ax3.contourf(x,y,total_precip, alpha = 0.7, cmap = 'cool',transform=zH5_crs, levels=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    tcprecip = ax3.contour(x,y,total_precip,colors=['b','darkblue','darkviolet'],levels=[0.5,1,1.5],linewidths=1.5)
    cbar3 = fig.colorbar(tprecip,orientation='vertical',pad=0.01,ax=ax3,aspect=50,extendrect=False,ticks=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    cbar3.set_label('Total Precipitation (inches)',fontsize=14)


    #refp3 = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs)
    #cbar4 = fig.colorbar(tccp,orientation='horizontal',pad=0.01,ax=ax4,aspect=50,extendrect=False)

    #ax1.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax2.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax3.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax4.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax5.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #north=50, south=15, east=-70, west=-115

    sub_w1 = 260
    sub_w = 262
    sub_e = 295
    sub_n = 50
    sub_s = 25

    ax1.set_extent((sub_w1-1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    #fig.canvas.draw()
    fig.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex_HoP/EC_fivepanelwinter9_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    ax1.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=6)
    ax1.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    for i in range(len(stations)):
        station = stations[i]
        lat = coords[i][0]
        lon = coords[i][1]
        ax3.text(360+lon,lat,str(station_qpf[i]),ha='center',transform=zH5_crs)
        ax1.text(360+lon,lat,str(station_temp[i]),ha='center',transform=zH5_crs)
    plt.savefig(output_dir+'/HRRR_ex_HoP/NE_fivepanelwinter9_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()

    plt.clf()

### SECOND PLOT ###
    fig2 = plt.figure(figsize=(15,15))
    ax6 = fig2.add_subplot(111, projection = zH5_crs)

    ax6.coastlines(resolution='10m')
    ax6.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax6.add_feature(cfeature.STATES.with_scale('10m'))

    tmp_2m = ax6.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax6.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr6 = fig2.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax6, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5))
    cbr6.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax6.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]

    try:
        ra = ax6.contourf(x,y,rain,cmap='Greens',levels=ref_levs,alpha=0.7)
    except:
        print('no rain')
    try:
        sn = ax6.contourf(x,y,snow,cmap='cool',levels=ref_levs,alpha=0.7)
    except:
        print('no snow')
    try:
        ip = ax6.contourf(x,y,sleet,cmap='autumn',levels=ref_levs,alpha=0.7)
    except:
        print('no sleet')
    try:
        zr = ax6.contourf(x,y,ice, cmap='RdPu',levels=ref_levs,alpha=0.7)
    except:
        print('no ice')

    ax6.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)
    ax6.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.set_title('HRRR Composite Forecast valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/HRRR_ex_HoP/EC_ptype_composite_'+dtfs+'_.png',bbox_inches='tight')
    ax6.set_extent((265, 277, 30, 40))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=6)
    plt.savefig(output_dir+'/HRRR_ex_HoP/NE_ptype_composite_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    ax6.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    wsl = slice(5,-5,5)
    ax6.barbs(x[wsl],y[wsl],u_10m[wsl,wsl],v_10m[wsl,wsl], length=6)
    ax6.set_extent((270,274,35,39))
    
    ax6.text(360+lon,lat,str(station_temp[i]),ha='center',transform=zH5_crs)
    
    plt.savefig(output_dir+'/HRRR_ex_HoP/local_ptype_composite_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.close()
    plt.clf()
    plt.close()
    plt.clf()
    ### END SECOND PLOT ###