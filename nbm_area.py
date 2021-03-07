import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Load dataset
now  = datetime.utcnow()-timedelta(hours=1)
data = xr.open_dataset(f'http://nomads.ncep.noaa.gov:80/dods/blend/blend{now:%Y%m%d}/blend_1hr_{str(now.hour).zfill(2)}z')

# Extract specifc variables
data = data[['tmp2m', 'dpt2m', 'wdir10m', 'wind10m', 'gust10m', 'apcpsfc', 'apcp254gtsfc', 'tstmsfc', 'tcdcsfc', 'asnow1016gtsfc', 'asnowsfc']]

# Slice dataset to bounds
data = data.sel(lat=slice(36.25, 36.54), lon=slice(-83.10, -82.65))

# Get mean over bounds
data = data.mean(dim=['lat', 'lon'])
data

# Setup Plot 
fig, axs = plt.subplots(ncols=1, nrows=6, constrained_layout=True, figsize=(15, 10))

# Add title
#plt.suptitle(f'weather.carterhumphreys.com', fontsize=16)

# Change datetiem64 to datetime
valid = datetime.utcfromtimestamp(data.time[0].values.astype('O')/1e9)

# Add plot headers
axs[0].set_title(f'NBM Meteogram for Hamblen County, TN', loc='left')
axs[0].set_title(f'Run: {valid.strftime("%a %Y-%m-%d %H:%M")} UTC', loc='right')

# 2m Temperature and Dewpoint
axs[0].plot(data.time, (data.tmp2m-273.15)*(9/5)+32, label='2m Temperature', color='red', linewidth=2, marker='o')
axs[0].plot(data.time, (data.dpt2m-273.15)*(9/5)+32, label='2m Dewpoint', color='green', linewidth=2, marker='o')

# Wind Speed and Wind Gust
axs[1].bar(data.time, data.wind10m*2.237, label='10m Wind Speed', color='royalblue', width=0.04)
axs[1].bar(data.time, (data.gust10m-data.wind10m)*2.237, bottom=data.wind10m, label='10m Wind Gust', color='navy', width=0.04)
axs[1].set_ylim(0)

# Wind Speed and Wind Gust
axs[2].scatter(data.time, data.wdir10m, label='10m Wind Direction', color='orange')
axs[2].set_yticks([0, 90, 180, 270, 360])

# Prob Snow, Prob Rain, Prob Thunder, Cloud Cover
axs[3].plot(data.time, data.tcdcsfc, label='Cloud Cover', color='gray', linewidth=2, marker='o')
axs[3].plot(data.time, data.apcp254gtsfc, label='Prob Precipitation', color='green', linewidth=2, marker='o')
axs[3].plot(data.time, data.tstmsfc, label='Prob Thunder', color='red', linewidth=2, marker='o')
axs[3].plot(data.time, data.asnow1016gtsfc, label='Prob Snow', color='blue', linewidth=2, marker='o')
axs[3].set_ylim(0, 100)

# Precipitation
axs[4].bar(data.time, np.around(data.apcpsfc/25.4, 2), label='1hr Precip', color='green', width=0.04)
axs[4].set_ylim(0)

# Snwofall
axs[5].bar(data.time, np.around(data.asnowsfc/39.37, 1), label='1hr Snowfall', color='grey', width=0.04)
axs[5].set_ylim(0)

for ax in axs:
    ax.grid(True)
    ax.legend(loc='upper left')
    ax.set_xlim(data.time[0], data.time[-1])

plt.show()
