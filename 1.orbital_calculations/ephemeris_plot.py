import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as feature
import matplotlib.pyplot as plt
import numpy as np

#--- Read data into a dataframe
file_path = 'EPHEMERIS.DAT'

#------ Skip first four lines of headers
df = pd.read_csv(file_path, sep='\s+', skiprows=4, header=None)

df.columns = [
    'Time (days)', 'Geocentric Latitude', 'Longitude', 'R (km)',
    'X (km)', 'Y (km)', 'Z (km)'
]

print(df.head())

#--- Plot the data
projection=ccrs.PlateCarree(central_longitude=0)
fig,ax=plt.subplots(1, figsize=(12,12),subplot_kw={'projection': projection})
ax.scatter(df['Longitude'], df['Geocentric Latitude'], color='#574AE2', s=30, transform=ccrs.PlateCarree())

#------ Set the title
ax.set_title('ISS Ground Track 2024 Dec 3 \n Centered on Colorado, US')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

#------ Option: Global or local
#ax.set_global()
#------ Colorado
ax.set_extent([-110, -100, 34, 43], crs=ccrs.PlateCarree())

#------ label for Fort Collins, CO
ax.text(-105.0844, 40.5853, 'Fort Collins', color='black', fontsize=12,
        ha='center', transform=ccrs.PlateCarree())

#ax.gridlines(draw_labels=True)
ax.coastlines(resolution='50m', color='black', linewidth=1)
ax.add_feature(feature.STATES, zorder=100, edgecolor='#000', facecolor='none', linewidth=0.5)

#------ Set the file name
fig.savefig("results/ephemeris_plot_CO", dpi=200, bbox_inches='tight')