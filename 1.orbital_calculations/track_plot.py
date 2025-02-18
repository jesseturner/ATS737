import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#--- Read data into a dataframe
file_path = 'Track.txt'
df = pd.read_csv(file_path, sep=r'\s{2,}', skiprows=2, header=None, engine='python')
df.columns = [
    'Date & Time (UTC)', 'Zenith', 'Azimuth', 'Elevation'
]

#--- Filter to within 5 minutes of maximum elevation

#------ Convert date and time to datetime
df['Date & Time (UTC)'] = pd.to_datetime(df['Date & Time (UTC)'])

#------ Find row with max elevation
max_elevation_time = df.loc[df['Elevation'].idxmax(), 'Date & Time (UTC)']

#------ Define the time window
time_window_start = max_elevation_time - pd.Timedelta(minutes=5)
time_window_end = max_elevation_time + pd.Timedelta(minutes=5)

#------ Filter the dataframe to this window
filtered_df = df[(df['Date & Time (UTC)'] >= time_window_start) & 
                 (df['Date & Time (UTC)'] <= time_window_end)]

print(filtered_df.head())

#--- Plot the results
fig,ax = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

ax[0].scatter(filtered_df['Date & Time (UTC)'], filtered_df['Elevation'], color='#574AE2', s=30, label='Elevation')
ax[0].set_ylabel('Elevation (deg)')
ax[0].set_title('ISS Elevation during Overpass')

ax[1].scatter(filtered_df['Date & Time (UTC)'], filtered_df['Azimuth'], color='#90C3C8', s=30, label='Azimuth')
ax[1].set_ylabel('Azimuth (deg)')
ax[1].set_title('ISS Azimuth during Overpass')
ax[1].set_xlabel('Date & Time (UTC)')

fig.suptitle('ISS Sky Track 2024 Dec 3\nFrom Fort Collins, Colorado (40.5853 N, 105.0844 W)', fontsize=21)

fig.savefig("results/track_plot_Fort_Collins", dpi=200, bbox_inches='tight')
