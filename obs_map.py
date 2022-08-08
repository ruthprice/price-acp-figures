# obs_map.py

# Creates a map of the Arctic with ASCOS, ATom
# and AO2018 routes plotted

# =====================================================================
import campaign_routes as cr
import numpy.ma as ma
import numpy as np
import matplotlib
import figure_props as fp
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature   # for land/sea shading
## =====================================================================
def plot_obs_map(atom_path, ascos_path, fig_filename):
    # Load AO2018 ship track
    ao2018_track = cr.get_ao2018_track()
    ao2018_lats = ao2018_track['lat']
    ao2018_lons = ao2018_track['lon']

    # Load ATom tracks and use ma to mask lats <60N
    atom1_t1, atom1_t2, atom1_lats, atom1_lons = cr.read_atom_nav_file(atom_path, 1)[:4]
    atom1_lats = ma.masked_where(atom1_lats <= 60, atom1_lats)
    atom1_lons = ma.masked_where(atom1_lats <= 60, atom1_lons)

    # Load ASCOS track
    ascos_flight_lons, ascos_flight_lats = cr.get_ascos_track(ascos_path)
    ascos_lons = []
    ascos_lats = []
    for i in np.arange(len(ascos_flight_lons)):
        ascos_lons.extend(ascos_flight_lons[i])
        ascos_lats.extend(ascos_flight_lats[i])

    # Make plot
    # plot on map
    cmap = matplotlib.cm.get_cmap('tab10')
    fig = plt.figure(figsize=(8.3*fp.cm,8.3*fp.cm),dpi=fp.dpi)
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,60,90],crs=ccrs.PlateCarree())
    ax.gridlines(linewidth=0.5)
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '110m',
                                            edgecolor='k',
                                            linewidth=0.5,
                                            facecolor='lightgrey')
    ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                            facecolor='lightblue')
    ax.add_feature(ocean_50m)
    ax.add_feature(land_50m)

    plt.plot(ao2018_lons, ao2018_lats, color=cmap(0),
             transform=ccrs.PlateCarree(), label='AO2018 ship',
             linestyle='None', marker='o', ms=1)
    plt.plot(atom1_lons, atom1_lats, color=cmap(1),
             transform=ccrs.PlateCarree(), label='ATom 1 flights',
             linestyle='None', marker='o', ms=1)
    plt.plot(ascos_lons, ascos_lats, color=cmap(4),
             transform=ccrs.PlateCarree(), label='ASCOS helicopter flights',
             linestyle='None', marker='o', ms=1)

    plt.legend(fontsize=6,facecolor='white',framealpha=1)
    plt.savefig(fig_filename, bbox_inches='tight', facecolor='white', format='pdf')
    plt.close()

    return