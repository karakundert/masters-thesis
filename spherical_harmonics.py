from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy.special as sp

long = -100
lat = 30

def plot(filename, m, n, longitude=0, latitude=0, inches=(1,1), 
         cmap='RdYlBu', points=100):

    figure, ax = plt.subplots(1,1)
    figure.set_size_inches(*inches)

    lon = np.linspace(0, 2*np.pi, points)
    lat = np.linspace(-np.pi / 2, np.pi / 2, points)
    colat = lat + np.pi / 2
    d = np.zeros((len(lon), len(colat)), dtype=np.complex64)

    meshed_grid = np.meshgrid(lon, lat)
    lat_grid = meshed_grid[1]
    lon_grid = meshed_grid[0]

    mp = Basemap(projection='ortho', lat_0=latitude, lon_0=longitude, ax=ax)
    mp.drawmapboundary()
    mp.drawmeridians(np.arange(0, 360, 30))
    mp.drawparallels(np.arange(-90, 90, 30))

    for j, yy in enumerate(colat):
        for i, xx in enumerate(lon):
            d[i,j] = sp.sph_harm(m, n, xx, yy)

    drm = np.round(np.real(d) / np.max(np.real(d)), 2)
    x, y = mp(np.degrees(lon_grid), np.degrees(lat_grid))
    mp.pcolor(x, y, np.transpose(drm), cmap=cmap)

    figure.savefig(filename, transparent=True)

i = 0
for l in xrange(4):
    for m in xrange(4):
        if m >= l:
            filename = 'sph%d.png' % i
            plot(filename, l, m, long, lat)
            i += 1
