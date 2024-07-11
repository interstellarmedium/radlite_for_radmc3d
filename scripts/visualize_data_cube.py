# ### interactive viewer of the radlite datacube
# ### based on https://matplotlib.org/stable/api/animation_api.html
# note that although this works out side of a jupyter notebook, it does not have the gui control

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits  
from astropy.visualization import (AsinhStretch, LogStretch, ImageNormalize)
from matplotlib.animation import FuncAnimation
#from IPython.display import HTML

path = "./radlite/"
filename = "radlite_image.fits"

f = fits.open(path+filename)
continuum = f[1].data
line = f[2].data
hd = f[2].header
wave = hd['WAVELEN']    # microns

v = hd['CRVAL3'] + (1 + np.arange(hd['NAXIS3']) - hd['CRPIX3']) * hd['CDELT3']
i1, i2 = 40, 60
j1, j2 = 40, 60
k1, k2 = 40, 60
frames = np.linspace(k1, k2, k2-k1+1, dtype=int)
norm = ImageNormalize(np.nanmax(line, axis=0), stretch=AsinhStretch())
slice = line[k1, j1:j2, i1:i2]

fig, ax = plt.subplots(figsize=(8,8))
implt = ax.imshow(slice, cmap='gnuplot2', origin='lower', interpolation='bilinear')
ax.set_xlabel("X")
ax.set_ylabel("Y")

def animate(k):
    implt.set_data(line[k, j1:j2, i1:i2])
    ax.text(0.03, 0.94, f"{v[k]:.2f} km/s", fontsize=14, color='white', backgroundcolor='black', transform = ax.transAxes)
    return implt,

ani = FuncAnimation(fig, animate, frames=frames)
#HTML(ani.to_jshtml())
plt.show()
