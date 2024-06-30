### calculate the SA signal for a given slit PA
### based on https://learn.astropy.org/tutorials/PVDiagramPlotting.html
### returns velocity, mean, and standard deviation of SA signal

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
from astropy import units as u
from astropy import wcs
from astropy.io import fits  

from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice, Path
from astropy.visualization import (AsinhStretch, LogStretch, ImageNormalize)

def measure(path='./radlite/', filename='radlite_image.fits', slitPA=0, plot=True):
    # slit position angle in degrees, E of N
    if slitPA < 0 or slitPA > 360:
        print('slitPA must be between 0 and 360')
        sys.exit()

    # slit width (hardwired in for IRTF)
    slit_width = 0.3*u.arcsec

    f = fits.open(path+filename)
    continuum = f[1].data
    line = f[2].data
    hd = f[2].header
    wave = hd['WAVELEN']    # microns

    # fill in the hole in the continuum image since this otherwise messes up the slit image and spectroastrometry
    # (effectively a star and/or convolution by the PSF would do this in practice)
    # just do this by making the central region constant; adjust the area manually
    icen = int(hd['CRPIX1'])
    jcen = int(hd['CRPIX2'])
    ifill = 2
    jfill = 2
    median_continuum = np.nanmedian(continuum[jcen-jfill:jcen+jfill, icen-ifill:icen+ifill])
    for i in range(icen-ifill, icen+ifill+1):
        for j in range(jcen-jfill, jcen+jfill+1):
            continuum[j, i] = median_continuum

    cube = SpectralCube(data=line+continuum, wcs=wcs.WCS(hd)) 
    flux = cube.moment0().value

    # define path and direction starting from north
    nx = cube.header['NAXIS1']
    ny = cube.header['NAXIS2']
    i0 = cube.header['CRPIX1']
    j0 = cube.header['CRPIX2']

    reverse = False
    if slitPA >= 180:
        slitPA -= 180
        reverse = True

    if np.abs(slitPA-90) < 1:
        path_start = (0, j0)
        path_end = (nx-1, j0)
    elif (slitPA <= 45) or (slitPA >= 135):
        t = np.tan(np.radians(slitPA))
        path_start = (i0-(ny-1-j0)*t, ny-1)
        path_end = (i0+j0*t, 0)
    else:
        t = np.tan(np.radians(slitPA))
        path_start = (0, j0+i0/t)
        path_end = (nx-1, j0-(nx-1-i0)/t)
    if reverse:
        path_start, path_end = path_end, path_start
    path = Path([path_start, path_end], width=slit_width)

    pv = extract_pv_slice(cube=cube, path=path, spacing=1)
    SA_image = pv.data.T
    ns = pv.header['NAXIS1']
    nv = pv.header['NAXIS2']
    is0 = pv.header['CRPIX1']
    iv0 = pv.header['CRPIX2']
    ds = pv.header['CDELT1'] * 3.6e6    # milli-arcsec
    dv = pv.header['CDELT2'] / 1000     # km/s
    s0 = pv.header['CRVAL1'] * 3.6e6
    v0 = pv.header['CRVAL2'] / 1000
    s = (np.arange(ns) - ns/2 + 0.5) * ds   # make this zero at the center
    v = v0 + (np.arange(nv) - iv0) * dv

    # calculate the mean slit position (and std) at each velocity = spectroastrometric signal
    SA_mean = np.zeros(nv)
    SA_stdev = np.zeros(nv)
    for i in range(nv):
        # is the velocity vector increating or decreasing with index?
        if v[1]-v[0] > 0:
            Iv = SA_image[:, i]
        else:
            Iv = SA_image[:, nv-1-i]
        Ivsum = np.nansum(Iv)
        mu = np.nansum(s * Iv) / Ivsum
        SA_mean[i] = mu
        SA_stdev[i] = np.sqrt((np.nansum(s**2 * Iv) / Ivsum) - mu**2)

    if plot:
        # turn off interactive plots (i.e., just make png files)
        plt.ioff()
        fig = plt.figure(figsize=(18, 6))

        # show moment 0 image and slit path
        ax0 = plt.subplot(131, projection=cube.wcs.celestial)
        norm0 = ImageNormalize(vmin=0, vmax=flux.max(), stretch=LogStretch())
        ax0.imshow(flux, norm=norm0, origin='lower')
        #path.show_on_axis(ax0, spacing=3, color='r', alpha=0.1)
        ax0.arrow(path_start[0], path_start[1], \
                  path_end[0]-path_start[0], path_end[1]-path_start[1], \
                  color='w', length_includes_head=True, width=1, alpha=0.5)
        ax0.coords[0].set_format_unit(u.arcsec)
        ax0.coords[1].set_format_unit(u.arcsec)
        ax0.set_xlabel(r"$\Delta\alpha$ [arcsec]")
        ax0.set_ylabel(r"$\Delta\delta$ [arcsec]")

        # plot the spectroastrometric image
        ax1 = plt.subplot(132)
        im_bounds = [v.min(), v.max(), s.min(), s.max()]
        v_lims = [v.min(), v.max()]
        slit_lims = [s.min(), s.max()]
        norm1 = ImageNormalize(vmin=0, vmax=SA_image.max(), stretch=LogStretch())
        p1 = ax1.imshow(SA_image, origin='lower', cmap='viridis', norm=norm1, extent=im_bounds, aspect='auto')
        ax1.set_ylim(slit_lims)
        ax1.set_xlim(v_lims)
        ax1.set_ylabel('Slit offset [milli-arcsec]')
        ax1.set_xlabel('Velocity [km/s]')
        #cb = plt.colorbar(p1, ax=ax1, pad=0.05)
        #cb.set_label(r'Intensity [Jy / pix]', rotation=90, labelpad=17)
        
        # plot the flux-weighted spectroscopic signal
        ax2 = plt.subplot(133)
        ax2.set_aspect('auto')
        ax2.plot(v, SA_mean, lw=3)
        ax2.fill_between(v, SA_mean-SA_stdev, SA_mean+SA_stdev, alpha=0.3)
        ax2.set_xlim(v_lims)
        ax2.set_xlabel('Velocity [km/s]')
        ax2.set_ylabel('Slit offset [milli-arcsec]')
        
        fig.tight_layout()
        fig.savefig('./figures/spectroastrometry.png')

    return v, SA_mean, SA_stdev