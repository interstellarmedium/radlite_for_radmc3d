# ### write out radlite image as a fitsfile
# #### https://docs.astropy.org/en/stable/io/fits/
# #### Note that, like C (and unlike Fortran), Python is 0-indexed and the indices have the slowest axis first and fastest changing axis last; that is, for a 2D image, the fast axis (X-axis) which corresponds to the FITS NAXIS1 keyword, is the second index.

import numpy as np
import astropy.constants as c
from astropy.io import fits

def read_and_write(path="./rundir_image/outputdir/", filename="lineposvel_moldata_0_1.dat", dist=100):
    # load the image
    # see RADLITE/telescope.F for the format aubroutine calc_write_line_posvel
    # note units are brightness temperature = lambda**2 * F_nu / 2*k

    with open(path+filename) as f:
        try:
            blank = f.readline()
            iformat = f.readline()
            molname1 = f.readline()
            molname2 = f.readline()
            dist0, v0, incl = [float(abc) for abc in f.readline().split()]
            trans = f.readline()
            nu0 = float(f.readline())
            nv = int(f.readline())
            nx, ny, dx, dy, rotang, xoff, yoff = [float(abc) for abc in f.readline().split()]
        except:
            print("Error reading header")
            exit()
        
        nx = int(nx)
        ny = int(ny)
        x = xoff + (np.arange(nx) - nx/2 + 0.5) * dx
        y = yoff + (np.arange(ny) - ny/2 + 0.5) * dy
        
        blank = f.readline()
        
        v = np.empty(nv)
        for k in range(nv):
            v[k] = float(f.readline())
        dv = v[1] - v[0]
    
        Tb = np.empty((nv, ny, nx))
        tau = np.empty((nv, ny, nx))
        for iv in range(nv):
            blank = f.readline()
            for iy in range(ny):
                for ix in range(nx):
                    Tb1, tau1 = [float(ab) for ab in f.readline().split()]
                    Tb[iv, ix, iy] = Tb1
                    tau[iv, ix, iy] = tau1

    # convert to flux in Jy at 1pc
    Inu = 1e23 * (2*c.k_B.cgs.value/c.c.cgs.value**2) * nu0**2 * Tb
    Fnu = Inu * np.abs(dx*dy)/c.pc.cgs.value**2

    # convert step size to au as in model.json
    dx /= c.au.cgs.value
    dy /= c.au.cgs.value

    # convert flux to Jy and step size to degrees at source distance
    # note that 1au at 1pc is 1 arcsecond
    Fnu /= dist**2
    dx /= (3600*dist)
    dy /= (3600*dist)

    # separate line from continuum
    #(linearly interpolate across endpoints as for spectra in radlite)
    # I prefer to order my arrays differently from the fortran code
    nv, nx, ny = Fnu.shape
    line = np.zeros((nx, ny, nv))
    continuum = np.zeros((nx, ny))

    iv = np.arange(nv)
    for iy in range(ny):
        for ix in range(nx):
            f0 = Fnu[0, ix, iy]
            f1 = Fnu[-1, ix, iy]
            #tempys = [Fnu[0, ix, iy], Fnu[-1, ix, iy]]
            #interpolate_spec = interp1d([0,nv-1], tempys)
            #continuum[ix, iy] = interpolate_spec(nv/2)
            continuum[ix, iy] = (f0+f1)/2
            for iv in range(nv):
                #line[ix, iy, iv] = Fnu[iv, ix, iy] - interpolate_spec(iv)
                line[ix, iy, iv] = Fnu[iv, ix, iy] - (f0 + (f1-f0)*iv/nv)

    # set up the headers and write them out as 2 extensions in a single fits file
    hdr = fits.Header()
    hdr['COMMENT'] = 'radlite image'
    hdr['COMMENT'] = 'continuum image is first extension'
    hdr['COMMENT'] = 'line datacube is second extension'
    hdr['COMMENT'] = 'created by image2fits notebook in radlite_for_radmc3d'

    hdu0 = fits.PrimaryHDU(header=hdr)

    hdu1 = fits.ImageHDU(continuum.T)
    hd1 = hdu1.header
    hd1.set('OBJECT', filename, 'continuum')
    hd1.set('CRPIX1', nx/2-0.5, 'image center')
    hd1.set('CRPIX2', ny/2-0.5, 'image center')
    hd1.set('CDELT1', dx, 'degrees')
    hd1.set('CDELT2', dy, 'degrees')
    hd1.set('CRVAL1', 0.0, 'degrees')
    hd1.set('CRVAL2', 0.0, 'degrees')
    hd1.set('CTYPE1', 'RA---SIN')
    hd1.set('CTYPE2', 'DEC--SIN')
    hd1.set('CUNIT1', 'deg')
    hd1.set('CUNIT2', 'deg')
    hd1.set('BUNIT', 'Jy/pix')
    hd1.set('RESTFRQ', nu0, 'line frequency in Hz')
    hd1.set('WAVELEN', 1e6*c.c.si.value/nu0, 'wavelength in microns')

    hdu2 = fits.ImageHDU(line.T)
    hd2 = hdu2.header
    hd2.set('OBJECT', filename, 'line')
    hd2.set('CRPIX1', nx/2-0.5, 'image center')
    hd2.set('CRPIX2', ny/2-0.5, 'image center')
    hd2.set('CRPIX3', nv/2-0.5, 'center of velocity axis')
    hd2.set('CDELT1', dx, 'degrees')
    hd2.set('CDELT2', dy, 'degrees')
    hd2.set('CDELT3', dv, 'step size in km/s')
    hd2.set('CRVAL1', 0.0, 'degrees')
    hd2.set('CRVAL2', 0.0, 'degrees')
    hd2.set('CRVAL3', v0, 'central velocity in km/s')
    hd2.set('CTYPE1', 'RA---SIN')
    hd2.set('CTYPE2', 'DEC--SIN')
    hd2.set('CTYPE3', 'VELO-LSRK')
    hd2.set('CUNIT1', 'deg')
    hd2.set('CUNIT2', 'deg')
    hd2.set('CUNIT3', 'km/s')
    hd2.set('BUNIT', 'Jy/pix')
    hd2.set('RESTFRQ', nu0, 'line frequency in Hz')
    hd2.set('WAVELEN', 1e6*c.c.si.value/nu0, 'wavelength in microns')

    hdu_list = fits.HDUList([hdu0, hdu1, hdu2])
    hdu_list.writeto('radlite_image.fits', overwrite=True)
