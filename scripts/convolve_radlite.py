# #### read in radlite image and convolve it!

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from tqdm import tqdm

def observe(fwhm=None):
    if fwhm == None:
        print("Usage: observe(fwhm) where fwhm is in arcseconds")
        return

    path = "./radlite/"
    filename = "radlite_image.fits"

    hdu = fits.open(path+filename)
    continuum = hdu[1].data
    hd1 = hdu[1].header
    line = hdu[2].data
    hd2 = hdu[2].header

    sigma = fwhm / np.sqrt(8*np.log(2))
    sigma_pix = sigma / (hd1['CDELT2']*3600)
    print(f"Convolution sigma = {sigma_pix:.1f} pixels")

    gauss = Gaussian2DKernel(sigma_pix)
    continuum_convolved = convolve_fft(continuum, gauss)
    nk = line.shape[0]
    line_convolved = np.zeros(line.shape)
    for k in tqdm(range(nk)):
        line_convolved[k,:,:] = convolve_fft(line[k,:,:], gauss)

    # set up the headers and write them out as 2 extensions in a single fits file
    hdr = fits.Header()
    hdr['COMMENT'] = 'convolved radlite image'
    hdr['COMMENT'] = str(f'FWHM = {fwhm} arcseconds')
    hdr['COMMENT'] = 'continuum image is first extension'
    hdr['COMMENT'] = 'line datacube is second extension'
    hdr['COMMENT'] = 'created by convolve_radlite in radlite_for_radmc3d'

    hdu0 = fits.PrimaryHDU(header=hdr)
    hdu1 = fits.ImageHDU(continuum_convolved)
    hdu1.header = hd1
    hdu2 = fits.ImageHDU(line_convolved)
    hdu2.header = hd2

    fileout = 'radlite_convolved.fits'
    hdu_list = fits.HDUList([hdu0, hdu1, hdu2])
    hdu_list.writeto(path+fileout, overwrite=True)
    print(f"Image written to {fileout}")
