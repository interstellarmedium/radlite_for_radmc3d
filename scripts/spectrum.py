import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
sys.path.append('/Users/jpw/py/spectools_ir/')
from spectools_ir.utils import extract_hitran_data

def plot(path='./radlite/', filename='radlite_spectrum.fits', normalized=False):
    # set up plot parameters
    mpl.rc('xtick.major', size=5, pad=3, width=2)
    mpl.rc('ytick.major', size=5, pad=3, width=2)
    mpl.rc('xtick.minor', size=2, pad=3, width=2)
    mpl.rc('ytick.minor', size=2, pad=3, width=2)
    mpl.rc('axes', linewidth=2, labelsize=18, titlesize=14)
    mpl.rc('legend', fontsize=14)
    mpl.rc('lines', markersize=5)
    mpl.rc('font', size=12)

    # read in line lists
    hitran_CO = extract_hitran_data('CO', 4.5, 5.5, vup=1)
    hitran_CO2 = extract_hitran_data('CO', 4.5, 5.5, vup=2)
    hitran_CO3 = extract_hitran_data('CO', 4.5, 5.5, vup=3)
    hitran_13CO = extract_hitran_data('CO', 4.5, 5.5, isotopologue_number=2, vup=1)

    # load the spectrum
    # 0 is metadata (units, inclination, distance), no data
    # 1 is the spectrum (lambda, nu, emission = line, spectrum = line+continuum, continuum)
    # 2 is metadata on the transition (molecule, degeneracy, energy, etc)
    hdu = fits.open(path+filename)
    hdu.info()
    data = hdu[1].data
    wave = data["wavelength"]
    wmin = np.min(wave)
    wmax = np.max(wave)

    if normalized:
        flux = 1 + data["emission"]/data["continuum"]
    else:
        flux = data["spectrum"]
    fmin = np.nanmin(flux)
    fmax = np.nanmax(flux)

    # narrow down onto one transition
    #wmin = 4.6325
    #wmax = 4.6340

    # plot the full spectrum but not to screen, only to file
    plt.ioff()
    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(111)

    dy = fmax - fmin
    ymax = 1.05*fmax
    ymin = fmin - 0.2*dy
    y1 = ymin + 0.05*dy
    y2 = ymin + 0.15*dy
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(wmin, wmax)

    ax.plot(wave, flux, lw=1, color='darkblue')
    ax.set_xlabel(r'$\lambda$ ($\mu$m)', fontsize = 20, labelpad=10)
    if normalized:
        ax.set_ylabel('Line-to-continuum ratio', fontsize = 20, labelpad=10)
        figfile = 'spectrum_normalized.png'
    else:
        ax.set_ylabel('Flux (Jy)', fontsize = 20, labelpad=10)
        figfile = 'spectrum.png'

    for w in hitran_CO['wave']:
        if w > wmin and w < wmax:
            ax.plot([w, w], [y1, y2], color='blue')
    for w in hitran_CO2['wave']:
        if w > wmin and w < wmax:
            ax.plot([w, w], [y1, y2], color='red')
    for w in hitran_CO3['wave']:
        if w > wmin and w < wmax:
            ax.plot([w, w], [y1, y2], color='red', linestyle=':')
    for w in hitran_13CO['wave']:
        if w > wmin and w < wmax:
            ax.plot([w, w], [y1, y2], color='green')

    fig.tight_layout()
    fig.savefig('./figures/'+figfile)
    return

