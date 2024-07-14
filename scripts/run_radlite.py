# ## Run radlite to create an image and/or spectrum based on the output of make_model.ipynb
# ## Includes a super-class to customize velocity field
# ## and plotting routines to visualize the density, temperature, and velocity structure
# ## 3/5/24 jpw

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.interpolate import RectBivariateSpline
import astropy.constants as c
from astropy.io import fits

#Import Python-RADLite package
sys.path.append("/Users/jpw/py/radlite_for_radmc3d")
import radlite as RDL

# create a super class to allow for more flexible velocity structure
# wind velocity prescription from Kurosawa 2006; collimation distance in au, terminal velocity in km/s
class Wind(RDL.RadliteModel):
    def __init__(self, infilename, hitranfilename, radmcfilepath):
        super().__init__(infilename, hitranfilename, radmcfilepath)
        
    def _calc_velocity(self):
        # this function supersedes the one in radlite which only does Keplerian motion

        # pass wind parameters into radlite through a json file - this somewhat clunky way was the best I could do
        with open('wind_parameters.json', 'r') as openfile:
            winddict = json.load(openfile)
            for key in winddict:
                self._attrdict[key] = winddict[key]["value"]

        # radlite spherical coordinates
        r = self.get_attr("radius") #Radii
        rlen = len(r)
        theta = self.get_attr("theta") #Co-lattiude in Dullemond parlance (radians)
        tlen = len(theta)//2 #Only take the top half since we assume symmetry
        theta = theta[0:tlen]

        mstar = self.get_attr("starinfo")["mstar"] #Stellar mass
        rstar = self.get_attr("starinfo")["rstar"] #Stellar radius
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating velocity field components...")
            print("Used starinfo.inp file for mstar and rstar...")

        G0 = 6.67259e-8
        au = 1.49597870e13

        if self.get_attr("verbose"):
            print("Kurosawa wind model")
            print('Collimation distance (au) =', self.get_attr("D_collim"))
            print('Terminal speed (km/s) =', self.get_attr("v_terminal"))
            print('Wind scale length (au) =', self.get_attr("Rs"))
            print('Power law index =', self.get_attr("beta"))
        D_collim = self.get_attr("D_collim") * au
        v_terminal = self.get_attr("v_terminal") * 1e5
        Rs = self.get_attr("Rs") * au
        beta = self.get_attr("beta")

        # sound speed in cm/s at 1000 K
        cs = 6e4

        # used for pressure scale height (from parameterized_disk.py)
        hrdisk = 0.1          # ratio of the pressure scale height over radius at hrpivot
        hrpivot = 100 * au    # reference radius at which Hp/R is taken
        plh = 1/7             # radial power law index

        # use Kurosawa notation where (w,z) are cylindrical coordinates
        rr, tt = np.meshgrid(r, theta)
        ww = rr * np.sin(tt)
        zz = rr * np.cos(tt)
        wcross = ww*D_collim/(zz+D_collim)
        S = np.sqrt(wcross**2 + D_collim**2)
        l = np.sqrt((ww-wcross)**2 + zz**2)
        vp = cs + (v_terminal-cs) * (1 - Rs/(l+Rs))**beta
        ksi = np.arcsin(D_collim * np.sin(np.pi-tt) / (S+l))

        # zero out the wind below some number of pressure scale heights (where the temperature changes)
        Hp = hrdisk * (ww/hrpivot)**plh
        no_wind = np.where(zz < 2*Hp*ww)[0]
        vp[no_wind] = 0

        vel_radial = vp * np.cos(ksi)
        vel_theta = vp * np.sin(ksi)
        vel_phi = np.sqrt(G0*mstar/rr) * wcross/ww

        ##Below Section: RECORD velocity components and EXIT
        self._set_attr(attrname="velocity_radial", attrval=vel_radial,
                            attrunit=r"cm s$^{-1}$")
        self._set_attr(attrname="velocity_theta", attrval=vel_theta,
                            attrunit=r"cm s$^{-1}$")
        self._set_attr(attrname="velocity_phi", attrval=vel_phi,
                            attrunit=r"cm s$^{-1}$")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating velocity components!\n")
        return

def image2fits(path="./rundir_image/outputdir/", filename="lineposvel_moldata_0_1.dat"):
    # write the radlite image as a fits file
    # see RADLITE/telescope.F for the format aubroutine calc_write_line_posvel
    # note units are brightness temperature = lambda**2 * F_nu / 2*k
    print("Writing radlite_image.fits")

    # read in header of radlite file
    # note that radmc calculations are made for a distance of 1pc
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
            print(f"Error reading header in {path}{filename}")
            exit()

        # read in model setup files to get physical distance and inclination
        with open("spectrum.json", "r") as openfile:
            dict = json.load(openfile)
        dist = dict["dist"]["value"]
        print(f"Reading distance from spectrum.json, d={dist} pc")

        with open("model_image.json", "r") as openfile2:
            dict = json.load(openfile2)
        incl = dict["incl"]["value"]
        print(f"Reading inclination from model_image.json, i={incl} degrees")
        
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
    # include inclination and distance as in radlite_spectrum.fits
    hdr = fits.Header()
    hdr['COMMENT'] = 'radlite image'
    hdr['COMMENT'] = 'continuum image is first extension'
    hdr['COMMENT'] = 'line datacube is second extension'
    hdr['COMMENT'] = 'created by image2fits notebook in radlite_for_radmc3d'
    hdr["INCL_deg"] = incl
    hdr["DIST_pc"] = dist

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
    return


def write_wind_parameters(D_collim=1, v_terminal=100, Rs=1, beta=2):
    winddict = {
        "D_collim": {"value":D_collim, "comment":"Collimation distance in au"},
        "v_terminal": {"value":v_terminal, "comment":"Terminal velocity in km/s"},
        "Rs": {"value":Rs, "comment":"Scale length in au"},
        "beta": {"value":beta, "comment":"power law index"}
    }
    with open("wind_parameters.json", "w") as fileout:
        json.dump(winddict, fileout, indent=4)

    return


def make_image(wind=False):
    infilename = "model_image.json"

    radmcfilepath = "./"
    hitranfilename = "data_hitran.json"
    if wind:
        if os.path.exists('wind_parameters.json'):
            myMod = Wind(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)
        else:
            print('Write wind parameters first')
            return
    else:
        myMod = RDL.RadliteModel(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)
    myMod.run_radlite()

    image2fits()
    return


def make_spectrum(wind=False):
    infilename = "model_spectrum.json"
    inspecfilename = "spectrum.json"

    radmcfilepath = "./"
    hitranfilename = "data_hitran.json"
    if wind:
        if os.path.exists('wind_parameters.json'):
            myMod = Wind(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)
        else:
            print('Write wind parameters first')
            return
    else:
        myMod = RDL.RadliteModel(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)
    myMod.run_radlite()

    mySpec = RDL.RadliteSpectrum(infilename=inspecfilename)
    mySpec.gen_spec()
    mySpec.write_fits('radlite_spectrum.fits', overwrite=True)
    return


def make_plots(wind=False):
    # turn off interactive plots (i.e., just make png files)
    plt.ioff()

    infilename = "model_image.json"
    radmcfilepath = "./"
    hitranfilename = "data_hitran.json"
    if wind:
        if os.path.exists('wind_parameters.json'):
            myMod = Wind(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)
        else:
            print('Write wind parameters first')
            return
    else:
        myMod = RDL.RadliteModel(infilename=infilename, hitranfilename=hitranfilename, radmcfilepath=radmcfilepath)

    # regrid from radlite spherical coordinates (r,theta) into (R,Z) cylindrical coordinates
    r = myMod.get_attr("radius") / 1.49598e13
    theta = myMod.get_attr("theta")
    nr = r.size
    nt = theta.size

    logdens = np.log10(myMod.get_attr("gasdensity"))
    Tgas = myMod.get_attr("gastemperature")
    dens_spline = RectBivariateSpline(theta, r, logdens)
    Tgas_spline = RectBivariateSpline(theta, r, Tgas)

    nR = 200
    R = np.logspace(-1.4, 0.7, nR)
    nZ = 100
    Z = np.logspace(-2, 0.7, nZ)

    logdens_RZ = np.zeros((nZ, nR))
    Tgas_RZ = np.zeros((nZ, nR))
    for j in range(nR):
        for i in range(nZ):
            r1 = np.sqrt(R[j]**2 + Z[i]**2)
            t1 = np.arctan(R[j]/Z[i])
            logdens_RZ[i,j] = dens_spline(t1, r1)
            Tgas_RZ[i,j] = Tgas_spline(t1, r1)
    
    Tmax = np.max(Tgas_RZ)
    nlevs = int(Tmax/100) - 4
    Tlevs = (3 + np.arange(nlevs)) * 100

    # spherical coordinate plot (not used)
    #fig = plt.figure(figsize=(10, 5))
    #ax = fig.add_subplot(111)
    #c1 = ax.contourf(r, np.degrees(theta), logdens, 30, cmap='jet')
    #ax.set_xlabel('r [au]')
    #ax.set_ylabel(r'theta [radians]')
    #ax.set_xscale('log')
    #ax.set_title('Gas density and temperature (spherical coordinates)')
    #cb = plt.colorbar(c1)
    #cb.set_label(r'$\log_{10}{\rho}$', rotation=270.)
    #c2 = ax.contour(r, np.degrees(theta), Tgas, Tlevs,  colors='w', linestyles='solid')
    #ax.clabel(c2, inline=1, fontsize=10)
    #plt.savefig('gas_density_temperature_spherical.png')
    #plt.close(fig)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    c1 = ax.contourf(R, Z, logdens_RZ, 30, cmap='jet')
    ax.set_xlabel('R [au]')
    ax.set_ylabel(r'Z [au]')
    ax.set_title('Gas density and temperature')
    ax.set_xscale('log')
    ax.set_yscale('log')
    cb = plt.colorbar(c1)
    cb.set_label(r'$\log_{10}{\rho}$', rotation=270.)
    c2 = ax.contour(R, Z, Tgas_RZ, Tlevs,  colors='w', linestyles='solid')
    ax.clabel(c2, inline=1, fontsize=10)
    #----------------------------------------------------------
    # calculate tau=1 surface at 0.5 microns
    # this is from midplane vertically up so in cyleindrical coordinates
    # BUT this does not agree with the getTau module using the new format so do not use
    # (also need to uncomment a line in radmc.py to make it work so a bit untested perhaps?)
    #import radmc
    #RM = radmc.radmc_model('./')
    #R_thick, Z_thick = RM.surface_at_tau(1, wave=0.5)
    #R_thick /= 1.5e13
    #Z_thick /= 1.5e13
    #ax.plot(R_thick, Z_thick,  color='black', ls='solid', lw=5, alpha=0.3)
    #----------------------------------------------------------
    plt.savefig('../figures/gas_density_temperature.png')
    plt.close(fig)

    # velocity field
    nt2 = theta.size // 2
    theta2 = theta[0:nt2]
    vrad = myMod.get_attr("velocity_radial") / 1e5
    vtheta = myMod.get_attr("velocity_theta") / 1e5
    vphi = myMod.get_attr("velocity_phi") / 1e5

    if np.max(vrad)-np.min(vrad) > 0 or np.max(vtheta)-np.min(vtheta) > 0:
        rr, tt2 = np.meshgrid(r, theta2)
        vR = vrad*np.sin(tt2) - vtheta*np.cos(tt2)
        vZ = vrad*np.cos(tt2) + vtheta*np.sin(tt2)
        vR_spline = RectBivariateSpline(theta2, r, vR)
        vZ_spline = RectBivariateSpline(theta2, r, vZ)
        vphi_spline = RectBivariateSpline(theta2, r, vphi)

        # sparse grid because otherwise the quiver plot is too crowded
        nR = 15
        Rv = np.logspace(-1.4, 0.7, nR)
        nZ = 10
        Zv = np.logspace(-2, 0.7, nZ)
        vrad_RZ = np.zeros((nZ, nR))
        vtheta_RZ = np.zeros((nZ, nR))
        vR_RZ = np.zeros((nZ, nR))
        vZ_RZ = np.zeros((nZ, nR))
        vphi_RZ = np.zeros((nZ, nR))
        for j in range(nR):
            for i in range(nZ):
                r1 = np.sqrt(Rv[j]**2 + Zv[i]**2)
                t1 = np.arctan(Rv[j]/Zv[i])
                vR_RZ[i,j] = vR_spline(t1, r1)
                vZ_RZ[i,j] = vZ_spline(t1, r1)
                vphi_RZ[i,j] = vphi_spline(t1, r1)

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        c1 = ax.contourf(R, Z, logdens_RZ, 30, cmap='jet')
        ax.set_xlabel('R [au]')
        ax.set_ylabel(r'Z [au]')
        ax.set_title('Gas density and non-azimuthal velocity field')
        ax.set_xscale('log')
        ax.set_yscale('log')
        cb = plt.colorbar(c1)
        cb.set_label(r'$\log_{10}{\rho}$', rotation=270.)
        ax.quiver(Rv, Zv, vR_RZ, vZ_RZ, color='white')
        plt.savefig('../figures/gas_non_azimuthal_velocity.png')
        plt.close(fig)

    return

def delete_rundir():
    # delete the rundir directories to save space
    # (only do once fits files have been created!)
    os.system('rm -r rundir_image')
    os.system('rm -r rundir_spectrum')
    return
