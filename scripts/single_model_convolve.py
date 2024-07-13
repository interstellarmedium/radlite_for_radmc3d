# ### runs the full suite of modeling scripts to go from radmc to radlife to producing an image, spectrum and SA signal

import os
import numpy as np
import make_model
import run_radlite
import spectrum
import spectroastrometry
import convolve_radlite

def run_model(cwd, working_dir, clean_space=False, **keys):
    print('-*-'*27)
    print('*** RUNNING MODEL IN ', working_dir)
    os.makedirs(working_dir)
    filelist = ['dustkappa_silicate.inp', 'line.inp', 'make_model.py', 'make_model_mctherm.py', 'model_image.json',
                'model_spectrum.json', 'molecule_co.inp', 'plot_radmc3d_model.ipynb', 
                'run_radlite.py', 'spectroastrometry.py', 'spectrum.json', 'spectrum.py', 'data_hitran.json']
    for file in filelist:
        os.system("cp "+file+" "+working_dir)
    os.chdir(working_dir)

    print('-*-'*27)
    print('*** MAKE MODEL')
    model_keys = keys['model']
    make_model.setup(**model_keys)
    #make_model.setup(mstar=mstar, tstar=tstar, rstar=rstar, mdust=mdust, dusttogas=dusttogas,
    #        rin=rin, rdisk=rdisk, gap_rin=gap_rin, gap_rout=gap_rout,
    #        Tmid=Tmid, Tatm=Tatm, Tmax=Tmax)
    os.makedirs(working_dir+'figures')
    if clean_space:
        os.system('rm -r radmc')

    print('-*-'*27)
    print('*** RUN RADLITE')
    os.chdir(working_dir+'radlite')
    wind = keys['wind']
    if wind['flag']:
        run_radlite.write_wind_parameters(D_collim=wind['D_collim'], v_terminal=wind['v_terminal'])
    run_radlite.make_plots(wind=wind['flag'])
    run_radlite.make_image(wind=wind['flag'])
    run_radlite.make_spectrum(wind=wind['flag'])  # this can take a few minutes and is not necessary if you only want the SA signal
    if clean_space:
        run_radlite.delete_rundir()
    os.chdir(working_dir)

    print('-*-'*27)
    print('*** PLOT SPECTRUM')
    spectrum.plot()

    print('-*-'*27)
    print('*** CONVOLUTION')
    convolve_radlite.observe(fwhm=0.05)

    print('-*-'*27)
    print('*** SPECTROASTROMETRY')
    v, SA, SA_err = spectroastrometry.measure(filename='radlite_convolved.fits', slitPA=0, outputfig='SA_minor.png')
    v, SA, SA_err = spectroastrometry.measure(filename='radlite_convolved.fits', slitPA=90, outputfig='SA_major.png')

    print('*** MODEL SCRIPT FINISHED ***')
    os.chdir(cwd)
    return
