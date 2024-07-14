# ### create a model to run through radlite
# ### model creates the (dust) temperature profile so no need to run radmc3d mctherm

import os
import numpy as np
import radmc3dPy

def setup(mstar=['1.0*ms'], tstar=[4000], rstar=['2.0*rs'], mdisk='1e-3*ms', dusttogas=0.01,
          rin='0.05*au', rdisk='5*au', gap_rin='[0.0*au]', gap_rout='[0.0*au]',
          Tmid=300, Tatm=700, Tmax=2500):

    keys = {'model':'parameterized_disk', 'mstar':mstar, 'rstar':rstar, 'tstar':tstar,
            'mdisk':mdisk, 'dusttogas':dusttogas,
            'rin':rin, 'rdisk':rdisk,
            'xbound':'[0.05*au, 0.5*au, 5.0*au]', 'nx':'[150, 150]',
            'nz':'0', 'binary':False,
            'Tmid':Tmid, 'Tatm':Tatm, 'Tmax':Tmax, 'writeDustTemp':True}

    # setup problem_params.log with default parameters
    radmc3dPy.analyze.writeDefaultParfile(keys['model'])

    # update the parameters in problem_params.inp and write out files in new radmc3d format
    radmc3dPy.setup.problemSetupDust(**keys)

    # translate temperature array to old radmc format (settings keys to old does not do this unfortunately)
    data = radmc3dPy.analyze.readData(dtemp=True)
    Tdust = data.dusttemp[:,:,0,0]
    nr, nt = Tdust.shape
    fname = "dusttemp_final.dat"
    with open(fname, "w") as wfile:
        print("Writing " + fname)
        wfile.write(f"   1    {nr:d}   {nt//2:d}   1\n")
        wfile.write(" \n")
        wfile.write("   1\n")
        for i in range(nr):
            for j in range(nt//2):
                wfile.write(f"{Tdust[i,j]:.7f}\n")

    fname = "dusttemp.info"
    with open(fname, "w") as wfile:
        print("Writing " + fname)
        wfile.write("  -2\n")
        wfile.write("   1\n")
        wfile.write("   1\n")
        wfile.write("   1\n")
        wfile.write("   1\n")

    # make a directory for the radmc new format outputs
    outputdir = "./radmc"
    print(f"Moving radmc files to output directory {outputdir}")
    if os.path.exists(outputdir):
        print("Will overwrite existing files")
    else:
        print("Directory does not exist; will create")
        os.makedirs(outputdir)

    filelist = ["problem_params.inp", "dustopac.inp", "wavelength_micron.inp", "amr_grid.inp", "stars.inp", \
                "dust_density.inp", "dust_temperature.dat", "radmc3d.inp", \
                "plot_radmc3d_model.ipynb"]
    for file in filelist:
        os.system("mv "+file+" radmc/")

    # copy rather than move this file as we still need it for the old format
    os.system("cp dustkappa_silicate.inp radmc/")

    # create problem_params.inp, update parameters, and write out files for radlite)
    # (note that radlite uses the old format of radmc3d)
    radmc3dPy.analyze.writeDefaultParfile(keys['model'])
    radmc3dPy.setup.problemSetupDust(**keys, old=True)

    # get rid of the new format dust temperature file which is not needed
    print('Removing dust_temperature.dat')
    os.remove('dust_temperature.dat')

    # make a directory for the radmc outputs for radlite to use
    outputdir = "./radlite"
    print(f"Moving radmc files to output directory {outputdir}")
    if os.path.exists(outputdir):
        print("Will overwrite existing files")
    else:
        print("Directory does not exist; will create")
        os.makedirs(outputdir)

    filelist = ["problem_params.inp", "dustdens.inp", "dustopac.inp", "dusttemp_final.dat", "dusttemp.info", \
                "frequency.inp", "dustopac_1.inp", "radius.inp", "theta.inp", "starinfo.inp", "starspectrum.inp", "radmc.inp"]
    for file in filelist:
        os.system("mv "+file+" radlite/")

    # move files over to radmc_outputs for radlite to run
    # (I tried and failed to hide any errors that come up if you run this a second time and the files have already been moved...)
    filelist = ["run_radlite.py", "data_hitran.json", "model*.json", "spectrum.json", "line.inp", "molecule_co.inp"]
    for file in filelist:
        cmd = "mv "+file+" radlite/"
        x = os.system(cmd + "> /dev/null")
        if x==0:
            print(cmd)
    return
