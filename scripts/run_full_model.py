# ### runs the full suite of modeling scripts to go from radmc to radlife to producing an image, spectrum and SA signal

import os
import numpy as np
import make_model
import run_radlite
import spectrum
import spectroastrometry


cwd = os.getcwd()

# model parameters
mstar = ['1.0*ms']
tstar = [4000]
rstar = ['2.0*rs']
mdust = '1e-5*ms'
dusttogas = 0.001,
rin = '0.05*au'
rdisk = '5*au'
gap_rin = '[0.0*au]'
gap_rout = '[0.0*au]',
Tmid = 300
Tatm = 700
Tmax = 2500

wind = True
D_collim = 1
v_terminal = 100

print('-*-'*25)
print('*** MAKE MODEL')
make_model.setup(mstar=mstar, tstar=tstar, rstar=rstar, mdust=mdust, dusttogas=dusttogas,
          rin=rin, rdisk=rdisk, gap_rin=gap_rin, gap_rout=gap_rout,
          Tmid=Tmid, Tatm=Tatm, Tmax=Tmax)

print('-*-'*25)
print('*** RUN RADLITE')
os.chdir(cwd+'/radlite')
if wind:
    run_radlite.write_wind_parameters(D_collim=D_collim, v_terminal=v_terminal)
run_radlite.make_plots(wind=wind)
run_radlite.make_image(wind=wind)
run_radlite.make_spectrum(wind=wind)
os.chdir(cwd)

print('-*-'*25)
print('*** PLOT SPECTRUM')
spectrum.plot()

print('-*-'*25)
print('*** SPECTROASTROMETRY')
v, SA, SA_err = spectroastrometry.measure(slitPA=45)


print('*** FULL MODEL SCRIPT FINISHED ***')
print('-*-'*25)
