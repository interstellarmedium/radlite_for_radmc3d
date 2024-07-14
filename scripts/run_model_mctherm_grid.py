# ### runs the full suite of modeling scripts to go from radmc to radlife to producing an image, spectrum and SA signal

import os
import single_model_mctherm

cwd = os.getcwd()

###------------------------------------------------
# model parameters
# customize values/ranges and how to name the directory

# fixed values
mstar = ['1.0*ms']
rstar = ['2.0*rs']
tstar = [4000]
dusttogas = 0.01,
rdisk = '4.9*au'
nphot = 1e6

# parameters with a range of values - place in a list
mdisk = ['1e-4*ms', '1e-3*ms']
rin = ['0.05*au', '0.1*au']

wind = True
D_collim = 1
v_terminal = 100

# recommended to save 75% disk space per model
clean_space = True
###------------------------------------------------

if wind:
    prefix = 'wind'
    wind_dict = {'flag':True, 'D_collim':D_collim, 'v_terminal':v_terminal}
else:
    prefix = 'keplerian'
    wind_dict = {'flag':False}

for mdisk1 in mdisk:
    mdisk_str = 'Mdisk'+mdisk1.split('*')[0]
    for rin1 in rin:
        rin_str = 'Rin'+rin1.split('*')[0]
        keys = {'model':{'mstar':mstar, 'rstar':rstar, 'tstar':tstar,
                         'mdisk':mdisk1, 'dusttogas':dusttogas,
                         'rin':rin1, 'rdisk':rdisk},
                'wind':wind_dict}
        working_dir = cwd+'/'+prefix+str(f'_{mdisk_str}_{rin_str}/')
        single_model_mctherm.run_model(cwd, working_dir, clean_space=clean_space, **keys)
