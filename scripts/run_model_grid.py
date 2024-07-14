# ### runs the full suite of modeling scripts to go from radmc to radlife to producing an image, spectrum and SA signal

import os
import single_model

cwd = os.getcwd()

###------------------------------------------------
# model parameters
# customize values/ranges and how to name the directory

# fixed values
mstar = ['1.0*ms']
rstar = ['2.0*rs']
tstar = [4000]
dusttogas = 0.01
rdisk = '4.0*au'
Tmax = 2500

# parameters with a range of values - place in a list
mdisk = ['1e-4*ms', '1e-3*ms']
rin = ['0.05*au', '0.1*au']
Tmid = [300, 500]
Tatm = [700, 1000]

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
        for Tmid1 in Tmid:
            for Tatm1 in Tatm:
                keys = {'model':{'mstar':mstar, 'rstar':rstar, 'tstar':tstar,
                                 'mdisk':mdisk1, 'dusttogas':dusttogas,
                                 'rin':rin1, 'rdisk':rdisk,
                                 'Tmid':Tmid1, 'Tatm':Tatm1, 'Tmax':Tmax},
                        'wind':wind_dict}
                working_dir = cwd+'/'+prefix+str(f'_{mdisk_str}_{rin_str}_Tmid{Tmid1}_Tatm{Tatm1}/')
                single_model.run_model(cwd, working_dir, clean_space=clean_space, **keys)
