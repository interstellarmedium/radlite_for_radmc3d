These are some notebooks and required files to use radmc3dPy to create a disk model and process through radlite. The procedure to use is as follows:

1. Edit the model parameters in make_radmc3d_model.ipynb and run it. This will create a series of files in both the new and old format for radmc3d (new format) and radlite (old format). The radlite files and template notebook are moved to radmc_ouputs.

2. (Optional). Plot the radmc density and temperature files, and make images using plot_radmc3d_model.

3. Go into the newly created directory radmc_outputs and run the notebook run_radlite.ipynb. With the default parameters, this will create a CO spectrum. Edit model.json and spectrum.json as required.

Note that there are files which must be included here for everything to work. Do not edit them! But add other dust opacity or molecular line data models as needed. The files are dustkappa_silicate.inp, line.inp, molecule_co.inp.
