These are some notebooks and required files to use radmc3dPy to create a disk model and process through radlite. The procedure to use is as follows:

1. Edit the model parameters in make_model.ipynb and run it. This will create a series of files in both the new and old format for radmc3d (new format) and radlite (old format). The radmc and radlite files are moved into newly created directories. Note that this calls on a model (which you can create or modify) in radmc3dPy that has a prescription for the dust temperature. If your model does not have that, you should run "radmc3d mctherm" and run make_radmc3d_model.ipynb instead (although I found that the temperature structure for the inner disk is slow/hard to calculate accurately using mctherm).

2. (Optional). Got into directory radmc and plot the density and temperature files, and make a scattered light image using plot_radmc3d_model.

3. Go into directory radlite and run the notebook run_radlite.ipynb. With the default parameters, this will create a CO spectrum. Edit model.json and spectrum.json as required.

Note that there are files which must be included here for everything to work. Do not edit them! But add other dust opacity or molecular line data models as needed. The files are dustkappa_silicate.inp, line.inp, molecule_co.inp.
