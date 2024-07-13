This is modification of Klaus Pontoppidan's pyradlite (https://github.com/pontoppi/radlite) with small edits to radlite.py and radmc.py to make it work with radmc3dPy, and implementation of making the imaging part work. Example notebooks and scripts are included to demonstrate utility (the former are more for use as a tutorial and the scripts are better for fitting to data).

Note that to allow for finer gridding, which is important for the inner disk, it is worth recompiling the fortran radlite code that these python scripts call after changing the maximum grid size defaults in configure.h to:

#define FRSIZE_X      420

#define FRSIZE_Y      420

#define FRSIZE_FREQ      401

#define FRSIZE_MU   101

(be careful about making them too big as it may compile but then crashes when executing).

parameterized_disk.py is similar to ppdisk.py in the model examples of radmc3dPy but calculates the dust temperature analytically rather than relying on radmc3d mctherm (which can be slow or noisy).
Move it over to where the other radmc3dPy models are. For me, this is $HOME/.local/lib/python3.9/site-packages/radmc3dPy/models

Instructions on running this version of radlite are in a notion page posted here
https://excellent-decade-3d5.notion.site/radlite-for-radmc3d-documentation-3cce4e0e195244a0b38f4d15226a99fa
