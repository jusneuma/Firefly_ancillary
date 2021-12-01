# Plotting scripts for the MaNGA FIREFLY VAC - DR17

We have written two small python scripts as an example of how to plot 2D maps and 1D profiles from the MaNGA FIREFLY VAC. The fits file of the VAC should be downloaded and by default placed in the same directory as this script. The script should then be run directly by typing 
```
python plot_maps_vac.py [plate] [ifu]
python plot_profiles_vac.py [plate] [ifu]
```

on the commandline. Replace [plate] and [ifu] with the numbers of the datacube that you want to plot.

By default, the scripts use the 'miles' version of the VAC, this can be changed in the scripts.
