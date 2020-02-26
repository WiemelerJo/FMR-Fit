# FMR-Fit
Small GUI for extracting parameters of ferromagnetic resonance measurements made in python 3.6.
Modules used:
matplotlib, numpy, lmfit(wrapper for scipy), PyQt5(Basicly PySide2), sympy and symengine(C wrapper of Symengine for Python)

# First tab: 'Fit'

At first open a measurement pressing ctrl + O, or anything that requires this file. 
Remember that only measurements made from Brukers Xepr converted to ascii will work! Alternatively any .dat file will work, if every column is arranged like Xepr does.

The next step will be to select a dataset (a measurement extracted from the ascii file will be called dataset in the following), by either typing the number or using your scrollwheel. You will see a change in the top right plot. There it is possible to select the viewport. Use the slider on the bottom first, to select the maximum viewing position, then proceed with the left slider to select the minimum. This will also set the fitting range. Faulty datasets can be detected by going trhu every measuremnt and wathcing the data plot. Once a faulty measuement is detected, it can be left out for fitting, by adding it to Dropped Points.

Then select a function for fitting. Till now the number of functions is restricted to only 1 line! And set the constraints as well as the inital values by again using your scrollwheel or typing into the box. 

Plot will first fit the measurement and plots it in a new window. The button below Plot, 'Set inital parameters to fit parameters', will take the fitted parameters and sets them as the new initial values.

For dynamic plotting it is advisable to set the initial values to fitted parameters using the dataset 0, to get the best result! By pushing the dynamic fit button, a window is opened asking for a place to save the fitted measurement. Select a folder type in a name for the file and dont forget to put .dat or .txt etc. at the end! 
The checkbox below does'nt work anymore!

# Going to the second tab 'Colour-Plotting'

Colour plotting works. You just have to have a C++ dist. installed, see Visual Studio.
Temporarily plotting of fitted or externaly loaded parameters is also placed on this tab. Be carefull while using 'load parameters from file', in order to plot them a measurements has to be opened using ctrl + O!.

One can select the angle to view the fit together with the experimentally obtained points.
By checking the box 'Change values' the parameters of the fit become changeable.
