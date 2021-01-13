# FMR-Fit
Small GUI for extracting parameters of ferromagnetic resonance measurements made in python.
Packages used:
matplotlib, numpy, lmfit (wrapper for scipy), PyQt5 (Basicly PySide2), symengine (C wrapper of Symengine for Python), PyQtGraph 


# Usage
First usage:
      You can download the code pressing the green "Code" button on the GitHub Page. Please unzip the files to your desired path. After unpacking go to  folder "Setup and             Install" and follow the instructions there.
      If everything is setup, you can start the programm via the "Run.bat" file in the main folder. The programm and a console will then start. Errors and generous output will         be presented in the console!
In order to manipulate this code a Python editor is needed. For big programming projects PyCharm is recommended, but if its just a small script sublime text 3 will also work.

# First tab: 'Fit'

At first open a measurement pressing ctrl + O, or anything that requires this file. 
Remember that only measurements made from Brukers Xepr converted to ascii will work for now! Alternatively any .dat file will work, if every column is arranged like Xepr does.

The next step will be to select a dataset (a measurement extracted from the ascii file will be called dataset in the following), by either typing the number or using your scrollwheel. You will see a change in the top right plot. There it is possible to select the viewport. Use the slider on the bottom first, to select the maximum viewing position, then proceed with the left slider to select the minimum. This will also set the fitting range. Faulty spectra can be detected by going through the spectra and watching the data plot. Once a faulty spectrum is detected, it can be left out for fitting, by adding it to Dropped Spectra.

Then select a function for fitting and the desired number of function, for now between 1 and 10.
Button "Plot" will first fit the measurement and plots it in a new window. The button below Plot, 'Set inital parameters to fit parameters', will take the fitted parameters and sets them as the new initial values (start values for fit).

For dynamic plotting it is advisable to set the initial values to fitted parameters using the spectrum 0 (or the first spectrum according to dropped Spectra), to get the best result! By pushing the dynamic fit button, a window is opened asking for a place to save the fitted measurement. Select a folder and type in a name for the file, if the file ending was not specified it will be saved as a .dat file.

The checkbox "Robust Fit" on the far bottom can be checked if the result of the normal dyn-Fit is not acceptable. This will then activate a global minimizer, that will find an optimum but takes more time to fully solve.

# Going to the second tab 'Colour-Plotting'

Colour plotting works and is done using the mayavi libary. You just have to have a C++ dist. installed, see Visual Studio. This is why mayavi is optional, therefore this option is also optional at the moment.

Temporarily plotting of fitted or externaly loaded parameters is also placed on this tab. Be carefull while using 'load parameters from file', in order to plot them a measurements has to be opened using ctrl + O!.

One can select the angle to view the fit together with the experimentally obtained points.
By checking the box 'Change values' the parameters of the fit become changeable. And "Save adjusted Params" will save your adjusted parameters.

# Anisotropy Fit Tabs

There are two algorithms implemented to solve for anisotropy constants. The firt one is written in Mathematica, the latter one in Python. 
On these Tabs you can specify which Free Energy density to use and what the fit parameters are. Also the resoultion of the Simulation/Fit is adjusteable via Anglestep. 
What is special for the Python algorithm is, that you can run a quick preview of your initial Parameters and a possibility to quickly determine the shift, by changeing the shift value. 
Hitting any of the "Submit" buttons will start either Mathematica or Python. 

The difference between these algorithms hardly depends on your experiment. The Mathematica approach in general is way more robust in solving the anisotropy constans and definetly will converge at some point, but will be slower compared to the Python approach.

For example you can run a quick Fit in Python ~ 40sec - 1-2 min (depending on problem,resolution,CPU speed and cores), if it does not converge use Mathematica.

Please think about your results, critically analyze them. A converged solution does not mean it is the best or even a correct solution of your problem! Change the free Energy denisity if neccessary!



# Additional things

If you find a bug please got to the issues page and add an entry there, to explain what happend and where it happend.
Also if you have new suggestions you want to add to the programm, feel free to write an entry in the implementations Todo list at the 'Projects' tab on github.
