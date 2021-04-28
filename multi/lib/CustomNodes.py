from pyqtgraph.flowchart import Flowchart, Node
import pyqtgraph.flowchart.library as fclib
from pyqtgraph.flowchart.library.common import CtrlNode
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import numpy as np

from scipy.signal import savgol_filter




## At this point, we need some custom Node classes since those provided in the library
## are not sufficient. Each node will define a set of input/output terminals, a 
## processing function, and optionally a control widget (to be displayed in the 
## flowchart control panel)
class ImageViewNode(Node):
    """Node that displays image data in an ImageView widget"""
    nodeName = 'ImageView'
    
    def __init__(self, name):
        self.view = None
        ## Initialize node with only a single input terminal
        Node.__init__(self, name, terminals={'data': {'io':'in'}})
        
    def setView(self, view):  ## setView must be called by the program
        self.view = view
        
    def process(self, data, display=True):
        ## if process is called with display=False, then the flowchart is being operated
        ## in batch processing mode, so we should skip displaying to improve performance.
        
        if display and self.view is not None:
            ## the 'data' argument is the value given to the 'data' terminal
            if data is None:
                self.view.setImage(np.zeros((1,1))) # give a blank array to clear the view
            else:
                self.view.setImage(data)



## We will define a savgol filter node as a subclass of CtrlNode.
## CtrlNode is just a convenience class that automatically creates its
## control widget based on a simple data structure.
class SavGol_Smooth(CtrlNode):
    """Return the input data passed through a filter."""
    nodeName = "Savitzky-Golay"
    uiTemplate = [
        ('range',  'spin', {'value': 23.0, 'step': 1.0, 'bounds': [0.0, None]}),
        ('polynomial order', 'spin', {'value': 2.0, 'step': 1.0,'bounds': [1.0 , None]}),
    ]

    def __init__(self, name):
        ## Define the input / output terminals available on this node
        terminals = {
            'dataIn': dict(io='in'),    # each terminal needs at least a name and
            'dataOut': dict(io='out'),  # to specify whether it is input or output
        }                              # other more advanced options are available
                                       # as well..
        
        CtrlNode.__init__(self, name, terminals=terminals)
        
    def process(self, dataIn, display=True):
        # CtrlNode has created self.ctrls, which is a dict containing {ctrlName: widget}
        range_savgol = int(self.ctrls['range'].value())
        # Make sure range is uneven:
        if (range_savgol % 2) == 0:
            # range is even, so make it uneven
            range_savgol += 1

        polynomial_order = int(self.ctrls['polynomial order'].value())

        output = savgol_filter(dataIn, range_savgol, polynomial_order)
        return {'dataOut': output} # Returns 1dim list


class Measurement_Select(CtrlNode):
    nodeName = "Spectrum select"
    uiTemplate = [('Index/Spectra',  'spin', {'value': 0.0, 'step': 1.0, 'bounds': [0.0, None]})]

    def __init__(self, name):
        ## Define the input / output terminals available on this node
        terminals = {
            'dataIn': dict(io='in'),    # each terminal needs at least a name and
            'dataOut': dict(io='out'),  # to specify whether it is input or output
        }                              # other more advanced options are available
                                       # as well..
        CtrlNode.__init__(self, name, terminals=terminals)
        
    def process(self, dataIn, display=True):
        # CtrlNode has created self.ctrls, which is a dict containing {ctrlName: widget}
        index = int(self.ctrls['Index/Spectra'].value())
        output = dataIn[index]
        return {'dataOut': output} # Returns 1dim list


class Average_Ang_Dep(CtrlNode):
    # Average everything to one single averaged spectra
    nodeName = "Average angular dependence"

    def __init__(self, name):
        ## Define the input / output terminals available on this node
        terminals = {
            'dataIn': dict(io='in'),    # each terminal needs at least a name and
            'dataOut': dict(io='out'),  # to specify whether it is input or output
        }                              # other more advanced options are available
                                       # as well..
        CtrlNode.__init__(self, name, terminals=terminals)
        
    def process(self, dataIn, display=True):
        # CtrlNode has created self.ctrls, which is a dict containing {ctrlName: widget}
        transposed = np.transpose(dataIn)
        averaged = []
        for index, val in enumerate(transposed):
            averaged.append(np.mean(val))

        output = np.asarray(averaged)
        return {'dataOut': output} # Returns 1dim list


class FMR_Subtract_Average(CtrlNode):
    # Average everything to one single averaged spectra and subtract it from another spectra
    nodeName = "FMR_Subtract"
    uiTemplate = [('Index/Spectra',  'spin', {'value': 0.0, 'step': 1.0, 'bounds': [0.0, None]})]

    def __init__(self, name):
        ## Define the input / output terminals available on this node
        terminals = {
            'dataIn': dict(io='in'),    # each terminal needs at least a name and
            'dataOut': dict(io='out'),  # to specify whether it is input or output
        }                              # other more advanced options are available
                                       # as well..
        CtrlNode.__init__(self, name, terminals=terminals)
        
    def process(self, dataIn, display=True):
        # CtrlNode has created self.ctrls, which is a dict containing {ctrlName: widget}
        index = int(self.ctrls['Index/Spectra'].value())
        spectra = dataIn[index]
        transposed = np.transpose(dataIn)
        averaged = []

        for index, val in enumerate(transposed):
            averaged.append(np.mean(val))

        output = np.asarray(spectra) - np.asarray(averaged)
        return {'dataOut': output} # Returns 1dim list


class FMR_Subtract_Average_Colour(CtrlNode):
    # Average everything to one single averaged spectra and subtract it from whole dependency
    nodeName = "FMR Subtract 2D"

    def __init__(self, name):
        ## Define the input / output terminals available on this node
        terminals = {
            'dataIn': dict(io='in'),    # each terminal needs at least a name and
            'dataOut': dict(io='out'),  # to specify whether it is input or output
        }                              # other more advanced options are available
                                       # as well..
        CtrlNode.__init__(self, name, terminals=terminals)
        
    def process(self, dataIn, display=True):
        # CtrlNode has created self.ctrls, which is a dict containing {ctrlName: widget}
        transposed = np.transpose(dataIn)
        averaged = []
        output = []

        for index, val in enumerate(transposed):
            averaged.append(np.mean(val))
        averaged = np.asarray(averaged)

        for array in dataIn:
            output.append( np.asarray(array) - averaged )

        return {'dataOut': np.array(output)} # Returns ndim list, for colourplot