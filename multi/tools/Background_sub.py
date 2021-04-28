import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets
from lib.BackSub import *
from lib.CustomNodes import *

from pyqtgraph.flowchart import Flowchart
#from pyqtgraph.Qt import QtGui, QtCore
#import pyqtgraph as pg
import numpy as np
#import pyqtgraph.metaarray as metaarray


class Background_sub(QtWidgets.QWidget):
    def __init__(self,data_In=None,*args,**kwargs):
        # data_In will be Adata
        if data_In == None:
            print("Please load data")
            # For testing, create noise data
            data_In = []
            for i in range(50):
                dat = np.random.normal(size=1024)
                dat[200:300] += 1
                dat += np.sin(np.linspace(0, 100, 1024))
                data_In.append(dat)
        data_In = np.array(data_In)
        super().__init__()
        self.ui = Ui_BackgroundSubtraction()
        self.ui.setupUi(self)

        # Change plot labels, because I was too lazy to create a new widget
        pw1 = self.ui.View_Orig
        pw1.plt.setLabel('bottom', 'Magnetic Field', units='mT')  # X-Axis
        pw1.plt.setLabel('left', "Amplitude (Arb. Units)", units='')  # Y-Axis

        pw2 = self.ui.View_Mod
        pw2.plt.setLabel('bottom', 'Magnetic Field', units='mT')  # X-Axis
        pw2.plt.setLabel('left',"Amplitude (Arb. Units)",units='')   # Y-Axis

        pw3 = self.ui.View_Colour

        # Create Flowchart with two IO Nodes
        fc = Flowchart(terminals={
            'dataIn': {'io': 'in'},
            'dataOut': {'io': 'out'}
        })
        w = fc.widget()
        self.ui.verticalLayout.addWidget(w)

        # Create metaarray using data_In and some information
        #data = metaarray.MetaArray(data_In, info=[{'name': 'Amplitude data', 'values': data_In}, {}])

        fc.setInput(dataIn=data_In) #Set data to Input Node

        ## If we want to make our custom node available only to this flowchart,
        ## then instead of registering the node type globally, we can create a new
        ## NodeLibrary:
        library = fclib.LIBRARY.copy()  # start with the default node set
        # Add the node to two locations in the menu to demonstrate
        # that we can create arbitrary menu structures
        library.addNodeType(SavGol_Smooth, [('FMR', 'Filter')])
        library.addNodeType(FMR_Subtract_Average, [('FMR', 'Filter')])
        library.addNodeType(FMR_Subtract_Average_Colour, [('FMR', 'Filter')])
        library.addNodeType(Measurement_Select, [('FMR', 'Data')])
        library.addNodeType(Average_Ang_Dep, [('FMR', 'Data')])
        library.addNodeType(ImageViewNode, [('Display',)])
        fc.setLibrary(library)

        plotList = {'Top Plot': pw1.plt, 'Bottom Plot': pw2.plt}

        pw1Node = fc.createNode('PlotWidget', pos=(0, -150))
        pw1Node.setPlotList(plotList)
        pw1Node.setPlot(pw1.plt)

        pw2Node = fc.createNode('PlotWidget', pos=(150, -150))
        pw2Node.setPlot(pw2.plt)
        pw2Node.setPlotList(plotList)

        pw3Node = fc.createNode('ImageView', pos=(300, -150))
        pw3Node.setView(pw3.img)

        fNode = fc.createNode('Spectrum select', pos=(0, 0))
        fc.connectTerminals(fc['dataIn'], fNode['dataIn']) # Use this Slice node to select the measurement that you want to smooth aka axis
        fc.connectTerminals(fNode['dataOut'], pw1Node['In'])
        fc.connectTerminals(fNode['dataOut'], pw2Node['In'])
        fc.connectTerminals(fNode['dataOut'], fc['dataOut'])
        #fc.connectTerminals(fc['dataIn'], fc['dataOut'])

        # To obtain your filtered data use:
        #print(fc.output()['dataOut'])

