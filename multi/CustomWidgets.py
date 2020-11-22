import numpy as np
import pyqtgraph as pg
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtWidgets, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas, NavigationToolbar2QT as NavigationToolbar



class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)

class PlotWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                  # Create canvas object
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.vbl = QtWidgets.QVBoxLayout()         # Set box for plotting
        self.vbl.addWidget(self.toolbar)         
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

class GradWidget(QtWidgets.QWidget):
    sigGradientChanged = QtCore.Signal(object)
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.vbl = QtWidgets.QVBoxLayout()
        self.w = pg.GradientWidget()
        self.w.loadPreset('magma') #Magma colormap set as default
        self.vbl.addWidget(self.w)
        self.setLayout(self.vbl)
        self.w.sigGradientChanged.connect(self.sigGradientChanged)  # If GradientWidget is changed this will receive a Signal

    def get_colormap(self):
        return self.w.colorMap()

    def print_colormap(self):
        print(self.w.colorMap())