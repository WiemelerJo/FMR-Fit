import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
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

class Popup_View(QtWidgets.QWidget):
    def __init__(self,*args):
        super().__init__()

        Z = np.asarray(args[0])
        self.Z = np.transpose(Z)

        chunk = args[1]
        winkel_max = args[2]

        layout = QtWidgets.QGridLayout()
        #win = pg.GraphicsLayoutWidget()
        win = pg.GraphicsLayoutWidget(show=True)
        view = win.addViewBox(0,1)

        xScale = pg.AxisItem(orientation='bottom', linkView=view)
        win.addItem(xScale, 1, 1)

        yScale = pg.AxisItem(orientation='left',linkView=view)
        win.addItem(yScale,0,0)

        self.label = pg.TextItem(text='Hover Event',anchor=(0,0))
        view.addItem(self.label,ignoreBounds=True)

        xScale.setLabel('Magnetic Field',units='mT')
        yScale.setLabel('Angle',units='Deg')

        img = pg.ImageItem(border='w')
        self.img = img
        data = np.array(self.Z)
        img.setImage(data)
        img.hoverEvent = self.imageHoverEvent

        img.setRect(QtCore.QRect(0, 0, chunk, winkel_max))
        view.addItem(img)

        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        win.addItem(hist,0,2)

        layout.addWidget(win, 0, 1, 3, 1)
        self.setLayout(layout)


    def imageHoverEvent(self, event):
        """Show the position, pixel, and value under the mouse cursor.
        """
        if event.isExit():
            self.label.setText("")
            return
        pos = event.pos()
        i, j = pos.y(), pos.x()
        i = int(np.clip(i, 0, self.Z.shape[0] - 1))
        j = int(np.clip(j, 0, self.Z.shape[1] - 1))
        val = self.Z[i, j]
        ppos = self.img.mapToParent(pos)
        x, y = ppos.x(), ppos.y()
        self.label.setText("Pos (Field,Angle): %0.1f, %0.1f     Amplitude: %g" % (x, y, val))