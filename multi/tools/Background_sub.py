import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets
from lib.BackSub import *


class Background_sub(QtWidgets.QWidget):
    def __init__(self,*args,**kwargs):
        super().__init__()
        self.ui = Ui_BackgroundSubtraction()
        self.ui.setupUi(self)

        # Change plot labels, because I was too lazy to create a new widget
        self.ui.View_Orig.plt.setLabel('bottom', 'Magnetic Field', units='mT')  # X-Axis
        self.ui.View_Orig.plt.setLabel('left',"Amplitude (Arb. Units)", units='')   # Y-Axis
        self.ui.View_Mod.plt.setLabel('bottom', 'Magnetic Field', units='mT')  # X-Axis
        self.ui.View_Mod.plt.setLabel('left',"Amplitude (Arb. Units)",units='')   # Y-Axis

    def test(self):
        print('GAgwagawg')

