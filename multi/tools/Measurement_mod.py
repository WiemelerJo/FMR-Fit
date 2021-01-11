import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets
from lib.Modifier import *
from scipy.stats import linregress

class Measurement_Mod(QtWidgets.QWidget):
    def __init__(self,*args,**kwargs):
        super().__init__()
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.ui.button_load_sec_meas.pressed.connect(self.load_second_data)
        self.ui.button_calc.pressed.connect(self.normalise_frequency)
        self.ui.button_save.pressed.connect(self.merge_data)
        #self.ui.lineEdit_freq_1.editingFinished
        #self.ui.lineEdit_freq_2.editingFinished
        self.ui.spinbox_offset.valueChanged.connect(self.offset_event)

        # args: Measurement 1 data
        #   [0] Z
        #   [1] H_range
        #   [2] WinkelMax
        # Let Z_1 be the main measurement, that is fixed
        # Z_2 will be shiftable

        self.Bdata1 = args[0]
        self.Z_1 = args[1] # Adata1
        self.Winkeldata1 = args[2]
        self.winkelstep = self.Winkeldata1[1][0] - self.Winkeldata1[0][0]

        self.H_range_1 = max(self.Bdata1[0]) * 1000
        self.WinkelMax_1 = args[3]

        self.i_min1 = args[4]
        self.i_max1 = args[5]
        self.plot(1)

    def plot(self,order):
        if order == 1:
            self.ui.combined_plot.img.setImage(np.asarray(self.Z_1))
            self.ui.combined_plot.img.setRect(QtCore.QRect(0, 0, self.WinkelMax_1, self.H_range_1))
        else:
            self.ui.combined_plot.img2.setImage(np.asarray(self.Z_2))
            self.ui.combined_plot.img2.setRect(QtCore.QRect(0, self.H_offset_2, self.winkel_diff_2, self.H_diff_2))

    def load_second_data(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')
        if fname[0]:
            try:
                D = np.loadtxt(fname[0],dtype='float') #Load Data
            except:
                D = np.loadtxt(fname[0],dtype='float',skiprows=2)   #Load Data with header file
            counts = np.unique(D[:,2], return_counts=True)  #get the int count of different numbers at row 3
            counts1 = counts[1]
            n = counts1[0]
            chunksize = int(len(D[:,0])/n)  #set chunk size for example 1024,2048,4096

            i2 = 0   # laufvariable = 0
            i_min2 = 0
            i_max2 = int(chunksize)
            self.i_min2 = i_min2
            self.i_max2 = i_max2

            self.Bdata2 = np.split(np.true_divide(np.array(D[:, 1]), 10000), chunksize)  # List of sublists, magnetic field
            self.Z_2 = np.split(np.array(D[:, 3]), chunksize)  # List of sublists, amplitude data, for colour plot

            self.BMax2 = max(self.Bdata2[0])
            self.H_range_2 = max(self.Bdata2[0]) * 1000
            self.H_diff_2 = self.BMax2 * 1000 - min(self.Bdata2[0]) * 1000
            self.H_offset_2 = min(self.Bdata2[0]) * 1000
            self.WinkelMin_2 = min(D[:, 2])
            self.WinkelMax_2 = max(D[:, 2])
            self.winkel_diff_2 = max(D[:, 2]) - min(D[:, 2])
            self.winkelstep2 = D[:, 2][n] - D[:, 2][n-1]
            print(self.winkelstep2)

            self.plot(2)
            self.calculate = True

    def merge_data(self):
        # Todo: VERY IMPORTANT: IF THE ANGLESTEPS OF MEASUREMENTS ONE AND TWO ARE DIFFERENT EG. 0,5 AND 1,0 DEG PER STEP THIS WILL FAIL!!!
        # merge and save data to file similar to bruker ascii

        temp_arr = self.Bdata1[0]
        Bdiff1 = temp_arr[1] - temp_arr[0]
        BMin1 = min(temp_arr)
        BMax1 = max(temp_arr)

        # Erstelle Array das selbe stepsize wie Bdata1 hat aber range von Ende Bdata1 bis Ende Bdata2
        # (Bdata2 sollte hier bereits Frequenz angepasst sein)
        Add_B_array = np.arange(BMax1, self.BMax2, Bdiff1)

        Bdata_fin = []
        Adata_fin = []
        Winkeldata_fin = []

        offset = self.ui.spinbox_offset.value()
        if self.winkelstep != self.winkelstep2:
            print('Winkelssteps different!')
            print(self.winkelstep,self.winkelstep2)
        self.winkelarray = np.arange(offset, offset + self.winkel_diff_2 + self.winkelstep, self.winkelstep)

        A_extrapol = np.vectorize(self.extrapol)

        for i in range(self.i_min1,self.i_max1):
            Bdata_fin.append(np.append(self.Bdata1[i], Add_B_array)) # Add to Bdata1
            Winkeldata_fin.append(np.append(self.Winkeldata1[i], np.full((1, len(Add_B_array)), self.Winkeldata1[i][0]))) # Add to Winkeldata1

            # Add to Adata1
            # A_extrapol = interp1d(self.Bdata1[i], self.Z_1[i], fill_value='extrapolate')
            # Adata_fin.append(np.append(self.Z_1[i],A_extrapol(Add_B_array)))

            slope, intercept, r, p, se = linregress(self.Bdata1[i],self.Z_1[i])
            Adata_fin.append(np.append(self.Z_1[i], A_extrapol(Add_B_array, slope, intercept)))
            #print(len(self.winkelarray),len(self.Z_2))
            if self.Winkeldata1[i][0] in self.winkelarray:
                index = np.where(self.winkelarray == self.Winkeldata1[i][0])[0]
                for k in range(len(self.Z_2[index[0]])):
                    Adata_fin[i][-k] = self.Z_2[index[0]][-k]

        Bdata_fin = np.asarray(Bdata_fin).flatten() * 10000
        Adata_fin = np.asarray(Adata_fin).flatten()
        Winkeldata_fin = np.asarray(Winkeldata_fin).flatten()

        # index = np.linspace(1,Bdata_fin.shape[0],num=Bdata_fin.shape[0])
        index = np.full(len(Bdata_fin),420)
        #print(Bdata_fin.shape[0],index.shape, Bdata_fin.shape, Adata_fin.shape, Winkeldata_fin.shape)
        save_array = [index,Bdata_fin, Winkeldata_fin, Adata_fin]
        np.savetxt('test.dat',np.transpose(save_array) ,delimiter='\t', newline='\n')

    def extrapol(self, B, slope, offset):
        return B * slope + offset

    def offset_event(self,d):
        # Offset Image by changeing QRect x value
        # QRect(x, y, Breite, Höhe)
        self.ui.combined_plot.img2.setRect(QtCore.QRect(d, self.H_offset_2, self.winkel_diff_2, self.H_diff_2))

    def normalise_frequency(self):
        # Use one of the two measurements and normalize it to the frequency of the other measurement
        if self.calculate:
            freq1 = float(self.ui.lineEdit_freq_1.text())
            freq2 = float(self.ui.lineEdit_freq_2.text())

            self.Bdata2 = np.asarray(self.Bdata2)
            self.Bdata2 *= freq2/freq1 # Resonance position is linear in frequency

            self.H_range_2 = max(self.Bdata2[0]) * 1000
            self.H_diff_2 = max(self.Bdata2[0]) * 1000 - min(self.Bdata2[0]) * 1000
            self.H_offset_2 = min(self.Bdata2[0]) * 1000
            self.plot(2)
            self.calculate = False