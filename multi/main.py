import sys
import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import math as m
import threading
from datetime import datetime
from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from PyQt5.QtWidgets import QMainWindow, QApplication, QAction,QFileDialog,QProgressBar,QDoubleSpinBox, QCheckBox
from PyQt5.QtCore import QThread, pyqtSignal
from Fitprogramm import *
from fbs_runtime.application_context.PyQt5 import ApplicationContext
#from matplotlib.colors import BoundaryNorm
#from matplotlib.ticker import MaxNLocator
from fitting import Fit
#from fitting import Worker
from parameter_plot import Params_Plot

'''def dyson_func(B, alpha, dB, R, A, slope, offset):
        return ((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))+slope*B+offset

def lorentz_func(B,dB, R, A, slope, offset):
    #return(-(A*(2*(B-dB)))/(m.pi*R**3*((B-dB**2)/(R**2)+1)**2))+slope*B+offset
    return(-64*A*(B-R)/(9*dB*(1+(B-R)**2)/(dB*m.sqrt(3)/4)**2)**2+slope*B+offset)'''


data_value = 'unloaded'
fit_value = 'unfitted'
dyn_value = 'single'
dyn_fit_value = 'unfitted'
parameter_value = 'unloaded'
p = 0
j = 350
j_min = 0
colormap = 'inferno'
tick_number = 256 #Max = 10000, bereits 1000 erstellt 200 MB file!
excep_append = []
exceptions = []
fit_num = 1

parameter_Table_names_L = ['dB','R','A']
parameter_table_names_L_final = []
parameter_Table_names_D = ['alpha','dB','R','A']
parameter_table_names_D_final = []

#Arrays to set the default values in SpinBox. #For linear: [slope, offset]
                                             #For L: [dB1, R1, A1, dB2, R2, A2, ...]
                                             #For D: [alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, ....]
default_linear = [0.5,1]
default_values_L = [0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7]
default_values_D = [0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7]

#Arrays for default Boundaries. For Linear : [slope_min,slope_max, offset_min,offset_max] 
                               #For L: [dB1_min,dB1_max, R1_min,R1_max, ...]
                               #For D: [alpha1_min,alpha1_max, dB1_min,dB1_max ...]
default_boundaries_linear_min = [-5, -5]
default_boundaries_linear_max = [5, 5]
default_boundaries_L_min = [0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0]
default_boundaries_L_max = [0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15]
default_boundaries_D_min = [0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001]
default_boundaries_D_max = [1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15]

#Arrays for step size for spinbox Linear: [slope, offset]
                                 #L: [dB, R, A]
                                 #D: [alpha,dB, R, A]
default_stepsize_linear = [0.1, 0.1]
default_stepsize_L = [0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1]
default_stepsize_D = [0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1]


class Worker(QThread):
    #This class is responsible for the dynamic fit by creating a so called worker to do the job
    i_signal = pyqtSignal() # signal for the processBar
    error_dynfit_signal = pyqtSignal() # error signal
    def __init__(self,index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,dyn_value,i_min,i_max,j,j_min,exceptions):
        QThread.__init__(self)

    def fit(self,l,params_table,model):
        Bdata2 = Bdata[l]
        Adata2 = Adata[l]
        result_dyn = model.fit(Adata2[j_min:j], params_table, B=Bdata2[j_min:j])
        return result_dyn

    def run(self):
        global dyn_fit_value
        names = []
        temp_paras = []
        params_table = Fit(index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).give_param_table(index_model)
        model = Fit(index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).give_model(index_model)

        #clean up, just in case
        names.clear()
        temp_paras.clear()

        for name, param in result.params.items():
            names.append('{}'.format(name))  #create an array of names for the header
            temp_paras.append(float(param.value))  #create the first parameters to save into the .dat file
        temp_paras.append(float(Winkeldata[0][0]))

        Parameter_list = np.zeros(shape=(i_max,len(temp_paras)))
        for l in range(i_min,i_max):
            if l not in exceptions:
                self.i_signal.emit() # update progressbar
                temp_paras.clear() # clear temp_paras for each iteration
                temp_result = self.fit(l,params_table,model) #this is the fit
                for name, param in temp_result.params.items():
                    temp_paras.append(float(param.value))  
                    params_table[name].set(value=param.value,min=None,max=None)
                temp_paras.append(float(Winkeldata[l][0]))
                Parameter_list[l] = temp_paras
        now = datetime.now()
        np.savetxt(fileName,Parameter_list,delimiter='\t',newline='\n', header='This data was fitted {}  \nDropped Points {}     \nData is arranged as follows {}'.format(now.strftime("%d/%m/%Y, %H:%M:%S"),exceptions,names))
        dyn_fit_value = 'fitted'

class ColourPlot(QThread):
    def __init__(self,X,Y,Z,colormap,tick_number,D_min,D_max,exceptions):
        QThread.__init__(self)
        print('Starting')

    def run(self):
        print('Starting to Plot')
        print(tick_number)
        #numrecs = np.size(D,axis=0)
        '''
        Alter Plot mit Matplotlib (Langsam und eher schlecht)
        x = Bdata
        y = Winkeldata
        z = Adata
        mp.ticker.Locator.MAXTICKS = 10000
        levels = MaxNLocator(nbins=tick_number).tick_values(D_min,D_max)
        cmap = plt.get_cmap(colormap)
        cf = plt.contourf(x, y, z, levels=levels,cmap=cmap)
        cbar = plt.colorbar(cf)
        MyForm.colour_plot_main(self,cbar)'''
        '''cbar.ax.get_yaxis().set_ticks([])
        for j, lab in range(0,numrecs):
            cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('# of contacts', rotation=270)'''
        #------------------------------------------------------------------------
        #Mayavi Plot in 3D (unterstützt OpenGL, ist daher schneller; mehr details)




class MyForm(QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.actionOpen.triggered.connect(self.openFileDialog)
        self.ui.actionSave.triggered.connect(self.saveFileDialog)
        self.ui.actionExit.triggered.connect(self.Exit)
        self.ui.Button_Plot.clicked.connect(self.plot)
        self.ui.Button_set_fit_params.clicked.connect(self.set_fit_params)
        self.ui.comboBox_fit_model.currentIndexChanged.connect(self.model_type)
        self.ui.comboBox_Fit_num.currentIndexChanged.connect(self.select_fit_number)
        self.ui.Button_dyn_fit.clicked.connect(self.dyn_fit)
        self.ui.colour_button.clicked.connect(self.plot_in_colour)
        self.ui.Drop_Scalebar.valueChanged.connect(self.data_range)
        self.ui.min_drop_Bar.valueChanged.connect(self.data_range)
        self.ui.select_datanumber.valueChanged.connect(self.set_datanumber)
        self.ui.Scroll_Bar_dropped_points.valueChanged.connect(self.set_dataslider)
        self.ui.Dropped_points_edit.editingFinished.connect(self.Exceptions)
        self.ui.Button_dropped_points.clicked.connect(self.append_dropped_points)
        self.ui.Dropped_points_edit.textEdited.connect(self.reset_drop)
        self.ui.plot_parameter.clicked.connect(self.parameter_plot)
        self.ui.load_params_from_file.clicked.connect(self.load_parameter_plot)
        self.ui.parameter_data_select.valueChanged.connect(self.change_parameter_angle)
        self.ui.param_change_alpha.valueChanged.connect(self.changeing_parameters)
        self.ui.param_change_R.valueChanged.connect(self.changeing_parameters)
        self.ui.param_change_dB.valueChanged.connect(self.changeing_parameters)
        self.ui.param_change_A.valueChanged.connect(self.changeing_parameters)
        self.ui.param_change_slope.valueChanged.connect(self.changeing_parameters)
        self.ui.param_change_offset.valueChanged.connect(self.changeing_parameters)
        self.ui.pushButton.clicked.connect(self.test_table)
        self.show()

    def reset_drop(self):
        #empties the exception array
        global excep_append
        if self.ui.Dropped_points_edit.text() == '':
            excep_append.clear()
        else:
            x = str(self.ui.Dropped_points_edit.text())
            excep_append = list(map(int, x.split(',')))
   
    def append_dropped_points(self):
        #as something is added to the exceptions, this is beeing called
        global excep_append
        excep_append.append(str(self.ui.Scroll_Bar_dropped_points.value()))
        text = ','.join(map(str, excep_append))
        self.ui.Dropped_points_edit.setText(str(text))

    def Exceptions(self):
        # takes the string input from the exception text field inside the GUI and seperates the entries into a list
        global exceptions
        global excep_append
        print('Hallo')
        x = str(self.ui.Dropped_points_edit.text())
        exceptions = list(map(int, x.split(',')))

    def change_parameter_angle(self):
        # has to be reworked 
        global W
        W = self.ui.parameter_data_select.value()
        self.select_fitted_params()

    def set_dataslider(self):
        self.ui.select_datanumber.setValue(int(self.ui.Scroll_Bar_dropped_points.value()))

    def set_datanumber(self):
        #Defines the dataset looked at
        global i
        global Bdata2
        global Adata2
        if data_value == 'loaded':

            if int(self.ui.select_datanumber.text()) > len(Adata)-1:
                i = len(Adata)-1
                self.ui.select_datanumber.setValue(int(i))
            else:
                i = int(self.ui.select_datanumber.value())
            Bdata2 = Bdata[i]
            Adata2 = Adata[i]
            self.plot_data()
        else:
            self.openFileDialog()

    def colour_plot_main(self,cbar):
        #unused at the moment. However it should create the contour plot
        cbar.set_label('Amplitude [Arb. U.]')
        plt.xlabel('Magnetic Field [T]')
        plt.ylabel('Angle [Degrees]')
        plt.show()
        plt.close()

    def update_bar(self):
        #updates the progressbar
        global p
        self.ui.label_7.setText(str(p))
        p +=1
        self.ui.progressBar.setValue(p)

    def dyn_fit(self):
        #starts the dynamic fitting routine
        global dyn_value
        dyn_value = 'dynamic'
        self.saveFileDialog()

    def model_type(self,Bool):
        #selects the type of function used in the fit
        global index_model
        index_model = self.ui.comboBox_fit_model.currentIndex()


    def select_fit_number(self):
        #define the number of function that are beeing fitted
        global fit_num
        fit_num = self.ui.comboBox_Fit_num.currentIndex()
        self.make_parameter_table()

    def make_parameter_table(self):
        global parameter_table_names_D_final
        global parameter_table_names_L_final
        #create the table to edit the parameters
        #uses for loops to create the rows

        #cleaning up
        self.ui.Parameter_table.clearContents()
        parameter_table_names_L_final = []
        parameter_table_names_D_final = []

        #function that catches errors!
        try:
            if index_model == 2:
                #table for lorentz
                self.ui.Parameter_table.setRowCount(fit_num*3+2)
                self.ui.Parameter_table.setCellWidget(0, 0, QDoubleSpinBox()) # initial_values slope
                self.ui.Parameter_table.setCellWidget(0, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(0, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(0, 3, QCheckBox()) # use for fitting ?  

                self.ui.Parameter_table.setCellWidget(1, 0, QDoubleSpinBox()) # initial_values offset
                self.ui.Parameter_table.setCellWidget(1, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(1, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(1, 3, QCheckBox()) # use for fitting ? 

                parameter_table_names_L_final.append('slope')
                parameter_table_names_L_final.append('offset')   

                for index in range(1,fit_num+1):
                    for list_index in parameter_Table_names_L:
                        parameter_table_names_L_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_L_final)

                for zähler in range(0,fit_num*3):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox()) # initial_values
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox()) # bound_min
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox()) # bound_max
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox()) # use for fitting ? 

                self.set_default_values()

            else:
                #table for dyson (at the moment!)
                self.ui.Parameter_table.setRowCount(fit_num*4+2)
                self.ui.Parameter_table.setCellWidget(0, 0, QDoubleSpinBox()) # initial_values
                self.ui.Parameter_table.setCellWidget(0, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(0, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(0, 3, QCheckBox()) # use for fitting ?  

                self.ui.Parameter_table.setCellWidget(1, 0, QDoubleSpinBox()) # initial_values offset
                self.ui.Parameter_table.setCellWidget(1, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(1, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(1, 3, QCheckBox()) # use for fitting ? 

                parameter_table_names_D_final.append('slope')
                parameter_table_names_D_final.append('offset')  

                for index in range(1,fit_num+1):
                    for list_index in parameter_Table_names_D:
                        parameter_table_names_D_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_D_final)

                for zähler in range(0,fit_num*4):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox())

                self.set_default_values()
        except:
                print('Error in make_parameter_table')
                #Assuming a lorentz func as long as index_model hasnt been set
                self.ui.Parameter_table.setRowCount(fit_num*3+2)
                self.ui.Parameter_table.setCellWidget(0, 0, QDoubleSpinBox()) # initial_values slope
                self.ui.Parameter_table.setCellWidget(0, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(0, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(0, 3, QCheckBox()) # use for fitting ?  

                self.ui.Parameter_table.setCellWidget(1, 0, QDoubleSpinBox()) # initial_values offset
                self.ui.Parameter_table.setCellWidget(1, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(1, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(1, 3, QCheckBox()) # use for fitting ? 

                parameter_table_names_L_final.append('slope')
                parameter_table_names_L_final.append('offset')   

                for index in range(1,fit_num+1):
                    for list_index in parameter_Table_names_L:
                        parameter_table_names_L_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_L_final)

                for zähler in range(0,fit_num*3):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox()) # initial_values
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox()) # bound_min
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox()) # bound_max
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox()) # use for fitting ? 

                self.set_default_values()

    def test_table(self):
        for name, param in result.params.items():
            print(float(param.value))


    def set_default_values(self):
        #sets default values into the spinbox according to arrays defined in the beginning
        #it basicly is just a for loop in order to cath every row of the table according to the numbver of lines
        #try:
        if index_model == 2:
            for zahl in range(0,2):
                self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-100)
                self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-100)
                self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-100)

                self.ui.Parameter_table.cellWidget(zahl,0).setSingleStep(default_stepsize_linear[zahl])
                self.ui.Parameter_table.cellWidget(zahl,1).setSingleStep(default_stepsize_linear[zahl])
                self.ui.Parameter_table.cellWidget(zahl,2).setSingleStep(default_stepsize_linear[zahl]) 

                self.ui.Parameter_table.cellWidget(zahl,0).setValue(default_linear[zahl]) # Inital value
                self.ui.Parameter_table.cellWidget(zahl,1).setValue(default_boundaries_linear_min[zahl]) # Boundary minimum
                self.ui.Parameter_table.cellWidget(zahl,2).setValue(default_boundaries_linear_max[zahl]) # Boundary maximum

            for zähler in range(0,fit_num*3):
                self.ui.Parameter_table.cellWidget(zähler+2,0).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zähler+2,1).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zähler+2,2).setDecimals(4)

                self.ui.Parameter_table.cellWidget(zähler+2,0).setSingleStep(default_stepsize_L[zähler])
                self.ui.Parameter_table.cellWidget(zähler+2,1).setSingleStep(default_stepsize_L[zähler])
                self.ui.Parameter_table.cellWidget(zähler+2,2).setSingleStep(default_stepsize_L[zähler]) 

                self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(default_values_L[zähler]) # Inital value
                self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(default_boundaries_L_min[zähler]) # Boundary minimum
                self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(default_boundaries_L_max[zähler]) # Boundary maximum

        else:
            for zahl in range(0,2):
                self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-100)
                self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-100)
                self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-100)

                self.ui.Parameter_table.cellWidget(zahl,0).setSingleStep(default_stepsize_linear[zahl])
                self.ui.Parameter_table.cellWidget(zahl,1).setSingleStep(default_stepsize_linear[zahl])
                self.ui.Parameter_table.cellWidget(zahl,2).setSingleStep(default_stepsize_linear[zahl]) 

                self.ui.Parameter_table.cellWidget(zahl,0).setValue(default_linear[zahl]) # Inital value
                self.ui.Parameter_table.cellWidget(zahl,1).setValue(default_boundaries_linear_min[zahl]) # Boundary minimum
                self.ui.Parameter_table.cellWidget(zahl,2).setValue(default_boundaries_linear_max[zahl]) # Boundary maximum

            for zähler in range(0,fit_num*4):
                self.ui.Parameter_table.cellWidget(zähler+2,0).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zähler+2,1).setDecimals(4)
                self.ui.Parameter_table.cellWidget(zähler+2,2).setDecimals(4)

                self.ui.Parameter_table.cellWidget(zähler+2,0).setSingleStep(default_stepsize_D[zähler])
                self.ui.Parameter_table.cellWidget(zähler+2,1).setSingleStep(default_stepsize_D[zähler])
                self.ui.Parameter_table.cellWidget(zähler+2,2).setSingleStep(default_stepsize_D[zähler]) 

                self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(default_values_D[zähler]) # Inital value
                self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(default_boundaries_D_min[zähler]) # Boundary minimum
                self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(default_boundaries_D_max[zähler]) # Boundary maximum
        #except:
         #   print('Error in set_default_values')



    def openFileDialog(self):
        global data_value
        global i
        global Bdata
        global Adata
        global Winkeldata
        global n
        global i_min
        global i_max
        global D_min
        global D_max
        global chunksize
        global X
        global Y
        global Z
        fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if fname[0]:
            try:
                D = np.loadtxt(fname[0],dtype='float') #Load Data 
            except:
                D = np.loadtxt(fname[0],dtype='float',skiprows=2)   #Load Data with header file
            counts = np.unique(D[:,2], return_counts=True)  #get the int count of different numbers at row 3
            counts1 = counts[1]
            n = counts1[0]
            chunksize = int(len(D[:,0])/n)  #set chunk size for example 1024,2048,4096
            self.ui.Drop_Scalebar.setMaximum(n-1)   #sets the maximum of the scroll bar for dropped points
            i = 0   # laufvariable = 0
            i_min = 0
            i_max = int(chunksize)
            D_min = min(D[:,3]) #for old colourplot
            D_max = max(D[:,3]) #for old colourplot

            Bdata = np.split(np.true_divide(np.array(D[:,1]),10000),chunksize)  #List of sublists, magnetic field
            Winkeldata = np.split(np.array(D[:,2]),chunksize)   #List of sublists, angle data
            Adata = np.split(np.array(D[:,3]),chunksize)    #List of sublists, amplitude data
            X = np.split(np.true_divide(np.array(D[:,1]),10),chunksize) #List of sublists, magnetic field, for colour plot
            Y = np.split(np.array(D[:,2]),chunksize)    #List of sublists, angle data, for colour plot
            Z = np.split(np.array(D[:,3]),chunksize)     #List of sublists, amplitude data, for colour plot
            if int(self.ui.select_datanumber.text()) > len(Adata)-1: #catches errors, slider for datanumber
                i = len(Adata)-1
                self.ui.select_datanumber.setText(str(i))
            else:
                i = int(self.ui.select_datanumber.text())
            self.ui.select_datanumber.setMaximum(i_max)
            self.ui.label_show_load.setText('Data loaded!')
            self.ui.progressBar.setMaximum(i_max)
            self.ui.Scroll_Bar_dropped_points.setMaximum(i_max)
            self.ui.parameter_data_select.setMaximum(i_max)
            data_value = 'loaded'
            self.set_datanumber()
    
    def data_range(self):
        global j
        global j_min
        global Bdata2
        global Adata2
        if data_value == 'loaded':
            '''Bdata2 = Bdata[i]
            Adata2 = Adata[i]'''
            j = self.ui.Drop_Scalebar.value()
            j_min = self.ui.min_drop_Bar.value()
            self.ui.min_drop_Bar.setMaximum(j) 
            self.ui.Drop_label.setText(str(j))
            self.ui.min_drop_value.setText(str(j_min))
            self.plot_data()
        else:
            self.openFileDialog()
        
    def error_msg(self):
        self.ui.label_params_output.setText('No data fitted yet!')
        print('Fehler irgendwo, bestimmt im dyn fit')

    def parameter_plot(self):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        if dyn_fit_value == 'fitted':
            Params_Plot(Parameter_list,'eigen')
        elif parameter_value == 'loaded':
            print('Unfitted parameters. Using from file instead!')
            Params_Plot(Parameter_from_text,'eigen')
        else:
            self.load_parameters_to_text()
            Params_Plot(Parameter_from_text,'eigen')

    def load_parameter_plot(self):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        print('Loading from file!')
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            Params_Plot(params_fname,'load')

    def load_parameters_to_text(self):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        global Parameter_from_text
        global parameter_value
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            Parameter_from_text = np.loadtxt(params_fname[0],dtype='float',skiprows=1)
        parameter_value = 'loaded'


    def select_fitted_params(self):
        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        global Parameter_from_text
        global parameter_value
        global data_value
        global dyn_fit_value

        if data_value == 'loaded':
            x_temp = Bdata[W]
            y_temp = Adata[W]
            x = x_temp[j_min:j]
            y = y_temp[j_min:j]

            if dyn_fit_value == 'fitted':
                Para = np.array(Parameter_list)
            elif  parameter_value == 'loaded':
                Para = Parameter_from_text
            else:
                self.load_parameters_to_text()
                Para = Parameter_from_text
            alpha_list = Para[:,0]
            db_list = Para[:,1]
            R_list = Para[:,2]
            A_list = Para[:,3]
            slope_list = Para[:,4]
            offset_list = Para[:,5]
            Winkeldata_list = Para[:,6]   
            self.plot_fitted_params(x,y,W,alpha_list,db_list,R_list,A_list,slope_list,offset_list)
            self.ui.param_change_alpha.setValue(alpha_list[W])
            self.ui.param_change_dB.setValue(db_list[W])            
            self.ui.param_change_R.setValue(R_list[W])
            self.ui.param_change_A.setValue(A_list[W])
            self.ui.param_change_slope.setValue(slope_list[W])
            self.ui.param_change_offset.setValue(offset_list[W])

        else:
            self.openFileDialog()
            print('Data is not yet assigned!!')

    def changeing_parameters(self):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        global Parameter_from_text
        global parameter_value
        global data_value
        global dyn_fit_value
        if self.ui.checkBox_change_values.isChecked()==True:
            if data_value == 'loaded':
                x_temp = Bdata[W]
                y_temp = Adata[W]
                x = x_temp[j_min:j]
                y = y_temp[j_min:j]

                if dyn_fit_value == 'fitted':
                    Para = np.array(Parameter_list)
                elif  parameter_value == 'loaded':
                    Para = Parameter_from_text
                else:
                    print('No data loaded!')

                alpha_list = Para[:,0]
                db_list = Para[:,1]
                R_list = Para[:,2]
                A_list = Para[:,3]
                slope_list = Para[:,4]
                offset_list = Para[:,5]
                Winkeldata_list = Para[:,6] 

                alpha_list[W] = self.ui.param_change_alpha.value()
                db_list[W] = self.ui.param_change_dB.value()
                R_list[W] = self.ui.param_change_R.value()
                A_list[W] = self.ui.param_change_A.value()
                slope_list[W] = self.ui.param_change_slope.value()
                offset_list[W] = self.ui.param_change_offset.value()
                self.plot_fitted_params(x,y,W,alpha_list,db_list,R_list,A_list,slope_list,offset_list)    

            else:       
                self.openFileDialog()
                print('Data is not yet assigned!!')



    def plot_fitted_params(self,x,y,W,alpha_list,db_list,R_list,A_list,slope_list,offset_list):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        self.ui.parameter_plot_widget.canvas.ax.clear()
        self.ui.parameter_plot_widget.canvas.ax.set_xlabel('Magnetic Field [T]')
        self.ui.parameter_plot_widget.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
        self.ui.parameter_plot_widget.canvas.ax.plot(x, y,'o', label='Experimental data')
        plotdys = dyson_func(x,alpha_list[W],db_list[W],R_list[W],A_list[W],slope_list[W],offset_list[W])
        self.ui.parameter_plot_widget.canvas.ax.plot(x,plotdys,'r--',label='Dyson')
        self.ui.parameter_plot_widget.canvas.ax.legend()
        self.ui.parameter_plot_widget.canvas.draw() 


    def plot(self):
        global fit_value
        global dyn_value
        global result
        dyn_value = 'single'
        self.model_type(False)
        if data_value == 'loaded':
            self.set_init_params()  #gets the initial parameters from the GUI
            result = Fit(index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).fit(index_model,Adata2,Bdata2, j_min,j) #fit and save it as result
            fit_value = 'fitted'
            self.ui.progressBar.setMaximum(i_max-len(exceptions))
            #self.ui.label_params_output.setText(self.fit(i).fit_report())
            plt.ion()   #enable user actions in plot
            plt.cla()   #plot clean up
            plt.clf()   #plot clean up
            plt.plot(Bdata2[j_min:j],Adata2[j_min:j], '-b', label='Experimental data')  #plot original data
            if index_model == 2:
                self.ui.label_params_output.setText(result.fit_report())    #print fit report Lorentz
                plt.plot(Bdata2[j_min:j],result.best_fit, 'r--', label='Best fit Lorentz')  #plot fit data Lorentz
            else:
                self.ui.label_params_output.setText(result.fit_report())    #print fit report Dyson
                plt.plot(Bdata2[j_min:j],result.best_fit, 'r--', label='Best fit Dyson')    #plot fit data Dyson
            plt.xlabel('Magnetic Field [T]',fontsize=12)
            plt.ylabel('Amplitude [Arb. U.]',fontsize=12)
            plt.legend(fontsize=12)
            plt.show()
        else:
            self.openFileDialog()
        if index_model == 0:
            self.ui.label_params_output.setText('Please select a fitting model first')

    def plot_data(self):
        #Plots data into a matplotlib canvas created in "pyqtgraph" skript, its the vieport for the measurement
        self.ui.plotView.canvas.ax.clear()
        x=Bdata2[j_min:j]
        y=Adata2[j_min:j]
        self.ui.plotView.canvas.ax.set_xlabel('Magnetic Field [T]')
        self.ui.plotView.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
        self.ui.plotView.canvas.ax.plot(x, y)
        self.ui.plotView.canvas.draw()

    def plot_in_colour(self):
        #the most recent colour plotting routine
        from mayavi import mlab
        global tick_number
        global colormap
        tick_number = self.ui.colour_tick_edit.text() #not needed for mayavi
        colormap = self.ui.colour_map_edit.text()   #can also be changed in mayavi UI
        mlab.mesh(X,Y,Z, colormap=colormap)
        mlab.xlabel('Magnetic Field')
        mlab.ylabel('Angle')
        mlab.colorbar(title='Amplitude [Arb. Units]')
        mlab.show()
        '''self.get_thread2 = ColourPlot(X,Y,Z,colormap,tick_number,D_min,D_max,exceptions)
        self.get_thread2.start()'''

    def Exit(self):
        sys.exit()  #Obviously not the start

    def saveFileDialog(self):
        #started when Dyn Fit button pressed
        #gets the filename for saving
        #then starts the worker responsible for dyn fit
        global p
        global fileName
        p=0
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.dat)",options=options)
        if fileName:
            self.get_thread = Worker(index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,dyn_value,i_min,i_max,j,j_min,exceptions)
            self.get_thread.start()
            self.get_thread.i_signal.connect(self.update_bar)
            self.get_thread.error_dynfit_signal.connect(self.error_msg)

    def set_init_params(self):
        #as the name say's
        global params
        global i
        global init_values
        global bound_min
        global bound_max
        '''
        global default_values_L
        global default_values_D
        global default_boundaries_L_min
        global default_boundaries_L_max
        global default_boundaries_D_min
        global default_boundaries_D_max
        '''

        #function that reloads the parameter arrays, grabbs the values from the widgets and writes them into the arrays which are used for fitting

        if index_model == 2:
            for zahl in range(0,2):
                init_lin = self.ui.Parameter_table.cellWidget(zahl,0).value() # Inital value linear
                bounds_min_lin = self.ui.Parameter_table.cellWidget(zahl,1).value() # Boundary minimum linear
                bounds_max_lin = self.ui.Parameter_table.cellWidget(zahl,2).value() # Boundary maximum linear

                default_linear[zahl] = init_lin
                default_boundaries_linear_min[zahl] = bounds_min_lin
                default_boundaries_linear_max[zahl] = bounds_max_lin

            for zähler in range(2,fit_num*3+2):
                init = self.ui.Parameter_table.cellWidget(zähler,0).value() # Inital value
                bounds_min = self.ui.Parameter_table.cellWidget(zähler,1).value() # Boundary minimum
                bounds_max =self.ui.Parameter_table.cellWidget(zähler,2).value() # Boundary maximum

                default_values_L[zähler] = init
                default_boundaries_L_min[zähler] = bounds_min
                default_boundaries_L_max[zähler] = bounds_max

            init_values = default_linear + default_values_L
            bound_min = default_boundaries_linear_min + default_boundaries_L_min
            bound_max = default_boundaries_linear_max + default_boundaries_L_max

        else:
            for zahl in range(0,2):
                init_lin = self.ui.Parameter_table.cellWidget(zahl,0).value() # Inital value linear
                bounds_min_lin = self.ui.Parameter_table.cellWidget(zahl,1).value() # Boundary minimum linear
                bounds_max_lin = self.ui.Parameter_table.cellWidget(zahl,2).value() # Boundary maximum linear

                default_linear[zahl] = init_lin
                default_boundaries_linear_min[zahl] = bounds_min_lin
                default_boundaries_linear_max[zahl] = bounds_max_lin

            for zähler in range(2,fit_num*4+2):
                init = self.ui.Parameter_table.cellWidget(zähler,0).value() # Inital value
                bounds_min = self.ui.Parameter_table.cellWidget(zähler,1).value() # Boundary minimum
                bounds_max = self.ui.Parameter_table.cellWidget(zähler,2).value() # Boundary maximum

                default_values_D[zähler-2] = init
                default_boundaries_D_min[zähler-2] = bounds_min
                default_boundaries_D_max[zähler-2] = bounds_max

            init_values = default_linear + default_values_D
            bound_min = default_boundaries_linear_min + default_boundaries_D_min
            bound_max = default_boundaries_linear_max + default_boundaries_D_max

        self.ui.label_show_init_params.setText('Parameter set!')

        if int(self.ui.select_datanumber.text()) > len(Adata)-1:
            i = len(Adata)-1
            self.ui.select_datanumber.setText(str(i))
        else:
            i = int(self.ui.select_datanumber.text())

    def set_fit_params(self):
        #reading names is adviced
        global temp_paras
        global default_values_D
        global default_values_L
        global default_linear
        if fit_value == 'fitted':
            if index_model == 2:
                temp_paras = Fit(index_model,fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).give_params(fit_num, parameter_table_names_L_final,index_model,Adata2,Bdata2, j_min,j) #grabs the params file from different class, Lorentz
            else:
                temp_paras = Fit(index_model,fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).give_params(fit_num,parameter_table_names_D_final,index_model,Adata2,Bdata2, j_min,j) #grabs the params file from different class, Dyson
            self.refresh_inital_parameter() #Well ... guess what it does
        else:
            self.plot() #....

    def refresh_inital_parameter(self):
        # Oh hey still here?
        if index_model == 2:
            for zähler in range(0,fit_num*3+2):
                self.ui.Parameter_table.cellWidget(zähler,0).setValue(temp_paras[zähler]) # Inital value
        else:
            for zähler in range(0,fit_num*4+2):
                self.ui.Parameter_table.cellWidget(zähler,0).setValue(temp_paras[zähler]) # Inital value


if __name__=="__main__":
    #appctxt = ApplicationContext() #
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    #exit_code = appctxt.app.exec_()  #
    #sys.exit(exit_code) #
    sys.exit(app.exec_())
