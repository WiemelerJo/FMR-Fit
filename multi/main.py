import sys
import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import math as m
import threading
import Functions
from datetime import datetime
from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from PyQt5.QtWidgets import QMainWindow, QApplication, QAction,QFileDialog,QProgressBar,QDoubleSpinBox, QCheckBox
from PyQt5.QtCore import QThread, pyqtSignal
from Fitprogramm import *
from arrays import *
#from fbs_runtime.application_context.PyQt5 import ApplicationContext
#from matplotlib.colors import BoundaryNorm
#from matplotlib.ticker import MaxNLocator
from fitting import Fit
#from fitting import Worker
from parameter_plot import Params_Plot


def define_value_opt():
    global value_opt
    value_opt = dict()
    value_opt['data'] = 'unloaded' #data_value
    value_opt['fit'] = 'unfitted' #fit_value
    value_opt['dyn'] = 'single' #dyn_value
    value_opt['dyn_fit'] = 'unfitted' #dyn_fit_value
    value_opt['parameter'] = 'unloaded' #parameter_value
    value_opt['params_copied'] = False #Value to estimate whether the dyn plot is changeable or not  
    value_opt['robust'] = False # Using Robust Fitting or not (Robust: Using Method "nelder")

p = 0
j = 350
j_min = 0
colormap = 'inferno'
tick_number = 256 #Max = 10000, bereits 1000 erstellt 200 MB file!
excep_append = []
exceptions = [] #Array to store excepted lines in the fit
fit_num = 1

class Worker(QThread):
    #This class is responsible for the dynamic fit by creating a so called worker to do the job
    i_signal = pyqtSignal() # signal for the processBar
    error_dynfit_signal = pyqtSignal() # error signal
    def __init__(self,index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,value_opt,i_min,i_max,j,j_min,exceptions):
        QThread.__init__(self)

    def fit(self,l,params_table,model,robust):
        Bdata2 = Bdata[l]
        Adata2 = Adata[l]
        if robust:
            result_dyn = model.fit(Adata2[j_min:j], params_table, B=Bdata2[j_min:j],method = 'ampgo') #lbfgsb is also a good alternative or ampgo
        else:
            result_dyn = model.fit(Adata2[j_min:j], params_table, B=Bdata2[j_min:j])
        return result_dyn

    def run(self):
        global value_opt
        global Parameter_list
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
        temp_paras.append(float(Winkeldata[0]))
        print(names)
        #exceptions = self.give_exceptions()
        print(exceptions)
        Parameter_list = np.zeros(shape=(i_max,len(temp_paras)))
        for l in range(i_min,i_max):
            if l not in exceptions:
                self.i_signal.emit() # update progressbar
                temp_paras.clear() # clear temp_paras for each iteration
                temp_result = self.fit(l,params_table,model,value_opt['robust']) #this is the fit
                for name, param in temp_result.params.items():
                    temp_paras.append(float(param.value))  
                    params_table[name].set(value=param.value,min=None,max=None)
                temp_paras.append(float(Winkeldata[l]))
                Parameter_list[l] = temp_paras
        now = datetime.now()
        np.savetxt(fileName,Parameter_list,delimiter='\t',newline='\n', header='FMR-Fit\nThis data was fitted {} using: $.{} Lineshape  \nDropped Points {}     \nData is arranged as follows {}'.format(now.strftime("%d/%m/%Y, %H:%M:%S"),index_model,exceptions,names))
        value_opt['dyn_fit'] = 'fitted'

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

        self.ui.Button_dropped_points.clicked.connect(self.button_dropped_points)

        #self.ui.Dropped_points_edit.textEdited.connect(self.reset_drop)


        self.ui.plot_parameter.clicked.connect(self.parameter_plot)
        self.ui.load_params_from_file.clicked.connect(self.load_parameter_plot)
        self.ui.parameter_data_select.valueChanged.connect(self.change_parameter_angle)
        self.ui.checkBox_dynPlot.stateChanged.connect(self.robust_fit)
        self.ui.comboBox_fit_model.activated.connect(self.make_parameter_table)
        self.show()

    def test(self):
        print('Ja moin')

    def robust_fit(self):
        global value_opt
        if self.ui.checkBox_dynPlot.checkState() == 2:
            value_opt['robust'] = True

        else:
            value_opt['robust'] = False

    '''def reset_drop(self):
                    #empties the exception array
                    if self.ui.Dropped_points_edit.text() == '':
                        self.ui.Dropped_points_edit.setText('')
                    else:
                        x = str(self.ui.Dropped_points_edit.text())
                        excep_append = list(map(int, x.split(',')))
                        print(excep_append)'''
   
    def button_dropped_points(self):
        #as something is added to the exceptions, this is beeing called
        excep_append = self.Exceptions()  #gets text of editLine
        excep_append.append(str(self.ui.Scroll_Bar_dropped_points.value())) #get value of scrollbar
        text = ','.join(map(str, excep_append))
        self.ui.Dropped_points_edit.setText(str(text))
        self.Exceptions()

    def Exceptions(self):
        # takes the string input from the exception text field inside the GUI and seperates the entries into a list
        x = str(self.ui.Dropped_points_edit.text()) #saves everything that is in the editLine as x
        if x == '':
            exceptions = [] #if y is empty excep should also be empty
        else:
            exceptions = list(map(int, x.split(','))) #else set excep as everything that is standing in editLine
        return exceptions

    def change_parameter_angle(self):
        '''W = self.ui.parameter_data_select.value()
                                if W not in exceptions:
                                    self.select_fitted_params(W)'''
        try:
            W_int = self.ui.parameter_data_select.value()
            W = int(New_W[0][W_int])
            W_real = New_W[1][W_int]
            self.ui.label_manual_edit_angle.setText(str(W_real))   
            self.select_fitted_params(W)
        except Exception as e:
            print("Error in change_parameter_angle",e)


    def set_dataslider(self):
        self.ui.select_datanumber.setValue(int(self.ui.Scroll_Bar_dropped_points.value()))

    def set_datanumber(self):
        #Defines the dataset looked at
        global i
        global Bdata2
        global Adata2
        global value_opt
        if value_opt['data'] == 'loaded':

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
        global value_opt
        value_opt['dyn'] = 'dynamic'
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
        except Exception as e:
                print('Error in make_parameter_table',e)
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

    def set_default_values(self):
        #sets default values into the spinbox according to arrays defined in the beginning
        #it basicly is just a for loop in order to cath every row of the table according to the numbver of lines
        try:
            if index_model == 2:
                for zahl in range(0,2):
                    self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-1000)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setMaximum(default_maximum_bound_spinbox_linear[zahl])

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

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setMaximum(default_maximum_bound_spinbox_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setMaximum(default_maximum_bound_spinbox_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setMaximum(default_maximum_bound_spinbox_L[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(default_values_L[zähler]) # Inital value
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(default_boundaries_L_min[zähler]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(default_boundaries_L_max[zähler]) # Boundary maximum

            else:
                for zahl in range(0,2):
                    self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-1000)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setMaximum(default_maximum_bound_spinbox_linear[zahl])

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

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setMaximum(default_maximum_bound_spinbox_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setMaximum(default_maximum_bound_spinbox_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setMaximum(default_maximum_bound_spinbox_D[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(default_values_D[zähler]) # Inital value
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(default_boundaries_D_min[zähler]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(default_boundaries_D_max[zähler]) # Boundary maximum
        except Exception as e:
            print("Error in set_default_values:",e)

    def openFileDialog(self):
        global value_opt
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
            Winkeldata_raw = np.split(np.array(D[:,2]),chunksize)   #List of sublists, angle data
            Adata = np.split(np.array(D[:,3]),chunksize)    #List of sublists, amplitude data
            Winkeldata = []
            for i in range(chunksize):
                Winkeldata.append(Winkeldata_raw[i][0])
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
            value_opt['data'] = 'loaded'
            self.set_datanumber()
    
    def data_range(self):
        global j
        global j_min
        global Bdata2
        global Adata2
        if value_opt['data'] == 'loaded':
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
        global value_opt
        try:
            if value_opt['dyn_fit'] == 'fitted':
                Params_Plot(Parameter_list,'eigen',index_model)
            elif value_opt['parameter'] == 'loaded':
                print('Unfitted parameters. Using from file instead!')
                Params_Plot(Parameter_from_text,'eigen',index_model)
            else:
                self.load_parameters_to_text()
                Params_Plot(Parameter_from_text,'eigen',index_model)
        except NameError:
            print('Please select the Lineshape first!!')

    def load_parameter_plot(self):
        print('Loading from file!')
        try:
            params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
            if params_fname[0]:
                Params_Plot(params_fname,'load',index_model)
        except NameError:
            print('Please select the Lineshape first!!')            

    def load_parameters_to_text(self):
        global Parameter_from_text
        global value_opt
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            Parameter_from_text = np.loadtxt(params_fname[0],dtype='float',skiprows=1)
        value_opt['parameter'] = 'loaded'


    def select_fitted_params(self,W):
        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------
        global Parameter_from_text
        global value_opt
        global Para_cp

        params = []
        if self.ui.checkBox_change_values.isChecked()==True: # In order to change parameters manual and save them afterwards

            if value_opt['data'] == 'loaded':
                #Setup basic Data xy for plotting
                x_temp = Bdata[W]
                y_temp = Adata[W]
                x = x_temp[j_min:j]
                y = y_temp[j_min:j]

                #Setup Parameters
                if value_opt['dyn_fit'] == 'fitted':
                    Para_orig = np.array(Parameter_list)
                elif  value_opt['parameter'] == 'loaded':
                    Para_orig = Parameter_from_text
                else:
                    print('Please select a Parameterset!!')
                    self.load_parameters_to_text()
                    Para_orig = Parameter_from_text

                #Copy Parameters in order to change them manualy, but without corrupting the orignal data
                if not value_opt['params_copied']:
                    Para_cp = np.copy(Para_orig)
                    value_opt['params_copied'] = True
                    slope = Para_cp[:,0][W]
                    offset = Para_cp[:,1][W]

                    Para_dim = Para_cp.shape[1]
                    for i in range(2,Para_dim): #loop that creates the params array for an angle W
                        params.append(Para_cp[:,i][W])

                    #------------------------------------Hier noch dringend die Labels der SpinBoxen ändern!-----------------------------------------------------------------------------------------------
                

                #Parameter_Plot_Table


                #self.ui.param_change_alpha.setValue(alpha_list[W])
                #self.ui.param_change_dB.setValue(db_list[W])            
                #self.ui.param_change_R.setValue(R_list[W])
                #self.ui.param_change_A.setValue(A_list[W])
                #self.ui.param_change_slope.setValue(slope_list[W])
                #self.ui.param_change_offset.setValue(offset_list[W])

                else:
                    #lesen der Variablen aus den Labels der SpinBoxen und dann formatieren zu params und plotten
                    print('fw')
                    


                self.plot_fitted_params(x,y,slope,offset,params) # params needs to be an array of fixed values for parameters (alpha1_value, db1_value , ...)



            else:
                print('Please select a Dataset first!')
                self.openFileDialog()            


        else:
            if value_opt['data'] == 'loaded':
                try:
                    x_temp = Bdata[W]
                    y_temp = Adata[W]
                    x = x_temp[j_min:j]
                    y = y_temp[j_min:j]

                    if value_opt['dyn_fit'] == 'fitted':
                        Para = np.array(Parameter_list)
                    elif  value_opt['parameter'] == 'loaded':
                        Para = Parameter_from_text
                    else:
                        print('Please select a Parameterset!!')
                        self.load_parameters_to_text()
                        Para = Parameter_from_text

                    slope = Para[:,0][W]
                    offset = Para[:,1][W]

                    Para_dim = Para.shape[1]
                    for i in range(2,Para_dim): #loop that creates the params array for an angle W
                        params.append(Para[:,i][W])

                    self.plot_fitted_params(x,y,slope,offset,params,) # params needs to be an array of fixed values for parameters (alpha1_value, db1_value , ...)

                    '''self.ui.param_change_alpha.setValue(alpha_list[W])
                                                        self.ui.param_change_dB.setValue(db_list[W])            
                                                        self.ui.param_change_R.setValue(R_list[W])
                                                        self.ui.param_change_A.setValue(A_list[W])
                                                        self.ui.param_change_slope.setValue(slope_list[W])
                                                        self.ui.param_change_offset.setValue(offset_list[W])'''
                except Exception as e:
                    print("Error in Func: select_fitted_params()",e) 
            else:
                print('Please select a Dataset first!')
                self.openFileDialog()
            

    def changeing_parameters(self,w):

        #----------------------------------------has to be reworked as of implementing multiple lines---------------------------------------------------

        global Parameter_from_text
        global value_opt

        if self.ui.checkBox_change_values.isChecked()==True:
            try:
                x_temp = Bdata[W]
                y_temp = Adata[W]
                x = x_temp[j_min:j]
                y = y_temp[j_min:j]

                if value_opt['dyn_fit'] == 'fitted':
                    Para = np.array(Parameter_list)
                elif  value_opt['parameter'] == 'loaded':
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

                self.plot_fitted_params(x,y,slope,offset,params) # params needs to be an array of fixed values for parameters (alpha1_value, db1_value , ...)     

            except Exception as e:
                print('Failed with Error',e)
                print('Please try to load a Dataset and a Parameterset!')

    def plot_fitted_params(self,x,y,slope, offset,params):
        self.ui.parameter_plot_widget.canvas.ax.clear()
        self.ui.parameter_plot_widget.canvas.ax.set_xlabel('Magnetic Field [T]')
        self.ui.parameter_plot_widget.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
        self.ui.parameter_plot_widget.canvas.ax.plot(x, y, color='black', marker='o', label='Experimental data') #Plot experimental Data
        
        '''if value_opt['dyn_fit'] == 'fitted':
                                    self.ui.parameter_plot_widget.canvas.ax.plot(x, result.best_fit, label='Best Fit') #Fitted Function'''
        try:
            plot_func = Functions.functions_value(x, slope, offset, index_model, params,fit_num) # params needs to be an array of fixed values for parameters (alpha1_value, db1_value , ...)
            self.ui.parameter_plot_widget.canvas.ax.plot(x,plot_func,'r--',label='Best Fit') # Best Fit

            for i_num in range(1,fit_num + 1):
                plt_single_func = Functions.single_func(x, slope, offset, index_model, params, i_num)
                self.ui.parameter_plot_widget.canvas.ax.plot(x,plt_single_func,label='Function'+str(i_num)) #Single lines
            
            self.ui.parameter_plot_widget.canvas.ax.legend()
            self.ui.parameter_plot_widget.canvas.draw()         
        except Exception as e:
            print(e)
            print('Please select a Model first!')

    def plot(self):
        global value_opt
        global result
        value_opt['dyn'] = 'single'
        self.model_type(False)
        if value_opt['data'] == 'loaded':
            try:
                self.set_init_params()  #gets the initial parameters from the GUI
                result = Fit(index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max).fit(index_model,Adata2,Bdata2, j_min,j) #fit and save it as result
                value_opt['fit'] = 'fitted'
                self.ui.progressBar.setMaximum(i_max-len(exceptions))
                #self.ui.label_params_output.setText(self.fit(i).fit_report())
                plt.ion()   #enable user actions in plot
                plt.cla()   #plot clean up
                plt.clf()   #plot clean up
                plt.plot(Bdata2[j_min:j],Adata2[j_min:j], '-b', label='Experimental data')  #plot original data
                self.evaluate_min_max(Bdata2[j_min:j],Adata2[j_min:j])
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
            except Exception as e:
                print(e)
        else:
            self.openFileDialog()
        if index_model == 0:
            self.ui.label_params_output.setText('Please select a fitting model first')

    def evaluate_min_max(self,Bdata,Adata):
        ind_min = Bdata[int(np.unravel_index(np.argmin(Adata, axis=None), Adata.shape)[0])]
        ind_max = Bdata[int(np.unravel_index(np.argmax(Adata, axis=None), Adata.shape)[0])]
        Bres = ind_max+(ind_min-ind_max)/2 #Asuming a symmetric line
        self.ui.label_test.setText('Bres from min/max: '+str(Bres))

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
        global  New_W
        global fileName
        global exceptions
        p=0

        New_W = [[],[]] #[0] = position, [1] = angle
        for pos,winkel in enumerate(Winkeldata):
            if pos not in exceptions:
                New_W[0].append(pos)
                New_W[1].append(winkel)
        New_W = np.asarray(New_W)
        #self.ui.parameter_data_select.setRange(np.amin(New_W[1]),np.amax(New_W[1]))

        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Please select the file to save to","","All Files (*);;Text Files (*.dat)",options=options)
        if fileName:
            exceptions = self.Exceptions()
            self.get_thread = Worker(index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,value_opt['dyn'],i_min,i_max,j,j_min,exceptions)
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

                default_values_L[zähler-2] = init
                default_boundaries_L_min[zähler-2] = bounds_min
                default_boundaries_L_max[zähler-2] = bounds_max

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
        if value_opt['fit'] == 'fitted':
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
    define_value_opt()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MyForm()
    w.show()
    #exit_code = appctxt.app.exec_()  #
    #sys.exit(exit_code) #
    sys.exit(app.exec_())
