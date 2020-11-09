import sys
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import math as m
import Functions
import ast
import py2mat
import time
import multiprocessing as multip

import datetime
from lmfit import Model
from lmfit import Parameters
from PyQt5.QtWidgets import QMainWindow, QApplication,QFileDialog,QDoubleSpinBox, QCheckBox, QLabel,QMessageBox
from PyQt5.QtCore import QThread, pyqtSignal, QSignalBlocker
from Fitprogramm import *
from arrays import *
from fitting import Fit
from parameter_plot import Params_Plot
from multiprocessing import Process
from ani_tools import *

# TODO: Use ExpressionModel from lmfit to generate custom function to Fit. Then from the Model parameter dict() generate the names to put into the parameterTable
# TODO: Maybe implement a performance mode using CuPy or MOT Module (CUDA and multiprocess optimized optimzation)
# TODO: Robust Fit mit Multithreading?

def define_value_opt():
    global value_opt
    value_opt = dict()
    value_opt['data'] = 'unloaded' #data_value
    value_opt['fit'] = 'unfitted' #fit_value
    value_opt['dyn'] = 'single' #dyn_value
    value_opt['dyn_fit'] = 'unfitted' #dyn_fit_value
    value_opt['parameter'] = 'unloaded' #parameter_value
    value_opt['params_copied'] = False #Value to estimate whether the dyn plot is changeable or not  
    value_opt['robust'] = False # Using Robust Fitting or not (Robust: Using Method "nelder" or "ampgo")
    value_opt['index_model_num'] = int
    value_opt['ani_pre_fit'] = True # used in get_shift

p = 0
j = 350
j_min = 0
colormap = 'magma'
tick_number = 50 #Max = 10000, bereits 1000 erstellt 200 MB file!
excep_append = []
exceptions = [] #Array to store excepted lines in the fit
fit_num = 1
Debug = 0

class Worker(QThread):
    #This class is responsible for the dynamic fit by creating a so called worker to do the job
    i_signal = pyqtSignal() # signal for the processBar
    error_dynfit_signal = pyqtSignal() # error signal
    def __init__(self,index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,value_opt,i_min,i_max,j,j_min,exceptions):
        self.index_model = index_model
        self.fit_num = fit_num
        self.Bdata = Bdata
        self.Adata = Adata
        self.Winkeldata = Winkeldata
        self.i = i
        self.fileName = fileName
        self.value_opt = value_opt
        self.i_min = i_min
        self.i_max = i_max
        self.j = j
        self.j_min = j_min
        self.exceptions = exceptions
        QThread.__init__(self)

    def fit(self,l,params_table,model,robust):
        Bdata2 = self.Bdata[l]
        Adata2 = self.Adata[l]
        if robust:
            result_dyn = model.fit(Adata2[self.j_min:self.j], params_table, B=Bdata2[self.j_min:self.j],method = 'ampgo') #lbfgsb is also a good alternative or ampgo
        else:
            result_dyn = model.fit(Adata2[self.j_min:self.j], params_table, B=Bdata2[self.j_min:self.j])
        return result_dyn

    def run(self):
        global value_opt
        global Parameter_list
        names = []
        temp_paras = []
        params_table = Fit(self.index_model, self.fit_num,Adata2,Bdata2, self.j_min,self.j,init_values,bound_min,bound_max).give_param_table(self.index_model)
        model = Fit(self.index_model, self.fit_num,Adata2,Bdata2, self.j_min,self.j,init_values,bound_min,bound_max).give_model(self.index_model)

        #clean up, just in case
        names.clear()
        temp_paras.clear()

        for name, param in result.params.items():
            names.append('{}'.format(name))  #create an array of names for the header
            temp_paras.append(float(param.value))  #create the first parameters to save into the .dat file
        temp_paras.append(float(self.Winkeldata[0]))
        print(names)
        #exceptions = self.give_exceptions()
        print(self.exceptions)
        Parameter_list = np.zeros(shape=(self.i_max-len(self.exceptions),len(temp_paras)))
        iteration = 0
        for l in range(i_min,i_max):
            if l not in self.exceptions:
                self.i_signal.emit() # update progressbar
                temp_paras.clear() # clear temp_paras for each iteration
                temp_result = self.fit(l,params_table,model,value_opt['robust']) #this is the fit
                for name, param in temp_result.params.items():
                    temp_paras.append(float(param.value))  
                    params_table[name].set(value=param.value,min=None,max=None)
                temp_paras.append(float(self.Winkeldata[l]))
                Parameter_list[iteration] = temp_paras
                iteration += 1  #cleaned iteration (iteration with exceptions removed)
        now = datetime.datetime.now()
        np.savetxt(fileName,Parameter_list,delimiter='\t',newline='\n', header='FMR-Fit\nThis data was fitted {} using: $.{} Lineshape  \nDropped Points {}     \nData is arranged as follows {}'.format(now.strftime("%d/%m/%Y, %H:%M:%S"),self.index_model,self.exceptions,names))
        value_opt['dyn_fit'] = 'fitted'

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
        self.ui.plot_parameter.clicked.connect(self.parameter_plot)
        self.ui.load_params_from_file.clicked.connect(self.load_parameter_plot)
        self.ui.parameter_data_select.valueChanged.connect(self.change_parameter_angle)
        self.ui.checkBox_dynPlot.stateChanged.connect(self.robust_fit)
        self.ui.comboBox_fit_model.activated.connect(self.make_parameter_table)
        self.ui.pushButton.clicked.connect(self.test)
        self.ui.Button_manual_save.clicked.connect(self.save_adjusted)
        self.ui.sumbit_mathematica.clicked.connect(self.mathematica_submit)
        self.ui.sumbit_mathematica_2.clicked.connect(self.python_submit)
        self.ui.shift_SpinBox.valueChanged.connect(self.get_shift)
        self.ui.shift_SpinBox_Py.valueChanged.connect(self.get_shift)
        self.ui.preview_button_Py.clicked.connect(self.make_preB_sim)
        #self.ui.LineEdit_Start_Val_Py.editingFinished.connect(self.make_preB_sim)
        self.show()

    def test(self):
        try:
            self.plot_ani_fit()
        except Exception as e:
            print(e)

    def make_preB_sim(self):
        global value_opt
        self.preB_Sim = self.python_submit(True)
        # if no pre is present, then ani_pre_fit is True
        # ani_pre_fit determines whether it is the first time fitting or just a shift
        if value_opt['ani_pre_fit']:
            value_opt['ani_pre_fit'] = False
            try:
                self.maxWinkel = max(self.preB_Sim[2])
                self.minWinkel = min(self.preB_Sim[2])
            except Exception as e:
                print('Error in main.make_preB_sim(): ',e)
                print('Try fitting the spectra first!')
        shift = self.ui.shift_SpinBox_Py.value()
        self.get_shift(shift)
        #preB_Sim: [0] = B_Sim; [1] = B_Exp; [2] = phi_RANGE_deg

    def get_shift(self,d):
        # d is returned value from spinbox
        global value_opt
        try:
            if value_opt['ani_pre_fit']:
                self.ani_fit_angle = False
                #print(d,value_opt['ani_pre_fit'])
                self.make_preB_sim()
                
                value_opt['ani_pre_fit'] = False
                self.update_canvas(self.preB_Sim,self.preB_Sim[2])
                self.changeing_phi = np.copy(self.preB_Sim[2])
                # call function to generate Sim to find good shift
            else:
                # update plot according to shift d,
                # by calculating new phi_shifted, then plot
                self.changeing_phi = np.copy(self.preB_Sim[2])
                angle = make_phirange(d, self.changeing_phi, 'deg', self.minWinkel, self.maxWinkel)
                self.update_canvas(self.preB_Sim, angle)
                self.ani_fit_angle = np.copy(angle)
        except Exception as e:
            print('Error in get_shift: ',e)

    def update_canvas(self,data:list,angle:list):
        #Data is 2D List: data[0] = B_Sim, data[1] = B_Exp
        # Mathematica Canvas
        self.ui.Ani_Const_Plot.canvas.ax.clear()
        self.ui.Ani_Const_Plot.canvas.ax.set_ylabel('Resonance Field [T]')
        self.ui.Ani_Const_Plot.canvas.ax.set_xlabel('Angle [Deg]')

        self.ui.Ani_Const_Plot.canvas.ax.scatter(angle, data[1], color='black', marker='o',
                                                 label='Experimental data')  # Plot experimental Data
        self.ui.Ani_Const_Plot.canvas.ax.plot(self.preB_Sim[2], data[0], 'r--', label='B_res Simulation')

        self.ui.Ani_Const_Plot.canvas.ax.legend()
        self.ui.Ani_Const_Plot.canvas.draw()

        # Python canvas
        self.ui.Ani_Const_Plot_Py.canvas.ax.clear()
        self.ui.Ani_Const_Plot_Py.canvas.ax.set_ylabel('Resonance Field [T]')
        self.ui.Ani_Const_Plot_Py.canvas.ax.set_xlabel('Angle [Deg]')

        self.ui.Ani_Const_Plot_Py.canvas.ax.scatter(angle, data[1], color='black', marker='o',
                                                 label='Experimental data')  # Plot experimental Data
        self.ui.Ani_Const_Plot_Py.canvas.ax.plot(self.preB_Sim[2], data[0], 'r--', label='B_res Simulation')

        self.ui.Ani_Const_Plot_Py.canvas.ax.legend()
        self.ui.Ani_Const_Plot_Py.canvas.draw()

    def robust_fit(self):
        global value_opt
        if self.ui.checkBox_dynPlot.checkState() == 2:
            value_opt['robust'] = True
        else:
            value_opt['robust'] = False

    def save_adjusted(self):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Please select the file to save to","","All Files (*);;Text Files (*.dat)",options=options)
        if fileName:
            now = datetime.now()
            np.savetxt(fileName, Para_cp, delimiter='\t', newline='\n',
                       header='FMR-Fit\nThis data was fitted {} using: $.{} Lineshape  \nDropped Points {}'.format(
                           now.strftime("%d/%m/%Y, %H:%M:%S"), index_model, exceptions))

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
        global Para_cp
        global W
        global W_real
        params = []
        try:
            #Temporarily block every signal from the Spinboxes. There is the Problem,
            # that everytime the value from these spinboxes is changed a signal is emited.
            # The bad thing now is, that there is no differentiation between user_event and machine_event,
            # therefore if the script changes the values to display the next angle,
            # the .connect() function will be called len(params) times, which has a min of 5 times with up to 42 times.
            for i in range(0,fit_num*value_opt['index_model_num']+2):
                self.ui.Parameter_Plot_Table.cellWidget(i, 0).blockSignals(True)
            W_int = self.ui.parameter_data_select.value()
            W = int(New_W[0][W_int])
            W_real = New_W[1][W_int]
            self.ui.label_manual_edit_angle.setText(str(W_real))

            #Define Big Array that includes every Parameter
            if value_opt['dyn_fit'] == 'fitted':
                Para_orig = np.array(Parameter_list)
            elif value_opt['parameter'] == 'loaded':
                Para_orig = Parameter_from_text
            else:
                print('Please select a Parameterset!!')
                self.load_parameters_to_text()
                Para_orig = Parameter_from_text
            if not value_opt['params_copied']:
                #Copy Parameterarray, in order to save the original file
                Para_cp = np.copy(Para_orig)
                value_opt['params_copied'] = True

            #Take just the Parameters representing the current angle W
            if self.ui.checkBox_change_values.isChecked() == False:
                params = Para_orig[W] #params has Entries like : [slope,offset,...,angle]
            else:
                params = Para_cp[W]
                #params = self.get_current_manual_params()
            self.update_parameter_display(params,False)
            self.plot_fitted_params(W,params)

            for i in range(0,fit_num*value_opt['index_model_num']+2):
                self.ui.Parameter_Plot_Table.cellWidget(i, 0).blockSignals(False)
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
        try:
            self.saveFileDialog()
        except Exception as e:
            print("Error in dyn_Fit: ",e)

    def set_model_type_number(self):
        #sets the integer number of index_model
        global value_opt
        global index_model
        if index_model == 2:
            value_opt['index_model_num'] = 3
        elif index_model == 3:
            value_opt['index_model_num'] = 4

    def model_type(self):
        # selects the type of function used in the fit
        global index_model
        index_model = self.ui.comboBox_fit_model.currentIndex()
        self.set_model_type_number()

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
                #table for lorentz (Fitting)
                self.ui.Parameter_table.setRowCount(fit_num*3+2)
                self.ui.Parameter_Plot_Table.setRowCount(fit_num * 3 + 2)
                self.ui.Parameter_table.setCellWidget(0, 0, QDoubleSpinBox()) # initial_values slope
                self.ui.Parameter_table.setCellWidget(0, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(0, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(0, 3, QCheckBox()) # use for fitting ?  

                self.ui.Parameter_table.setCellWidget(1, 0, QDoubleSpinBox()) # initial_values offset
                self.ui.Parameter_table.setCellWidget(1, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(1, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(1, 3, QCheckBox()) # use for fitting ?

                for i in range(2):
                    ch = QDoubleSpinBox(parent=self.ui.Parameter_Plot_Table)
                    #ch.valueChanged.connect(self.select_fitted_params) #Set signal for Spinbox in Table. Function test() is called
                    self.ui.Parameter_Plot_Table.setCellWidget(i, 0, ch)  # Value

                parameter_table_names_L_final.append('slope')
                parameter_table_names_L_final.append('offset')   

                for index in range(1,fit_num+1):
                    for list_index in parameter_Table_names_L:
                        parameter_table_names_L_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_L_final)
                self.ui.Parameter_Plot_Table.setVerticalHeaderLabels(parameter_table_names_L_final)

                for zähler in range(0,fit_num*3):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox()) # initial_values
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox()) # bound_min
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox()) # bound_max
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox()) # use for fitting ?
                #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    ch = QDoubleSpinBox(parent=self.ui.Parameter_Plot_Table)
                    #ch.valueChanged.connect(self.select_fitted_params) #Set signal for Spinbox in Table. Function test() is called
                    self.ui.Parameter_Plot_Table.setCellWidget(zähler+2, 0, ch) # Value

                self.set_default_values()

            else:
                #table for dyson (at the moment!) (Fitting)
                self.ui.Parameter_table.setRowCount(fit_num*4+2)
                self.ui.Parameter_Plot_Table.setRowCount(fit_num * 4 + 2)
                self.ui.Parameter_table.setCellWidget(0, 0, QDoubleSpinBox()) # initial_values
                self.ui.Parameter_table.setCellWidget(0, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(0, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(0, 3, QCheckBox()) # use for fitting ?  

                self.ui.Parameter_table.setCellWidget(1, 0, QDoubleSpinBox()) # initial_values offset
                self.ui.Parameter_table.setCellWidget(1, 1, QDoubleSpinBox()) # bound_min
                self.ui.Parameter_table.setCellWidget(1, 2, QDoubleSpinBox()) # bound_max
                self.ui.Parameter_table.setCellWidget(1, 3, QCheckBox()) # use for fitting ?

                for i in range(2):
                    ch = QDoubleSpinBox(parent=self.ui.Parameter_Plot_Table)
                    #ch.valueChanged.connect(self.select_fitted_params) #Set signal for Spinbox in Table. Function test() is called
                    self.ui.Parameter_Plot_Table.setCellWidget(i, 0, ch)  # Value

                parameter_table_names_D_final.append('slope')
                parameter_table_names_D_final.append('offset')  

                for index in range(1,fit_num+1):
                    for list_index in parameter_Table_names_D:
                        parameter_table_names_D_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_D_final)
                self.ui.Parameter_Plot_Table.setVerticalHeaderLabels(parameter_table_names_D_final)

                for zähler in range(0,fit_num*4):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox())

                #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    ch = QDoubleSpinBox(parent=self.ui.Parameter_Plot_Table)
                    #ch.valueChanged.connect(self.select_fitted_params)  #Set signal for Spinbox in Table. Function test() is called
                    self.ui.Parameter_Plot_Table.setCellWidget(zähler+2, 0, ch) # Value

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

                    #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setSingleStep(default_stepsize_linear[zahl])
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setValue(default_linear[zahl]) # Inital value

                    self.ui.Parameter_Plot_Table.cellWidget(zahl, 0).valueChanged.connect(self.signal_spinbox_manual_params) #different approach to connect a signal to the spinbox


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

                    #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setDecimals(4)
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setMaximum(default_maximum_bound_spinbox_L[zähler])
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setSingleStep(default_stepsize_L[zähler])
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setValue(default_values_L[zähler]) # Inital value
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox

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

                    #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setMaximum(default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setSingleStep(default_stepsize_linear[zahl])
                    self.ui.Parameter_Plot_Table.cellWidget(zahl,0).setValue(default_linear[zahl]) # Inital value
                    self.ui.Parameter_Plot_Table.cellWidget(zahl, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox

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

                    #--------------------------------------Fitting Table ends. Below is Table for parameter Plot-----------
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setDecimals(4)
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setMaximum(default_maximum_bound_spinbox_D[zähler])
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setSingleStep(default_stepsize_D[zähler])
                    self.ui.Parameter_Plot_Table.cellWidget(zähler+2,0).setValue(default_values_D[zähler])
                    self.ui.Parameter_Plot_Table.cellWidget(zähler + 2, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox
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
        except Exception as e:
            print('Error in parameter_plot:',e)

    def get_fit_options_from_file(self,fname):
        global index_model
        global fit_num
        global value_opt

        f = open(fname)
        Title = f.readline()

        header_parts = []

        if Title.split('# ')[-1] == 'FMR-Fit\n':
            header_parts.append(Title)
            for i in range(3):
                header_parts.append(f.readline())

            Lineshape = int(header_parts[1].split('$.')[1].split(' Lineshape')[0])  # Extract lineshape out of params file

            try:
                # Takes Element 3 of header_parts and converts it into an usable array
                Dropped_points = np.asarray(ast.literal_eval(header_parts[2].split('[')[1].split(']')[0]))
            except Exception as e:
                print('Error trying to get Dropped_Points, assuming Dropedpoints as empty:', e)
                Droppen_points = None

            # Takes Element 4 of header_parts and converts it into an usable array: ['A','dB','R',....]
            Params_name = np.asarray(ast.literal_eval(header_parts[3].split('[')[1].split(']')[0]))

            index_model = Lineshape
            print(value_opt['index_model_num'])
            self.set_model_type_number()
            print(value_opt['index_model_num'])
            fit_num = int((len(Params_name)-2)/value_opt['index_model_num'])
        else:
            print('File was not created by this script!')

    def load_parameter_plot(self):
        print('Loading from file!')
        try:
            self.load_parameters_to_text()
            Params_Plot(Parameter_from_text,'load',index_model)
        except Exception as e:
            print('Error in load_parameters_to_text',e)

    def load_parameters_to_text(self):
        global Parameter_from_text
        global value_opt
        #Todo: Remove as much global variables as possible!
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            self.get_fit_options_from_file(params_fname[0])
            Parameter_from_text = np.loadtxt(params_fname[0],dtype='float',skiprows=0)
        value_opt['parameter'] = 'loaded'

    def update_parameter_display(self,params,spinbox_bool):
        if not spinbox_bool:
            for zaehler in range(0, fit_num * value_opt['index_model_num']+2):
                self.ui.Parameter_Plot_Table.cellWidget(zaehler, 0).setValue(params[zaehler])

    def update_current_params(self):
        #-------------------------------------Todo: Find solution without global parameter
        global Para_cp
        params = []
        multiplikator = value_opt['index_model_num']
        #Params array with slope and offset inside!! So : [slope,offset,....]
        for zaehler in range(0, fit_num * multiplikator+2):
            params.append(self.ui.Parameter_Plot_Table.cellWidget(zaehler, 0).value())
        params.append(W_real)
        self.plot_fitted_params(W,params)
        Para_cp[W] = params
        #return params

    def signal_spinbox_manual_params(self):
        global Debug
        #spinbox_bool = True
        # used to estimate if the next function call comes from changed values in spinbox or changed value of angle slider
        #self."insert func name here"(spinbox_bool)
        #self.update_parameter_display(Para_cp[W],True)
        #self.update_current_params()
        if self.ui.checkBox_change_values.isChecked() == True:
            Debug += 1
            #if Debug == 5:
            self.update_current_params()
            Debug = 0

    def plot_fitted_params(self,W,params):
        global Debug
        Debug = 0
        x = Bdata[W][j_min:j]
        y = Adata[W][j_min:j]
        slope = params[0]
        offset = params[1]
        params = params[2:len(params)]
        self.ui.parameter_plot_widget.canvas.ax.clear()
        self.ui.parameter_plot_widget.canvas.ax.set_xlabel('Magnetic Field [T]')
        self.ui.parameter_plot_widget.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
        self.ui.parameter_plot_widget.canvas.ax.scatter(x, y, color='black', marker='o', label='Experimental data') #Plot experimental Data

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
        self.model_type()
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
        global tick_number
        global colormap
        '''#colour plotting routine using mayavi (OpenGL renderer)
        try:
            from mayavi import mlab
            tick_number = self.ui.colour_tick_edit.text() #not needed for mayavi
            colormap = self.ui.colour_map_edit.text()   #can also be changed in mayavi UI
            mlab.mesh(X,Y,Z, colormap=colormap,extent=[0.4,1.2,90,120,-300,300])
            mlab.xlabel('Magnetic Field')
            mlab.ylabel('Angle')
            mlab.colorbar(title='Amplitude [Arb. Units]')
            mlab.show()
        except:
            print('Mayavi not isntalled!')'''

        font = {'family': 'DejaVu Sans',
                'weight': 'bold',
                'size': '25'
                }
        mpl.rc('font', **font)

        tick_number = self.ui.colour_tick_edit.text()
        colormap = self.ui.colour_map_edit.text()
        fig, ax = plt.subplots()
        #contour = ax.contourf(X,Y,Z,30,cmap='magma')
        contour = ax.contourf(X, Y, Z, int(tick_number), cmap=str(colormap))
        cbar = fig.colorbar(contour)
        cbar.set_label('Amplitude [Arb. U.]')
        plt.xlabel('Magnetic Field [mT]')
        plt.ylabel('Angle [deg]')

        plt.show()


    def Exit(self):
        sys.exit()  #Obviously not the start

    def saveFileDialog(self):
        #started when Dyn Fit button pressed
        #gets the filename for saving
        #then starts the worker responsible for dyn fit
        global p
        global New_W
        global fileName

        p=0
        exceptions = self.Exceptions()
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
            self.get_thread = Worker(index_model,fit_num,Bdata,Adata,Winkeldata,i,fileName,value_opt['dyn'],i_min,i_max,j,j_min,exceptions)
            self.get_thread.start()
            self.get_thread.i_signal.connect(self.update_bar)
            self.get_thread.error_dynfit_signal.connect(self.error_msg)

    def python_submit(self,*args):
        # Python AniFit routine
        #try:
        # FreeEnergy from GUI has to be converted to Python synatx
        F = self.ui.LineEdit_free_E_den_Py.text()
        F = self.convert_freeE_to_py(F) # This F is now in Python Syntax!

        #Now create dicts for fit_Params and fixed_params
        fit_params = self.ui.LineEdit_fitted_params_Py.text().replace(' ','').split(',')
        fit_params_value = self.ui.LineEdit_Start_Val_Py.text().replace(' ','').split(',')
        fit_params_value = [float(i) for i in fit_params_value]
        ani_fit_params = dict(zip(fit_params,fit_params_value))     #dict of fit params

        fixed_params = self.ui.LineEdit_fixed_Param_Py.text().replace(' ','').split(',')
        fixed_params_value = self.ui.LineEdit_fixed_values_Py.text().replace(' ','').replace('^','**').replace('Pi','m.pi').split(',')
        fixed_params_value = [eval(i) for i in fixed_params_value]
        ani_fixed_params = dict(zip(fixed_params,fixed_params_value))   #dict of fixed params

        anglestep = m.pi / self.ui.spinBox_Anglestep_Py.value()
        shift = float(self.ui.shift_SpinBox_Py.value())
        #Then call init_load from ani_tool.py

        try:
            if not args[0]:
                init_load(fileName, F, ani_fit_params, ani_fixed_params, shift, anglestep, True, True)
            else:
                result = init_load(fileName, F, ani_fit_params, ani_fixed_params, shift, anglestep, False, False)
                return result
        except Exception as e:
            print('Error in main.python_submit(): ',e)
            print('Try fitting the spectra first!')
        #except Exception as e:
        #    if len(fit_params) != len(fit_params_value):
        #        print('Length of fit parameters and start values is unequal!')
        #    elif len(fixed_params) != len(fixed_params_value):
        #        print('Length of fixed parameters is unequal!')
        #    print('Error in python_submit: ',e)



    def convert_freeE_to_py(self,FreeE):

        #Convert FreeE to python syntax
        #FreeE is string of FreeE from GUI

        #I'm not satisfied with this approach. Open for suggestions, on how to improve this!

        FreeE = FreeE.replace('[', '(')
        FreeE = FreeE.replace(']', ')')
        FreeE = FreeE.replace(') ', ')')
        FreeE = FreeE.replace(' (', '(')
        FreeE = FreeE.replace('^', '**')
        FreeE = FreeE.replace('S', 's')
        FreeE = FreeE.replace(' s', 's')
        FreeE = FreeE.replace('C', 'c')
        FreeE = FreeE.replace(' c', 'c')
        FreeE = FreeE.replace(' - ', '-')
        FreeE = FreeE.replace('- ', '-')
        FreeE = FreeE.replace(' -', '-')
        FreeE = FreeE.replace(' + ', '+')
        FreeE = FreeE.replace('+ ', '+')
        FreeE = FreeE.replace(' +', '+')
        FreeE = FreeE.replace('Theta', 'theta')
        FreeE = FreeE.replace('Phi', 'phi')
        FreeE = FreeE.replace('phib', 'phiB')
        FreeE = FreeE.replace('thetab', 'thetaB')
        FreeE = FreeE.replace('phiu', 'phiU')
        FreeE = sympify(FreeE)
        return FreeE

    def mathematica_submit(self):
        #(fileName,F,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,outputPath,BresColumn,WinkelColumn,Shift)
        F = self.ui.LineEdit_free_E_den.text()
        Params = "fitParameters = " + self.create_mathematica_array(self.ui.LineEdit_fitted_params.text()) #'fitParameters = {K2p, K2s, K4p, phiu}'
        StartVal = "startvalues = " + self.create_mathematica_array(self.ui.LineEdit_Start_Val.text()) #'startvalues = {863.25, 261345, 13720.6, 5.0756}'
        Ranges = self.create_mathematica_array(self.ui.LineEdit_Range.text())  #'{0.5, 0.5, 0.5, 0.5}'
        Steps = self.create_mathematica_array(self.ui.LineEdit_Steps.text())   #'{0.5, 0.5, 0.5, 0.5}'
        fixedParams = self.create_mathematica_array(self.ui.LineEdit_fixed_Param.text()) #'{omega, g, M, K4s}'
        fixedValues = self.create_mathematica_array(self.ui.LineEdit_fixed_values.text()) #'{2 Pi*9.8782*10^9, 2.05, 1.53*10^6, 0}'
        anglestep = 'Pi/'+str(self.ui.spinBox_Anglestep.value())
        iterations = self.ui.spinBox_Iterations.value()
        BresColumn = 4
        WinkelColumn = 6
        Shift = self.ui.shift_SpinBox.value()
        #print(Params,StartVal)
        choice = QMessageBox.question(self, 'Sumbitting!',"This fitting can take a while!\nThe GUI will be unresponsive after submission, are you sure to continue?")
        try:
            if choice == QMessageBox.Yes:
                #p = Process(target=py2mat.submit(fileName,F,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,BresColumn,WinkelColumn,Shift))
                self.thread = py2mat.Py2Mat(fileName,F,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,BresColumn,WinkelColumn,Shift)
                self.thread.save_path_signal.connect(self.get_output_path)
                self.thread.start()
                #py2mat.submit(fileName,F,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,BresColumn,WinkelColumn,Shift)
                '''while not self.thread.isFinished():
                    time.sleep(1)
                    print(save_path)'''
            else:
                print('Submission aborted')
        except Exception as e:
            print(e)

    def get_output_path(self,arg):
        #This function is called from the Py2Mat Worker after the ANi Fit finished. It will emit the Signal "arg",
        # which is the path to the output folder
        global save_path
        save_path = arg
        self.plot_ani_fit()

    def shift_points(self,array,shift):
        for i,l in enumerate(array):
            array[i] += shift
            if array[i] > 355.0:
                array[i] -= 360
        return array

    def plot_ani_fit(self):
        try:
            ani_out = np.loadtxt(save_path + 'outputFunktion.dat',dtype='float')
            Shift = 37  #Temporary Variable

            Angle = np.multiply(ani_out[:,0], 180 / m.pi)
            B_res = ani_out[:,1]

            if index_model == 2:
                Angle_exp = self.shift_points(Parameter_list[:,5],Shift)
                B_res_exp = Parameter_list[:,3]
            else:
                Angle_exp = self.shift_points(Parameter_list[:,6],Shift)
                B_res_exp = Parameter_list[:,4]

            self.ui.Ani_Const_Plot.canvas.ax.clear()
            self.ui.Ani_Const_Plot.canvas.ax.set_ylabel('Resonance Field [T]')
            self.ui.Ani_Const_Plot.canvas.ax.set_xlabel('Angle [Deg]')

            self.ui.Ani_Const_Plot.canvas.ax.scatter(Angle_exp, B_res_exp, color='black', marker='o', label='Experimental data') #Plot experimental Data
            self.ui.Ani_Const_Plot.canvas.ax.plot(Angle, B_res, 'r--', label='Fitted Anisotropy')

            self.ui.Ani_Const_Plot.canvas.ax.legend()
            self.ui.Ani_Const_Plot.canvas.draw()
        except Exception as e:
            print(e)


    def create_mathematica_array(self,text):
        #text is the text given by the lineEdits
        mat_array = "{" + text + "}"
        return mat_array

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
    mutlip.set_start_method("spawn")
    define_value_opt()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MyForm()
    w.show()
    #exit_code = appctxt.app.exec_()  #
    #sys.exit(exit_code) #
    sys.exit(app.exec_())
