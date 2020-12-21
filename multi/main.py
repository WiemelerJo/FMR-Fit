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
import pyqtgraph as pg
import datetime

from lmfit import Model
from lmfit import Parameters
from matplotlib.colors import ListedColormap
from PyQt5.QtWidgets import QMainWindow, QApplication,QFileDialog,QDoubleSpinBox, QCheckBox, QLabel, QMessageBox, QShortcut
from PyQt5.QtCore import QThread, pyqtSignal, QSignalBlocker
from PyQt5.Qt import Qt
from PyQt5.QtGui import QKeySequence
from Fitprogramm import *
#from arrays import *
from fitting import Fit
from parameter_plot import Params_Plot
from multiprocessing import Process
from ani_tools import *
from CustomWidgets import Popup_View, Fit_Log
from func_gen import Gen_Lorentz, Gen_Dyson
from array_gen import *


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
#fit_num = 1
Debug = 0

class Worker(QThread):
    #This class is responsible for the dynamic fit by creating a so called worker to do the job
    i_signal = pyqtSignal() # signal for the processBar
    error_dynfit_signal = pyqtSignal() # error signal

    job_done = pyqtSignal(list)

    def __init__(self,Bdata,Adata,i,value_opt,i_min,i_max,j,j_min,dyn_params_table):
        self.Bdata = Bdata
        self.Adata = Adata
        self.i = i
        self.value_opt = value_opt
        self.i_min = i_min
        self.i_max = i_max
        self.j = j
        self.j_min = j_min
        self.dyn_params_table = dyn_params_table
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

        for index in self.dyn_params_table[1]:  #Find first element that is not None or False
            if index != None and index != False:
                spectrum = index
                break

        clean_init = []
        clean_min = []
        clean_max = []
        for i in range(2,len(spectrum)):
            clean_init.append(spectrum[i].value)
            clean_min.append(spectrum[i].min)
            clean_max.append(spectrum[i].max)

        Fit_obj = Fit(spectrum[0], spectrum[1],Adata2,Bdata2, self.j_min,self.j,clean_init,clean_min,clean_max)
        params_table = Fit_obj.give_param_table(spectrum[0])
        model = Fit_obj.give_model()

        for l in range(i_min,i_max):
            temp_spectra = self.dyn_params_table[1][l]
            if temp_spectra == None:    # Only use spaces that are equal to None
                self.i_signal.emit()  # update progressbar
                temp_result = self.fit(l, params_table, model, value_opt['robust']) # Fit
                temp_paras.clear()  # clear temp_paras for each iteration
                for name, param in temp_result.params.items():
                    temp_paras.append(param)
                    params_table[name].set(value=param.value, min=None, max=None)

                self.dyn_params_table[1][l] = [spectrum[0],spectrum[1]] # Set index_model and fit_num
                for i in temp_paras:
                    self.dyn_params_table[1][l].append(i)  # Set parameters

        self.job_done.emit(self.dyn_params_table) # emit finished self.dyn_param_table to main Class

class MyForm(QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Setup every Button, Slider, ....
        self.ui.actionOpen.triggered.connect(self.openFileDialog)
        self.ui.actionSave.triggered.connect(self.saveFileDialog)
        self.ui.actionExit.triggered.connect(self.Exit)
        self.ui.Button_Plot.clicked.connect(self.plot)
        self.ui.comboBox_fit_model.currentIndexChanged.connect(self.model_type)
        self.ui.Button_dyn_fit.clicked.connect(self.start_worker)
        self.ui.select_datanumber.valueChanged.connect(self.set_datanumber)
        self.ui.Scroll_Bar_dropped_points.valueChanged.connect(self.set_dataslider)
        self.ui.Dropped_points_edit.editingFinished.connect(self.Exceptions)
        self.ui.Button_dropped_points.clicked.connect(self.button_dropped_points)
        self.ui.plot_parameter.clicked.connect(self.parameter_plot)
        self.ui.load_params_from_file.clicked.connect(self.load_parameter_plot)
        self.ui.checkBox_dynPlot.stateChanged.connect(self.robust_fit)
        self.ui.comboBox_fit_model.activated.connect(self.make_parameter_table)
        self.ui.pushButton.clicked.connect(self.test)
        self.ui.Button_angu_view.clicked.connect(self.doit)
        self.ui.sumbit_mathematica.clicked.connect(self.mathematica_submit)
        self.ui.sumbit_mathematica_2.clicked.connect(self.python_submit)
        self.ui.shift_SpinBox.valueChanged.connect(self.get_shift)
        self.ui.shift_SpinBox_Py.valueChanged.connect(self.get_shift)
        self.ui.preview_button_Py.clicked.connect(self.make_preB_sim)
        self.ui.spinBox_fit_num.valueChanged.connect(self.select_fit_number)
        self.ui.checkBox_fit_log.stateChanged.connect(self.fit_log)
        self.ui.Button_reset_fit.clicked.connect(self.reset_fit)

        # Set up every Key press Event

        self.shortcut_next = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_F2), self) # Goto next Spectra
        self.shortcut_next.activated.connect(self.keyboard_next_spectra)
        self.shortcut_next = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_Right), self) # Goto next spectra
        self.shortcut_next.activated.connect(self.keyboard_next_spectra)

        self.shortcut_prev = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_F3), self) # Goto previous Spectra
        self.shortcut_prev.activated.connect(self.keyboard_prev_spectra)
        self.shortcut_prev = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_Left), self)   # Goto previous Spectra
        self.shortcut_prev.activated.connect(self.keyboard_prev_spectra)

        self.shortcut_drop = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_Down), self)   # Add Spectra to dropped Points
        self.shortcut_drop.activated.connect(self.button_dropped_points)

        self.shortcut_fit = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_Up), self)  # Fit Spectra
        self.shortcut_fit.activated.connect(self.plot)

        self.fit_num = 1
        self.increment = 0.001
        self.sanity = True

        self.show()

    def test(self):
        print("Debug Funktion")
        #print(self.dyn_params_table[1])
        #print(len(self.dyn_params_table[1]))
        self.plot_params_to_plot_tab()

    def doit(self):
        self.w = Popup_View(Z,self.H_range,self.WinkelMax)
        self.w.setWindowTitle("Colourplot")
        self.w.setGeometry(0, 0, 1280, 720)
        self.w.show()

    def wheelEvent(self, event):
        modifiers = QApplication.keyboardModifiers()
        if modifiers == Qt.ControlModifier:
            angle_delta = event.angleDelta().y()
            # print('Wheel', angle_delta)
            if angle_delta > 0:
                # Positiver Step
                step = 1
                self.increment *= 10
            else:
                # Negativer Step
                step = -1
                self.increment /= 10

            if self.increment > 1000.0:
                self.increment = 1000.0
            elif self.increment < 1e-05:
                self.increment = 1e-05

            self.ui.spinBox_increment.setValue(self.increment)
            self.set_increment()
            #print(step, self.increment)

    def keyboard_next_spectra(self):
        value = self.ui.select_datanumber.value()
        self.ui.select_datanumber.setValue(value + 1)

    def keyboard_prev_spectra(self):
        value = self.ui.select_datanumber.value()
        self.ui.select_datanumber.setValue(value - 1)

    def fit_log(self,*args):
        if args[0] == 2:    # if chekbox is checked, open new window
            self.flw = Fit_Log()
            self.flw.setWindowTitle("Fit Log")
            self.flw.setGeometry(0, 0, 600, 720)
            self.flw.show()
            try:
                self.flw.setText(self.fit_report_log)
            except:
                pass
        else:   # close window with uncheking of checkbox
            self.flw.close()

    def set_increment(self):
        try:
            if self.index_model == 2:
                print(self.increment)
                for zähler in range(0,self.fit_num*3+2):
                    self.ui.Parameter_table.cellWidget(zähler, 0).setSingleStep(self.increment)
                    self.ui.Parameter_table.cellWidget(zähler, 1).setSingleStep(self.increment)
                    self.ui.Parameter_table.cellWidget(zähler, 2).setSingleStep(self.increment)
            elif self.index_model == 3:
                for zähler in range(0,self.fit_num*4+2):
                    self.ui.Parameter_table.cellWidget(zähler, 0).setSingleStep(self.increment)
                    self.ui.Parameter_table.cellWidget(zähler, 1).setSingleStep(self.increment)
                    self.ui.Parameter_table.cellWidget(zähler, 2).setSingleStep(self.increment)
                    self.ui.Parameter_table.cellWidget(zähler, 3).setSingleStep(self.increment)
        except:
            return 0

    def make_preB_sim(self):
        global value_opt
        self.preB_Sim = self.python_submit(True)
        # if no pre is present, then ani_pre_fit is True
        # ani_pre_fit determines whether it is the first time fitting or just a shift
        if not self.preB_Sim == None:
            if value_opt['ani_pre_fit']:
                #value_opt['ani_pre_fit'] = False
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
        #try:
        if value_opt['ani_pre_fit']:
            self.ani_fit_angle = False
            #print(d,value_opt['ani_pre_fit'])
            #self.make_preB_sim()
            self.update_canvas(self.preB_Sim,self.preB_Sim[2])
            self.changeing_phi = np.copy(self.preB_Sim[2])
            value_opt['ani_pre_fit'] = False
            # call function to generate Sim to find good shift
        else:
            # update plot according to shift d,
            # by calculating new phi_shifted, then plot
            self.changeing_phi = np.copy(self.preB_Sim[2])
            angle = make_phirange(d, self.changeing_phi, True, self.minWinkel, self.maxWinkel)
            self.update_canvas(self.preB_Sim, angle)
            self.ani_fit_angle = np.copy(angle)
        #except Exception as e:
        #    print('Error in get_shift: ',e)

    def update_canvas(self,data:list,angle:list):
        #Data is 2D List: data[0] = B_Sim, data[1] = B_Exp
        # Mathematica Canvas
        self.ui.Ani_Const_Plot.canvas.ax.clear()
        self.ui.Ani_Const_Plot.canvas.ax.set_ylabel('Resonance Field [T]')
        self.ui.Ani_Const_Plot.canvas.ax.set_xlabel('Angle [Deg]')

        #print(angle)

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

    def button_dropped_points(self):
        #as something is added to the exceptions, this is beeing called
        excep_append = self.Exceptions()  #gets text of editLine
        excep_append.append(str(self.ui.select_datanumber.value())) #get value of scrollbar
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

    def block_spinbox_signal(self,block):
        # block: bool

        # Temporarily block every signal from the Spinboxes. There is the Problem,
        # that everytime the value from these spinboxes is changed a signal is emited.
        # The bad thing now is, that there is no differentiation between user_event and machine_event,
        # therefore if the script changes the values to display the next angle,
        # the .connect() function will be called len(params) times, which has a min of 5 times with up to 42 times.

        try:
            for i in range(0, self.fit_num * value_opt['index_model_num'] + 2):
                self.ui.Parameter_table.cellWidget(i, 0).blockSignals(block)
        except:
            for i in range(0, self.fit_num * (self.index_model + 1) + 2):
                self.ui.Parameter_table.cellWidget(i, 0).blockSignals(block)

    def set_dataslider(self):
        self.ui.select_datanumber.setValue(int(self.ui.Scroll_Bar_dropped_points.value()))

    def set_datanumber(self):
        #Defines the dataset looked at
        global i
        global Bdata2
        global Adata2
        global value_opt
        if value_opt['data'] == 'loaded':
            if int(self.ui.select_datanumber.value()) > len(Adata)-1:
                i = len(Adata)-1
                self.ui.select_datanumber.setValue(int(i))
            else:
                i = int(self.ui.select_datanumber.value())
            self.i = i
            Bdata2 = Bdata[i]
            Adata2 = Adata[i]
            self.plot_data(i)
            spectra = self.dyn_params_table[1][i]
            if spectra != None and spectra != False:
                self.set_default_values(False)
        else:
            self.openFileDialog()

    def update_bar(self):
        #updates the progressbar
        self.ui.label_7.setText(str(self.progress))
        self.progress += 1
        self.ui.progressBar.setValue(self.progress)

    def set_model_type_number(self):
        # sets the integer number of index_model
        global value_opt
        if self.index_model == 2:
            value_opt['index_model_num'] = 3
        elif self.index_model == 3:
            value_opt['index_model_num'] = 4

    def model_type(self):
        # selects the type of function used in the fit
        self.index_model = self.ui.comboBox_fit_model.currentIndex()
        self.set_model_type_number()

    def select_fit_number(self,d):
        # define the number of function that are beeing fitted
        # Parameter d is value of the SpinBoxand is passed when SpinBox calls this function
        self.fit_num = d
        if self.sanity:
            if d >= 10:
                choice = QMessageBox.question(self, 'Sanity Check',
                                              "You are about to create " + str(d) + " Fit functions!\n" +
                                              "Was it your intention, or do you want to see your Computer burn?!\n " +
                                              "Continue?")
                if choice == 65536: # False
                    print(choice)
                    print('Reset Fit_Num')
                    self.d = 1
                    self.ui.spinBox_fit_num.setValue(1)
                else: # True
                    self.sanity = False

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
            if self.index_model == 2:
                #table for lorentz (Fitting)
                self.ui.Parameter_table.setRowCount(self.fit_num*3+2)
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

                for index in range(1,self.fit_num+1):
                    for list_index in parameter_Table_names_L:
                        parameter_table_names_L_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_L_final)

                for zähler in range(0,self.fit_num*3):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox()) # initial_values
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox()) # bound_min
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox()) # bound_max
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox()) # use for fitting ?

                self.set_default_values(True)

            else:
                #table for dyson (at the moment!) (Fitting)
                self.ui.Parameter_table.setRowCount(self.fit_num*4+2)
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

                for index in range(1,self.fit_num+1):
                    for list_index in parameter_Table_names_D:
                        parameter_table_names_D_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_D_final)

                for zähler in range(0,self.fit_num*4):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox())
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox())

                self.set_default_values(True)
        except Exception as e:
                print('Error in make_parameter_table',e)
                #Assuming a lorentz func as long as index_model hasnt been set
                self.ui.Parameter_table.setRowCount(self.fit_num*3+2)
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

                for index in range(1,self.fit_num+1):
                    for list_index in parameter_Table_names_L:
                        parameter_table_names_L_final.append(list_index+str(index))

                self.ui.Parameter_table.setVerticalHeaderLabels(parameter_table_names_L_final)

                for zähler in range(0,self.fit_num*3):
                    self.ui.Parameter_table.setCellWidget(zähler+2, 0, QDoubleSpinBox()) # initial_values
                    self.ui.Parameter_table.setCellWidget(zähler+2, 1, QDoubleSpinBox()) # bound_min
                    self.ui.Parameter_table.setCellWidget(zähler+2, 2, QDoubleSpinBox()) # bound_max
                    self.ui.Parameter_table.setCellWidget(zähler+2, 3, QCheckBox()) # use for fitting ? 

                self.set_default_values(True)

    def set_default_values(self,init):
        #print("DEBUG_SET_DEFAULT_VALUES")
        # sets default values into the spinbox according to arrays defined in the beginning
        # its basicly is just a for loop in order to catch every row of the table according to the number of lines

        # init: bool
        try:
            self.block_spinbox_signal(True)
            if init:
                # Called after functions were generated (function above)
                incr = self.ui.spinBox_incr.value()
                if incr > 0 :
                    self.Arrays = Gen_array(self.fit_num, increase=True, increment=incr)
                else:
                    self.Arrays = Gen_array(self.fit_num)
            else:
                # Called, when fitted and angle is changed
                self.Arrays = Gen_array(self.fit_num, increase=False, params=self.dyn_params_table[1][self.i])

            if self.index_model == 2:    #Lorentz
                for zahl in range(0,2): # Linear Funktion, goto next for loop, for lorentz parameters
                    self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-1000)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])

                    self.ui.Parameter_table.cellWidget(zahl,0).setSingleStep(self.Arrays.default_stepsize_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setSingleStep(self.Arrays.default_stepsize_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setSingleStep(self.Arrays.default_stepsize_linear[zahl])

                    self.ui.Parameter_table.cellWidget(zahl,0).setValue(self.Arrays.default_linear[zahl]) # Inital value
                    self.ui.Parameter_table.cellWidget(zahl,1).setValue(self.Arrays.default_boundaries_linear_min[zahl]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zahl,2).setValue(self.Arrays.default_boundaries_linear_max[zahl]) # Boundary maximum

                    if init:
                        self.ui.Parameter_table.cellWidget(zahl, 0).valueChanged.connect(self.signal_spinbox_manual_params) #different approach to connect a signal to the spinbox


                for zähler in range(0,self.fit_num*3): # for loop for Lorentz parameters
                    self.ui.Parameter_table.cellWidget(zähler+2,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setSingleStep(self.Arrays.default_stepsize_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setSingleStep(self.Arrays.default_stepsize_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setSingleStep(self.Arrays.default_stepsize_L[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setMaximum(self.Arrays.default_maximum_bound_spinbox_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setMaximum(self.Arrays.default_maximum_bound_spinbox_L[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setMaximum(self.Arrays.default_maximum_bound_spinbox_L[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(self.Arrays.default_values_L[zähler]) # Inital value
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(self.Arrays.default_boundaries_L_min[zähler]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(self.Arrays.default_boundaries_L_max[zähler]) # Boundary maximum
                    if init:
                        self.ui.Parameter_table.cellWidget(zähler+2, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox

            else:   # Dyson
                for zahl in range(0,2): # Linear Funktion
                    self.ui.Parameter_table.cellWidget(zahl,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zahl,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,1).setMinimum(-1000)
                    self.ui.Parameter_table.cellWidget(zahl,2).setMinimum(-1000)

                    self.ui.Parameter_table.cellWidget(zahl,0).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setMaximum(self.Arrays.default_maximum_bound_spinbox_linear[zahl])

                    self.ui.Parameter_table.cellWidget(zahl,0).setSingleStep(self.Arrays.default_stepsize_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,1).setSingleStep(self.Arrays.default_stepsize_linear[zahl])
                    self.ui.Parameter_table.cellWidget(zahl,2).setSingleStep(self.Arrays.default_stepsize_linear[zahl])

                    self.ui.Parameter_table.cellWidget(zahl,0).setValue(self.Arrays.default_linear[zahl]) # Inital value
                    self.ui.Parameter_table.cellWidget(zahl,1).setValue(self.Arrays.default_boundaries_linear_min[zahl]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zahl,2).setValue(self.Arrays.default_boundaries_linear_max[zahl]) # Boundary maximum
                    if init:
                        self.ui.Parameter_table.cellWidget(zahl, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox

                for zähler in range(0,self.fit_num*4):  # for loop for Dyson Parameter
                    self.ui.Parameter_table.cellWidget(zähler+2,0).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setDecimals(4)
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setDecimals(4)

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setSingleStep(self.Arrays.default_stepsize_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setSingleStep(self.Arrays.default_stepsize_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setSingleStep(self.Arrays.default_stepsize_D[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setMaximum(self.Arrays.default_maximum_bound_spinbox_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setMaximum(self.Arrays.default_maximum_bound_spinbox_D[zähler])
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setMaximum(self.Arrays.default_maximum_bound_spinbox_D[zähler])

                    self.ui.Parameter_table.cellWidget(zähler+2,0).setValue(self.Arrays.default_values_D[zähler]) # Inital value
                    self.ui.Parameter_table.cellWidget(zähler+2,1).setValue(self.Arrays.default_boundaries_D_min[zähler]) # Boundary minimum
                    self.ui.Parameter_table.cellWidget(zähler+2,2).setValue(self.Arrays.default_boundaries_D_max[zähler]) # Boundary maximum
                    if init:
                        self.ui.Parameter_table.cellWidget(zähler + 2, 0).valueChanged.connect(self.signal_spinbox_manual_params)  # different approach to connect a signal to the spinbox
        except Exception as e:
            print("Error in set_default_values:",e)
        self.block_spinbox_signal(False)


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

            i = 0   # laufvariable = 0
            i_min = 0
            i_max = int(chunksize)
            self.i_min = i_min
            self.i_max = i_max

            #D_min = min(D[:,3]) #for old colourplot
            #D_max = max(D[:,3]) #for old colourplot

            Bdata = np.split(np.true_divide(np.array(D[:,1]),10000),chunksize)  #List of sublists, magnetic field
            Winkeldata_raw = np.split(np.array(D[:,2]),chunksize)   #List of sublists, angle data
            Adata = np.split(np.array(D[:,3]),chunksize)    #List of sublists, amplitude data
            Winkeldata = []
            for i in range(chunksize):
                Winkeldata.append(Winkeldata_raw[i][0])
            X = np.split(np.true_divide(np.array(D[:,1]),10000),chunksize) #List of sublists, magnetic field, for colour plot
            Y = np.split(np.array(D[:,2]),chunksize)    #List of sublists, angle data, for colour plot
            Z = np.split(np.array(D[:,3]),chunksize)     #List of sublists, amplitude data, for colour plot
            #X = D[:,1].reshape(72,1024)
            #Y = D[:,2].reshape(72,1024)
            #Z = D[:,3].reshape(72,1024)

            self.H_range = max(Bdata[0]) * 1000 # Value in mT
            #self.H_range = max(Bdata[0])
            self.WinkelMin = min(D[:, 2])
            self.WinkelMax = max(D[:, 2])
            self.B_min = min(D[:,1])/10000 # in Tesla
            self.B_max = max(D[:,1])/10000
            self.B_ratio = (n-1)/self.B_max

            try:
                if int(self.ui.select_datanumber.value()) > len(Adata)-1: #catches errors, slider for datanumber
                    i = len(Adata)-1
                    self.ui.select_datanumber.setValue(int(i))
                else:
                    i = int(self.ui.select_datanumber.value())
            except:
                print("Weird Bug, I know, just look at the first Tab while loading data :D")


            # B_min bis B_max


            self.ui.Plot_Indi_View.lr.setBounds([round(self.B_min),self.B_max]) # Set Range Boundaries for linearregionitem

            self.ui.select_datanumber.setMaximum(i_max)
            self.ui.progressBar.setMaximum(i_max - 1)
            self.ui.Scroll_Bar_dropped_points.setMaximum(i_max)

            value_opt['data'] = 'loaded'
            self.dyn_params_table = [[], []]  # New Parameter Table: [0] is angle, [1] corresponding fitted params
            for i in range(i_max):
                self.dyn_params_table[0].append([i,Winkeldata[i]])
                self.dyn_params_table[1].append(None)
            self.set_datanumber()


    def error_msg(self):
        self.ui.label_params_output.setText('No data fitted yet!')
        print('Fehler irgendwo, bestimmt im dyn fit')

    def parameter_plot(self):
        global value_opt
        try:
            if value_opt['dyn_fit'] == 'fitted':
                Params_Plot(Parameter_list,'eigen',self.index_model)
            elif value_opt['parameter'] == 'loaded':
                print('Unfitted parameters. Using from file instead!')
                Params_Plot(Parameter_from_text,'eigen',self.index_model)
            else:
                self.load_parameters_to_text()
                Params_Plot(Parameter_from_text,'eigen',self.index_model)
        except Exception as e:
            print('Error in parameter_plot: ',e)

    def get_fit_options_from_file(self,fname):
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
                print('Error trying to get Dropped_Points, assuming Droppedpoints as empty:', e)
                Droppen_points = None

            # Takes Element 4 of header_parts and converts it into an usable array: ['A','dB','R',....]
            Params_name = np.asarray(ast.literal_eval(header_parts[3].split('[')[1].split(']')[0]))

            self.index_model = Lineshape
            print(value_opt['index_model_num'])
            self.set_model_type_number()
            print(value_opt['index_model_num'])
            self.fit_num = int((len(Params_name)-2)/value_opt['index_model_num'])
        else:
            print('File was not created by this script!')

    def load_parameter_plot(self):
        print('Loading from file!')
        try:
            self.load_parameters_to_text()
            Params_Plot(Parameter_from_text,'load',self.index_model)
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

    def signal_spinbox_manual_params(self,*args):
        #print("DEBUG_SPINBOX_SIGNAL")
        # args is the value of the SpinBox, that was changed
        # Params array with slope and offset inside!! So : [slope,offset,....]

        #self.block_spinbox_signal(True)
        param = self.dyn_params_table[1][self.i]
        if param != None and param != False:
            for zaehler in range(0, self.fit_num * self.index_model + 3):
                param[zaehler+2].value = self.ui.Parameter_table.cellWidget(zaehler, 0).value()
            self.plot_data(self.i)
        #self.block_spinbox_signal(False)

    def plot(self):
        global value_opt
        global result
        value_opt['dyn'] = 'single'
        self.model_type()
        if value_opt['data'] == 'loaded':
            try:
                j_min,j = self.get_fit_region() # get fit region
                self.set_init_params()  #gets the initial parameters from the GUI
                result = Fit(self.index_model, self.fit_num,Adata2,Bdata2, j_min,j,self.init_values,bound_min,bound_max).fit(self.index_model,Adata2,Bdata2, j_min,j) #fit and save it as result
                value_opt['fit'] = 'fitted'
                self.ui.progressBar.setMaximum(i_max-len(exceptions))
                self.evaluate_min_max(Bdata2[j_min:j],Adata2[j_min:j])
                self.set_fit_params()
                self.add_params_to_table(self.i,self.index_model,self.fit_num,[ param for name, param in result.params.items()])
                self.plot_data(self.i,result.best_fit, 'Best fit')
                self.fit_report_log = result.fit_report()
                if self.ui.checkBox_fit_log.isChecked():
                    self.flw.setText(self.fit_report_log)
                #print(result.fit_report())
                #print(result.best_values)
            except Exception as e:
                print('Error in main.plot: ',e)
        else:
            self.openFileDialog()

    def evaluate_min_max(self,Bdata,Adata):
        ind_min = Bdata[int(np.unravel_index(np.argmin(Adata, axis=None), Adata.shape)[0])]
        ind_max = Bdata[int(np.unravel_index(np.argmax(Adata, axis=None), Adata.shape)[0])]
        Bres = ind_max+(ind_min-ind_max)/2 #Asuming a symmetric line
        self.ui.label_test.setText('Bres from min/max: ' + str(Bres))

    def get_fit_region(self):
        region = self.ui.Plot_Indi_View.lr.getRegion()
        return int(float(region[0]) * self.B_ratio), int(float(region[1]) * self.B_ratio)

    def add_params_to_table(self, angle_index, function_index, function_number, params):
        # angle: int; function_number == index_model: int; function_number == fit_num: int; params: dict
        self.dyn_params_table[1][angle_index] = [function_index,function_number]
        for i in params:
            self.dyn_params_table[1][angle_index].append(i)

    def plot_data(self,angle_index,*args):
        #print("DEBUG_PLOT_DATA")
        # Plots data into a matplotlib canvas created in "pyqtgraph" skript, its the viewport for the measurement

        # args[0] = best_fit
        # args[1] = label name

        j_min, j = self.get_fit_region()
        #self.block_spinbox_signal(True)
        self.ui.Plot_Indi_View.plt.clear()  # Delete previous data
        self.ui.Plot_Indi_View.plt_range.clear() # Delete previous data

        self.ui.Plot_Indi_View.plt.plot(Bdata2, Adata2,name='Experiment', pen=(255,255,255))    # Plot Experiment data
        self.ui.Plot_Indi_View.plt_range.plot(Bdata2, Adata2, pen=(255,255,255))

        self.ui.Plot_Indi_View.plt.addItem(self.ui.Plot_Indi_View.lr) # Plot Linear Range Select Item
        self.ui.Plot_Indi_View.plt.addItem(self.ui.Plot_Indi_View.label, ignoreBounds=True)

        spectra = self.dyn_params_table[1][angle_index]
        if spectra != None and spectra != False:
        #if spectra != None:
            try:
                # --------------------------self.dyn_params_table[1][angle_index] Structure:----------------------------
                        # This is a List mixed with dict:
                              # raw = self.dyn_params_table[1][angle_index]
                              # raw[0] = index_model ====> The Function that has been used e.g. Dyson/Lorentz
                              # raw[1] = fit_num ====> Number of functions used
                              # Rest will be dict: raw[i>1].value or raw[i>1].bounds or raw[i>1].name
                                    # Ordering is the same as everywhere, see "array_gen.py" file
                raw = spectra
                slope = raw[2].value
                offset = raw[3].value
                temp_param = []
                for param in range(4, len(raw)):
                    temp_param.append(raw[param].value)

                if raw[0] == 2: # Lorentz
                    label = 'Lorentz'
                elif raw[0] == 3: # Dyson
                    label = 'Dyson'
                else:
                    label = 'How did you do this?!'

                for i_num in range(1, raw[1] + 1): # Plot individual Lines
                    plt_single_func = Functions.single_func(Bdata2, slope, offset, self.index_model, temp_param, i_num)
                    self.ui.Plot_Indi_View.plt.plot(Bdata2, plt_single_func, name=label + str(i_num),pen=(255 - 10 * i_num, 10 * i_num, 5 * i_num))
                    self.ui.Plot_Indi_View.plt_range.plot(Bdata2, plt_single_func, pen=(255 - 10 * i_num, 10 * i_num, 5 * i_num))

                #self.data_plot.setData(x, args[0], name=args[1], pen=(255,0,0))
                #self.data_plot_range.setData(x,  args[0], pen=(255,0,0))
                #plt_single_func = Functions.single_func(Bdata2, slope, offset, index_model, params, i_num)
                if raw[1] > 1 and args:
                    x = Bdata2[j_min:j]
                    self.ui.Plot_Indi_View.plt.plot(x, args[0], name='Result', pen=(0,255,0)) # Plot Fit data
                    self.ui.Plot_Indi_View.plt_range.plot(x,  args[0], pen=(0,255,0)) # Plot Fit data

            except Exception as e:
                print('Error in main.plot_data: ', e)
       # self.block_spinbox_signal(False)

    def plot_params_to_plot_tab(self):
        # Todo: Find way to change colours of plots different functions
        params_raw = self.dyn_params_table[1]
        angles_raw = self.dyn_params_table[0]

        params = []
        angles = []

        for l, i in enumerate(params_raw):
            if i != None and i != False:
                params.append(i)
                angles.append(angles_raw[l][1])
        params = np.asarray(params)
        angle = np.asarray(angles)
        # [:, i], (i,2,...) will give a parameter as an array
        # i=2:slope; i=3:offset; ....

        pen = None
        symbol = 'd'
        symbolpen = None
        symbolsize = 10
        symbolBrush = (255, 255, 255, 100)

        self.ui.plot_params.plt_slope.plot(angle, params[:, 2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
        self.ui.plot_params.plt_offset.plot(angle, params[:, 3],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
        print(params[0])
        if self.index_model == 2: # Lorentz
            for i_num in range(1,self.fit_num + 1):
                self.ui.plot_params.plt_db.plot(angle,params[:,2+3*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
                self.ui.plot_params.plt_R.plot(angle,params[:,3+3*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
                self.ui.plot_params.plt_A.plot(angle,params[:,4+3*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
        elif self.index_model == 3: # Dyson
            for i_num in range(1, self.fit_num + 1):
                self.ui.plot_params.plt_alpha.plot(angle, params[:,2+4*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
                self.ui.plot_params.plt_db.plot(angle, params[:,3+4*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
                self.ui.plot_params.plt_R.plot(angle, params[:,4+4*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)
                self.ui.plot_params.plt_A.plot(angle, params[:,5+4*(i_num-1)+2],pen=pen, symbol=symbol, symbolPen=symbolpen, symbolSize=symbolsize, symbolBrush=symbolBrush)


    def Exit(self):
        sys.exit()  #Obviously not the start

    def saveFileDialog(self):
        #gets the filename for saving
        exceptions = self.Exceptions()
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Please select the file to save to","","All Files (*);;Text Files (*.dat)",options=options)
        if fileName:
            self.save_filename = fileName
            file, names = self.create_param_file()
            now = datetime.datetime.now()
            np.savetxt(self.save_filename, file, delimiter='\t', newline='\n', header='FMR-Fit\nThis data was fitted {} using: $.{} Lineshape  \nDropped Points {}     \nData is arranged as follows {}'.format(now.strftime("%d/%m/%Y, %H:%M:%S"), self.index_model,exceptions, names))

    def create_param_file(self):
        param_table = []
        name_table = []

        # create list of names
        for l in self.dyn_params_table[1]:
            if l != None and l != False:
                params = l[2:len(l)]
                for raw_name in params:
                    name_table.append(raw_name.name)
                name_table.append('angle')
                break

        # create list of params
        for list ,i in enumerate(self.dyn_params_table[1]):
            if i != None  and i != False:
                params = i[2:len(i)]
                temp = []
                for raw in params:
                    temp.append(raw.value)
                winkel = self.dyn_params_table[0][list]
                temp.append(winkel[1])
                param_table.append(temp)
        return param_table, name_table

    def start_worker(self):
            # started when Dyn Fit button pressed
            # then starts the worker responsible for dyn fit
            self.progress = 0
            j_min, j = self.get_fit_region()
            exceptions = self.Exceptions()

            # set every spectra that is excluded to False
            for l in range(self.i_min, self.i_max):
                if l in exceptions:
                    self.dyn_params_table[1][l] = False
                    self.dyn_params_table[0][l] = [l,False]
            self.ui.progressBar.setMaximum(self.i_max - 1 - len(exceptions))
            # Start second Thread to fit angular dependence
            self.get_thread = Worker(Bdata, Adata, self.i, value_opt['dyn'], i_min, i_max, j, j_min,
                                     self.dyn_params_table)
            self.get_thread.start()
            self.get_thread.i_signal.connect(self.update_bar)
            self.get_thread.error_dynfit_signal.connect(self.error_msg)
            self.get_thread.job_done.connect(self.dyn_fit_job_done)

    def dyn_fit_job_done(self,*args):
        # is called when worker is finished
        # refreshes dyn_params_table with new data
        self.dyn_params_table = args[0]
        self.plot_params_to_plot_tab()

    def reset_fit(self):
        # reset self.dyn_params_table to None
        for i in range(self.i_max):
            self.dyn_params_table[1][i] = None
        self.set_datanumber()
        self.progress = 0
        self.ui.progressBar.setValue(0)
        self.ui.label_7.setText('')
        # reset the plots
        self.ui.plot_params.plt_slope.clear()
        self.ui.plot_params.plt_offset.clear()
        self.ui.plot_params.plt_alpha.clear()
        self.ui.plot_params.plt_db.clear()
        self.ui.plot_params.plt_R.clear()
        self.ui.plot_params.plt_A.clear()

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
            if hasattr(self,"save_filename"):
                print("Python Submit")
                if not args[0]:
                    init_load(self.save_filename, F, ani_fit_params, ani_fixed_params, shift, anglestep, True, True)
                else:
                    result = init_load(self.save_filename, F, ani_fit_params, ani_fixed_params, shift, anglestep, False, False)
                    return result
            else:
                fileName, _ = QFileDialog.getOpenFileName(self, "Please select the file to load Data")
                self.save_filename = fileName
                self.make_preB_sim()
        except Exception as e:
            print('Error in main.python_submit(): ',e)
            #print('Try fitting the spectra first!')
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

        if self.index_model == 2:    #Lorentz
            BresColumn = 4
            WinkelColumn = 6
        elif self.index_model == 3:      #Dyson
            BresColumn = 5
            WinkelColumn = 7
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
            print('Error in main.mathematica_submit: ',e)

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

            if self.index_model == 2:
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
            print('Error in main.plot_ani_fit: ',e)


    def create_mathematica_array(self,text):
        #text is the text given by the lineEdits
        mat_array = "{" + text + "}"
        return mat_array

    def set_init_params(self):
        #as the name say's
        global params
        global i
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

        if self.index_model == 2:
            for zahl in range(0,2):
                init_lin = self.ui.Parameter_table.cellWidget(zahl,0).value() # Inital value linear
                bounds_min_lin = self.ui.Parameter_table.cellWidget(zahl,1).value() # Boundary minimum linear
                bounds_max_lin = self.ui.Parameter_table.cellWidget(zahl,2).value() # Boundary maximum linear

                self.Arrays.default_linear[zahl] = init_lin
                self.Arrays.default_boundaries_linear_min[zahl] = bounds_min_lin
                self.Arrays.default_boundaries_linear_max[zahl] = bounds_max_lin

            for zähler in range(2,self.fit_num*3+2):
                init = self.ui.Parameter_table.cellWidget(zähler,0).value() # Inital value
                bounds_min = self.ui.Parameter_table.cellWidget(zähler,1).value() # Boundary minimum
                bounds_max =self.ui.Parameter_table.cellWidget(zähler,2).value() # Boundary maximum

                self.Arrays.default_values_L[zähler-2] = init
                self.Arrays.default_boundaries_L_min[zähler-2] = bounds_min
                self.Arrays.default_boundaries_L_max[zähler-2] = bounds_max

            self.init_values = np.append(self.Arrays.default_linear, self.Arrays.default_values_L)
            bound_min = np.append(self.Arrays.default_boundaries_linear_min, self.Arrays.default_boundaries_L_min)
            bound_max = np.append(self.Arrays.default_boundaries_linear_max, self.Arrays.default_boundaries_L_max)

        else:
            for zahl in range(0,2):
                init_lin = self.ui.Parameter_table.cellWidget(zahl,0).value() # Inital value linear
                bounds_min_lin = self.ui.Parameter_table.cellWidget(zahl,1).value() # Boundary minimum linear
                bounds_max_lin = self.ui.Parameter_table.cellWidget(zahl,2).value() # Boundary maximum linear

                self.Arrays.default_linear[zahl] = init_lin
                self.Arrays.default_boundaries_linear_min[zahl] = bounds_min_lin
                self.Arrays.default_boundaries_linear_max[zahl] = bounds_max_lin

            for zähler in range(2,self.fit_num*4+2):
                init = self.ui.Parameter_table.cellWidget(zähler,0).value() # Inital value
                bounds_min = self.ui.Parameter_table.cellWidget(zähler,1).value() # Boundary minimum
                bounds_max = self.ui.Parameter_table.cellWidget(zähler,2).value() # Boundary maximum

                self.Arrays.default_values_D[zähler-2] = init
                self.Arrays.default_boundaries_D_min[zähler-2] = bounds_min
                self.Arrays.default_boundaries_D_max[zähler-2] = bounds_max

            self.init_values = np.append(self.Arrays.default_linear, self.Arrays.default_values_D)
            bound_min = np.append(self.Arrays.default_boundaries_linear_min, self.Arrays.default_boundaries_D_min)
            bound_max = np.append(self.Arrays.default_boundaries_linear_max, self.Arrays.default_boundaries_D_max)

        if int(self.ui.select_datanumber.value()) > len(Adata)-1:
            i = len(Adata)-1
            self.ui.select_datanumber.setValue(int(i))
        else:
            i = int(self.ui.select_datanumber.value())
        #self.i = i

    def set_fit_params(self):
        j_min, j = self.get_fit_region()
        if value_opt['fit'] == 'fitted':
            if self.index_model == 2:
                temp_paras = Fit(self.index_model,self.fit_num,Adata2,Bdata2, j_min,j,self.init_values,bound_min,bound_max).give_params(self.fit_num, parameter_table_names_L_final,self.index_model,Adata2,Bdata2, j_min,j) #grabs the params file from different class, Lorentz
            else:
                temp_paras = Fit(self.index_model,self.fit_num,Adata2,Bdata2, j_min,j,self.init_values,bound_min,bound_max).give_params(self.fit_num,parameter_table_names_D_final,self.index_model,Adata2,Bdata2, j_min,j) #grabs the params file from different class, Dyson
            self.refresh_inital_parameter(temp_paras)
        else:
            self.plot() #....

    def refresh_inital_parameter(self,temp_paras):
        self.block_spinbox_signal(True)
        if self.index_model == 2:
            for zähler in range(0,self.fit_num*3+2):
                self.ui.Parameter_table.cellWidget(zähler,0).setValue(temp_paras[zähler]) # Inital value
        else:
            for zähler in range(0,self.fit_num*4+2):
                self.ui.Parameter_table.cellWidget(zähler,0).setValue(temp_paras[zähler]) # Inital value
        self.block_spinbox_signal(False)

if __name__=="__main__":
    #appctxt = ApplicationContext() #
    multip.set_start_method("spawn")
    define_value_opt()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MyForm()
    w.show()
    #exit_code = appctxt.app.exec_()  #
    #sys.exit(exit_code) #
    sys.exit(app.exec_())
