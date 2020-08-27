import sys
import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import math as m
import threading
from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from PyQt5.QtWidgets import QMainWindow, QApplication, QAction,QFileDialog,QProgressBar
from PyQt5.QtCore import *
from Fitprogramm import *
#from fbs_runtime.application_context.PyQt5 import ApplicationContext
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from fitting import Fit
from parameter_plot import Params_Plot

def dyson_func(B, alpha, dB, R, A, slope, offset):
        return ((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))+slope*B+offset

def lorentz_func(B,dB, R, A, slope, offset):
    #return(-(A*(2*(B-dB)))/(m.pi*R**3*((B-dB**2)/(R**2)+1)**2))+slope*B+offset
    return(-64*A*(B-R)/(9*dB*(1+(B-R)**2)/(dB*m.sqrt(3)/4)**2)**2+slope*B+offset)


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
Parameter_list = []
fit_num = 1

class Worker(QThread):
    i_signal = pyqtSignal()
    error_dynfit_signal = pyqtSignal()
    def __init__(self,Bdata2,Adata2,Winkeldata,i,fileName,dyn_value,params,i_min,i_max,dmodel,j,j_min,exceptions):
        QThread.__init__(self)

    def fit(self,l):
        Bdata2 = Bdata[l]
        Adata2 = Adata[l]
        result = dmodel.fit(Adata2[j_min:j], params, B=Bdata2[j_min:j])
        return result

    def run(self):
        global dyn_fit_value
        f = open(fileName,'w')
        if dyn_value == 'single':
            text = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(float(self.fit(i).params['alpha'].value),
                float(self.fit(i).params['dB'].value),
                float(self.fit(i).params['R'].value),
                float(self.fit(i).params['A'].value),
                float(self.fit(i).params['slope'].value),
                float(self.fit(i).params['offset'].value),
                float(Winkeldata[i][0])
                )
            f.write('Alpha\tLinewidth-dB[T]\tResonance-field[T]\tAmplitude[Arb.Units]\tSlope\tOffset\tAngle[Degree]\n' + text)
        else:
            Palpha = params['alpha']
            PdB = params['dB']
            PR = params['R']
            PA = params['A']
            Poffset = params['offset']
            Pslope = params['slope']
            f.write('Alpha\tLinewidth-dB[T]\tResonance-field[T]\tAmplitude[Arb.Units]\tSlope\tOffset\tAngle[Degree]\n')
            Parameter_list.clear()
            for l in range(i_min,i_max):
                if l not in exceptions:    
                    self.i_signal.emit()
                    a = float(self.fit(l).params['alpha'].value)
                    db = float(self.fit(l).params['dB'].value)
                    R = float(self.fit(l).params['R'].value)
                    A = float(self.fit(l).params['A'].value)
                    slope = float(self.fit(l).params['slope'].value)
                    offset = float(self.fit(l).params['offset'].value)
                    text = a,db,R,A,slope,offset,float(Winkeldata[l][0])
                    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(a,db,R,A,slope,offset,float(Winkeldata[l][0])))
                    Parameter_list.append(text)
                    Palpha.set(value=a,min=None,max=None)
                    PdB.set(value=db,min=None,max=None)
                    PR.set(value=R,min=None,max=None)
                    PA.set(value=A,min=None,max=None)
                    Pslope.set(value=slope,min=None,max=None)
                    Poffset.set(value=offset,min=None,max=None)
        dyn_fit_value = 'fitted'
        f.close()

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
        self.show()

    def reset_drop(self):
        global excep_append
        if self.ui.Dropped_points_edit.text() == '':
            excep_append.clear()
        else:
            x = str(self.ui.Dropped_points_edit.text())
            excep_append = list(map(int, x.split(',')))
   
    def append_dropped_points(self):
        global excep_append
        excep_append.append(str(self.ui.Scroll_Bar_dropped_points.value()))
        text = ','.join(map(str, excep_append))
        self.ui.Dropped_points_edit.setText(str(text))

    def Exceptions(self):
        global exceptions
        global excep_append
        print('Hallo')
        x = str(self.ui.Dropped_points_edit.text())
        exceptions = list(map(int, x.split(',')))

    def change_parameter_angle(self):
        global W
        W = self.ui.parameter_data_select.value()
        self.select_fitted_params()


    def set_dataslider(self):
        self.ui.select_datanumber.setValue(int(self.ui.Scroll_Bar_dropped_points.value()))

    def set_datanumber(self):
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
        cbar.set_label('Amplitude [Arb. U.]')
        plt.xlabel('Magnetic Field [T]')
        plt.ylabel('Angle [Degrees]')
        plt.show()
        plt.close()

    def update_bar(self):
        global p
        self.ui.label_7.setText(str(p))
        p +=1
        self.ui.progressBar.setValue(p)

    def dyn_fit(self):
        global dyn_value
        dyn_value = 'dynamic'
        self.saveFileDialog()

    def model_type(self,Bool):
        global dmodel
        global index_model
        index_model = self.ui.comboBox_fit_model.currentIndex()
        if index_model == 1:
            if Bool == True:
                print('Gaussian deriv. selected')
            dmodel = Model(dyson_func)
        elif index_model == 2:
            if Bool == True:
                print('Lorentzian deriv. selected')
            dmodel = Model(lorentz_func)
        elif index_model == 3:
            if Bool == True:
                print('Dyson deriv. selected')
            dmodel = Model(dyson_func)
        else:
            if Bool == True:
                print('Not supported')
            dmodel = Model(dyson_func)

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
                D = np.loadtxt(fname[0],dtype='float')
            except:
                D = np.loadtxt(fname[0],dtype='float',skiprows=2)
            counts = np.unique(D[:,2], return_counts=True)
            counts1 = counts[1]
            n = counts1[0]
            chunksize = int(len(D[:,0])/n)
            self.ui.Drop_Scalebar.setMaximum(n-1)
            i = 0
            i_min = 0
            i_max = int(chunksize)
            D_min = min(D[:,3])
            D_max = max(D[:,3])

            Bdata = np.split(np.true_divide(np.array(D[:,1]),10000),chunksize)
            Winkeldata = np.split(np.array(D[:,2]),chunksize)
            Adata = np.split(np.array(D[:,3]),chunksize)
            X = np.split(np.true_divide(np.array(D[:,1]),10),chunksize)
            Y = np.split(np.array(D[:,2]),chunksize)
            Z = np.split(np.array(D[:,3]),chunksize)
            if int(self.ui.select_datanumber.text()) > len(Adata)-1:
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

    def fit(self,i):
        result = dmodel.fit(Adata2[j_min:j], params, B=Bdata2[j_min:j])
        return result
        
    def error_msg(self):
        self.ui.label_params_output.setText('No data fitted yet!')
        print('Fehler irgendwo, bestimmt im dyn fit')

    def parameter_plot(self):
        if dyn_fit_value == 'fitted':
            Params_Plot(Parameter_list,'eigen')
        elif parameter_value == 'loaded':
            print('Unfitted parameters. Using from file instead!')
            Params_Plot(Parameter_from_text,'eigen')
        else:
            self.load_parameters_to_text()
            Params_Plot(Parameter_from_text,'eigen')

    def load_parameter_plot(self):
        print('Loading from file!')
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            Params_Plot(params_fname,'load')

    def load_parameters_to_text(self):
        global Parameter_from_text
        global parameter_value
        params_fname = QFileDialog.getOpenFileName(self, 'Open file','/home')
        if params_fname[0]:
            Parameter_from_text = np.loadtxt(params_fname[0],dtype='float',skiprows=1)
        parameter_value = 'loaded'


    def select_fitted_params(self):
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
        self.ui.parameter_plot_widget.canvas.ax.clear()
        self.ui.parameter_plot_widget.canvas.ax.set_xlabel('Magnetic Field [T]')
        self.ui.parameter_plot_widget.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
        self.ui.parameter_plot_widget.canvas.ax.plot(x, y,color='black',marker='o', label='Experimental data')
        plotdys = dyson_func(x,alpha_list[W],db_list[W],R_list[W],A_list[W],slope_list[W],offset_list[W])
        print(plotdys)
        self.ui.parameter_plot_widget.canvas.ax.plot(x,plotdys,'r--',label='Dyson')
        self.ui.parameter_plot_widget.canvas.ax.legend()
        self.ui.parameter_plot_widget.canvas.draw() 



    def plot(self):
        global fit_value
        global dyn_value
        dyn_value = 'single'
        self.model_type(False)
        if data_value == 'loaded':
            fit_value = 'fitted'
            self.set_init_params()
            self.ui.label_params_output.setText(self.fit(i).fit_report())
            plt.ion()
            plt.cla()
            plt.clf()
            plt.plot(Bdata2[j_min:j],Adata2[j_min:j], '-b', label='Experimental data')
            plt.plot(Bdata2[j_min:j],self.fit(i).best_fit, 'r--', label='Best fit')
            plt.xlabel('Magnetic Field [T]',fontsize=12)
            plt.ylabel('Amplitude [Arb. U.]',fontsize=12)
            plt.legend(fontsize=12)
            plt.show()
        else:
            self.openFileDialog()
        if index_model == 0:
            self.ui.label_params_output.setText('Please select a fitting model first')

    def plot_data(self):
            self.ui.plotView.canvas.ax.clear()
            x=Bdata2[j_min:j]
            y=Adata2[j_min:j]
            self.ui.plotView.canvas.ax.set_xlabel('Magnetic Field [T]')
            self.ui.plotView.canvas.ax.set_ylabel('Amplitude [Arb. U.]')
            self.ui.plotView.canvas.ax.plot(x, y)
            self.ui.plotView.canvas.draw()

    def plot_in_colour(self):
        from mayavi import mlab
        global tick_number
        global colormap
        tick_number = self.ui.colour_tick_edit.text()
        colormap = self.ui.colour_map_edit.text()
        mlab.mesh(X,Y,Z, colormap=colormap)
        mlab.show()
        '''self.get_thread2 = ColourPlot(X,Y,Z,colormap,tick_number,D_min,D_max,exceptions)
        self.get_thread2.start()'''

    def Exit(self):
        sys.exit()

    def saveFileDialog(self):
        global p
        global fileName
        p=0
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.dat)",options=options)
        if fileName:

            self.get_thread = Worker(Bdata2,Adata2,Winkeldata,i,fileName,dyn_value,params,i_min,i_max,dmodel,j,j_min,exceptions)
            self.get_thread.start()
            self.get_thread.i_signal.connect(self.update_bar)
            self.get_thread.error_dynfit_signal.connect(self.error_msg)

    def set_init_params(self):
        global params
        global i
        params = Parameters()
        # add with tuples:(NAME,VALUE,VARY,MIN ,MAX,EXPR,BRUTE_STEP)
        params.add_many(
        ('alpha', float(self.ui.init_alpha.text()), True, float(self.ui.alpha_con_min.text()), float(self.ui.alpha_con_max.text()), None, None),
        ('dB',float(self.ui.init_dB.text()), True, float(self.ui.dB_con_min.text()), float(self.ui.dB_con_max.text()), None, None),
        ('R',float(self.ui.init_R.text()), True, float(self.ui.R_con_min.text()), float(self.ui.R_con_max.text()), None, None),
        ('A',float(self.ui.init_A.text()), True, float(self.ui.A_con_min.text()), float(self.ui.A_con_max.text()), None, None),
        ('offset',float(self.ui.init_offset.text()), True, float(self.ui.offset_con_min.text()), float(self.ui.offset_con_max.text()), None, None),
        ('slope',float(self.ui.init_slope.text()), True, float(self.ui.slope_con_min.text()), float(self.ui.slope_con_max.text()), None, None)
        )
        self.ui.label_show_init_params.setText('Parameter set!')
        if int(self.ui.select_datanumber.text()) > len(Adata)-1:
            i = len(Adata)-1
            self.ui.select_datanumber.setText(str(i))
        else:
            i = int(self.ui.select_datanumber.text())

    def set_fit_params(self):
        if fit_value == 'fitted':
            self.ui.init_alpha.setValue(float(self.fit(i).params['alpha'].value))
            self.ui.init_dB.setValue(float(self.fit(i).params['dB'].value))
            self.ui.init_R.setValue(float(self.fit(i).params['R'].value))
            self.ui.init_A.setValue(float(self.fit(i).params['A'].value))
            self.ui.init_slope.setValue(float(self.fit(i).params['slope'].value))
            self.ui.init_offset.setValue(float(self.fit(i).params['offset'].value))
            self.set_init_params()
        else:
            self.plot()


if __name__=="__main__":
    #appctxt = ApplicationContext()
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    #exit_code = appctxt.app.exec_()
    #sys.exit(exit_code)
    sys.exit(app.exec_())
