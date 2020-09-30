import os
import datetime
from PyQt5.QtCore import QThread, pyqtSignal
from subprocess import call

class Py2Mat(QThread):
    save_path_signal = pyqtSignal(str)
    def __init__(self,filename,F,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,BresColumn,WinkelColumn,Shift):
        self.filename = filename
        self.F = F
        self.Params = Params
        self.StartVal = StartVal
        self.Ranges = Ranges
        self.Steps = Steps
        self.fixedParams = fixedParams
        self.fixedValues = fixedValues
        self.anglestep = anglestep
        self.iterations = iterations
        self.BresColumn = BresColumn
        self.WinkelColumn = WinkelColumn
        self.Shift = Shift
        QThread.__init__(self)

    def init_folder(self):
        # This function should get the filename string from the main script and then checks if folders are present, if not they will be created
        #pfad is the path to the folder the user saved his parameters to
        pfad = os.path.dirname(self.filename)
        if not os.path.exists(pfad + '/Anisotropy'):
            print('Folder Anisotropy does not exist')
            print('Creating /Anisotropy ....')
            os.mkdir(pfad + '/Anisotropy')
        return pfad + '/Anisotropy/'

    def run(self):
        now = datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        pfad = self.init_folder() + str(now)
        os.mkdir(pfad)
        outputPath = pfad + "/"
        print(outputPath)
        cwd = os.getcwd()  # current working directory
        '''F = 'F[B_, thetaB_, phiB_, M_, theta_, phi_] = M*B (Sin[theta]*Sin[thetaB]*Cos[phi - phiB] + Cos[theta]*Cos[thetaB])' \
            ' - (1/2*mu0*M^2 - K2s) Sin[theta]^2 - K2p*Sin[theta]^2*Cos[phi - phiu]^2 - 1/2*K4s*Cos[theta]^4' \
            ' - 1/8*K4p (3 + Cos[4*phi]) Sin[theta]^4'
        FilePath = filename     #Kann von GUI gesetzt werden
        Params = 'fitParameters = {K2p, K2s, K4p, phiu}'    #LineEdit
        StartVal = 'startvalues ={863.25, 261345, 13720.6, 5.0756}' #LineEdit
        Ranges = '{0.5, 0.5, 0.5, 0.5}' #LineEdit
        Steps = '{0.5, 0.5, 0.5, 0.5}' #LineEdit
        fixedParams = '{omega, g, M, K4s}' #LineEdit
        fixedValues = '{2 Pi*9.8782*10^9, 2.05, 1.53*10^6, 0}' #LineEdit
        anglestep = 'Pi/10'
        iterations = 1
        outputPath = pfad
        BresColumn = 3
        WinkelColumn = 7
        Shift = 37'''

        print('Starting Mathematica Kernel')
        call(['wolframscript', '-file', cwd+'/Fit.wls', '{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}'.format(self.F,self.filename,self.Params,self.StartVal,self.Ranges,self.Steps,self.fixedParams,self.fixedValues,self.anglestep,self.iterations,outputPath,self.BresColumn,self.WinkelColumn,self.Shift)])
        self.save_path_signal.emit(outputPath)