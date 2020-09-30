from os import getcwd
from subprocess import call

cwd = getcwd() #current working directory

F = 'F[B_, thetaB_, phiB_, M_, theta_, phi_] = M*B (Sin[theta]*Sin[thetaB]*Cos[phi - phiB] + Cos[theta]*Cos[thetaB]) - (1/2*mu0*M^2 - K2s) Sin[theta]^2 - K2p*Sin[theta]^2*Cos[phi - phiu]^2 - 1/2*K4s*Cos[theta]^4 - 1/8*K4p (3 + Cos[4*phi]) Sin[theta]^4'
FilePath = cwd+'/parameter.dat'     #Kann von GUI gesetzt werden
Params = 'fitParameters = {K2p, K2s, K4p, phiu}'    #LineEdit
StartVal = 'startvalues ={863.25, 261345, 13720.6, 5.0756}' #LineEdit
Ranges = '{0.5, 0.5, 0.5, 0.5}' #LineEdit
Steps = '{0.5, 0.5, 0.5, 0.5}' #LineEdit
fixedParams = '{omega, g, M, K4s}' #LineEdit
fixedValues = '{2 Pi*9.8782*10^9, 2.05, 1.53*10^6, 0}' #LineEdit
anglestep = 'Pi/10'
iterations = 1
outputPath = cwd+'/output/'
BresColumn = 3
WinkelColumn = 7
Shift = 37



call(['wolframscript', '-file', cwd+'/Fit-test.wls', '{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}'.format(F,FilePath,Params,StartVal,Ranges,Steps,fixedParams,fixedValues,anglestep,iterations,outputPath,BresColumn,WinkelColumn,Shift)])
