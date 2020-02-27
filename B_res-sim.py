from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from symengine.lib.symengine_wrapper import solve
from sympy import sin,cos,sqrt,Eq,Symbol,Function,Matrix,re
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import numpy as np
import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import math as m


K2p = Symbol('K2p')
K2s = Symbol('K2s')
K4p = Symbol('K4p')
K4s = Symbol('K4s')
M = Symbol('M')
mu0 = Symbol('mu0')
mub = Symbol('mub')
B = Symbol('B')
theta = Symbol('theta')
thetaB = Symbol('thetaB')
phi = Symbol('phi')
phiU = Symbol('phiU')
phiB = Symbol('phiB')
gamma = Symbol('gamma')
w = Symbol('w')
g = Symbol('g')
hbar = Symbol('hbar')
f = Symbol('f')
F = Function('F')

K2s_WERT = -110000
K2p_WERT = 100
K4s_WERT = 0
K4p_WERT = 8000
phiU_WERT = 2.2
#phiB_WERT = 2*m.pi

maxBresDelta = 0.01
               

#consts Einträge
consts = np.zeros(8)
consts[0] = 9.8856*10**9 #Frequenz
consts[1] = 2.1 #g
consts[2] = 1.27*10**6 #M
consts[3] = 4*m.pi*10**(-7) #mu0
consts[4] = 9.27*10**(-24) #mub
consts[5] = (6.62606957*10**(-34))/(2*m.pi) #hbar
consts[6] = consts[1]*consts[3]/consts[5]
consts[7] = m.pi/2

#params Einträge
params = np.zeros(5)
params[0] = K2s_WERT
params[1] = K2p_WERT
params[2] = K4s_WERT
params[3] = K4p_WERT
params[4] = phiU_WERT

#Laden der gefitteten Resonanzpositionen
D = np.loadtxt("params.dat",dtype='float',skiprows=1)
R_raw = D[:,2]
Winkel = D[:,6]
R_inter = interp1d(Winkel, R_raw, kind='cubic')

#Freie Energy Formel + Smit-Beljers-Suhl approach
#Evtl Baselgia approach verwenden!!!!-------------------------------
F = M*B*(sin(theta)*sin(thetaB)*cos(phi-phiB)+cos(theta)*cos(thetaB))-(1/2*mu0*M**2-K2s)*sin(theta)**2-K2p*sin(theta)**2*cos(phi-phiU)**2-1/2*K4s*cos(theta)**4-1/8*K4p*(3+cos(4*phi))*sin(theta)**4
halb_res = gamma**2/(M*sin(theta))**2*(F.diff(theta,2)*F.diff(phi,2) - F.diff(theta,phi)**2)
eq_halb_res = Eq(w**2,halb_res)
#Formel der Resonanzfrequenz B_RES, Lösung [0] ist positive Lösung der Wurzelfunktion
B_RES = solve(eq_halb_res,w)


def replace():
    P = B_RES.args[0]
    a = P.subs({f:consts[0],g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],thetaB:consts[7]})
    b = a.subs({K2s:params[0],K2p:params[1],K4s:params[2],K4p:params[3],phiU:params[4]})
    #In b fehlen noch alle Winkel: PhiB, Theta, Phi
    return b

def calcB_things(a,b,c,B_type):
    #a ist Array vom Winkel PhiB, b ist Array vom Winkel ThetaB, B_type gibt später an ob B1 oder B2 berechnet werden soll 
    X = B_origin
    if B_type == B1:
        #Das B_FELD wurde hier gesetzt, erlaubt ?
        Y = B_origin.subs({phiB:a,phi:a,theta:b,B:c})
        return Y
    else:
        return(Y)

def create_B_RANGE(a):
    return R_inter(a*180/m.pi)

def zu_minimierende_funktion(x,):
    a = F.subs({f:consts[0],g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],thetaB:consts[7]})
    b = a.subs({K2s:params[0],K2p:params[1],K4s:params[2],K4p:params[3],phiU:params[4]})
    Y = b.subs({phiB:it,phi:x[1],theta:x[0],B:float(R_inter(it*180/m.pi))})
    return abs(Y)
#Beginn der Simulation

phi_min = min(Winkel)*m.pi/180
phi_max = max(Winkel)*m.pi/180
#phi_min = 0
#phi_max = 359*m.pi/180
phi_step = 2*m.pi/90
phi_RANGE = np.arange(phi_min,phi_max,phi_step,dtype='float')
B_RANGER = np.vectorize(create_B_RANGE)
B_RANGE = B_RANGER(phi_RANGE)
#Array für Equilibrium Angles:
EqAngles = np.zeros(shape=(len(phi_RANGE),2),dtype='float')
B_origin = replace()
#B1 = np.vectorize(calcB_things)


#------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!
#Jetzt ist die erste Resonanzposition berechnet und es kann iterativ weiter gehen. ABER: Das B_Feld wurde auf einen speziellen Wert gesetzt, nämlich B=9, hat also keinerlei bezug zur realität !!!!!!!!!!!
#------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!------WICHTIG!


#Erster Schritt B1 = B_RES(eqAngles = thetaB = pi/2,PhiB = array aus vielen Winkeln)
#Zweiter Schritt B2(theta,thetaB,phi,phiB) aus B1 berechnen, das heißt man bekommt 0=B1-....sin(theta)...
#Dritter Schritt dieses 0=... zu den Winkeln minimalisieren
#Vierter Schritt Aus neuen EqAngles B2 berechnen
#Fünfter Schritt B1-B2 bestimmen, wenn deltaB größer als limit, dann wiederhole iteration mit B1=B2

#Start Parameter der Winkel zur Minimalisierung
Angles_Start = np.zeros(2)
Angles_Start[0] = m.pi/2#theta
Angles_Start[1] = m.pi#phi


B2 = 1
fitted = 0
for it in phi_RANGE: 
    B1 = re(B_origin.subs({phiB:it,phi:it,theta:m.pi/2,B:float(R_inter(it*180/m.pi))}))
    start_Paras = Angles_Start
    #print(re(B1)-re(B2))
    while abs(re(B1)-re(B2))>maxBresDelta:
        B1 = B2
        Bounds = [(0,2*m.pi),(0,2*m.pi)]
        #cons=({'type':'ineq','fun':lambda x:abs(B1-B2)-maxBresDelta})
        #results = minimize(minimized_B_func,start_Paras,method = 'SLSQP',bounds=Bounds,constraints=cons)
        results = minimize(zu_minimierende_funktion,start_Paras,method = 'SLSQP',bounds=Bounds)
        B2 = re(B_origin.subs({phiB:it,phi:results.x[1],theta:results.x[0],B:float(R_inter(it*180/m.pi))}))
        Angles_Start[0] = results.x[0]
        Angles_Start[1] = results.x[1]
        print(B1-B2)
    #print(Angles_Start)




#Letzte Änderung: B_Range wurde hinzugefügt, als interpoliertes B-Feld. Der Plot sieht ungut aus!
#Zudem wird B1 bei manchen B-Values Imaginär, deswegen das abs()

#print(abs(B1(phi_RANGE,m.pi/2,B_RANGE,B1)))
'''
plt.plot(phi_RANGE,abs(B1(phi_RANGE,m.pi/2,B_RANGE,B1)), '-b', label='Daten')
plt.xlabel('Angle [Degrees]')
plt.ylabel('B-Field [Arb. Units]')
plt.show()
'''
