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
import time

def shifting(arschloch,shift,deg):
    if deg == 'deg':
        for l,i in enumerate(arschloch):
            arschloch[l] += shift
        return arschloch
    else:
        shift = shift*m.pi/180
        for l,i in enumerate(arschloch):
            arschloch[l] += shift
        return arschloch

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

K2s_WERT = 275220#120000#
K2p_WERT = 277#1200#
K4s_WERT = 0#0
K4p_WERT = 7500#2000#
phiU_WERT = 1.8#2.13#
#phiB_WERT = 2*m.pi

maxBresDelta = 0.01
               
shift = 26 # in degree

#consts Einträge
consts = np.zeros(9)
consts[0] = 9.4567*10**9 #Frequenz
consts[1] = 2.05 #g
consts[2] = 1.282*10**6 #M
consts[3] = 4*m.pi*10**(-7) #mu0
consts[4] = 9.27*10**(-24) #mub
consts[5] = (6.62606957*10**(-34))/(2*m.pi) #hbar
consts[6] = consts[1]*consts[4]/consts[5]
consts[7] = m.pi/2
consts[8] = 2*m.pi*consts[0]#omega bzw. w in der Formel

#params Einträge
params = np.zeros(5)
params[0] = K2s_WERT
params[1] = K2p_WERT
params[2] = K4s_WERT
params[3] = K4p_WERT
params[4] = phiU_WERT

#Laden der gefitteten Resonanzpositionen
D = np.loadtxt("parameter-IP.dat",dtype='float',skiprows=4)
R_raw = D[:,2]
Winkel = D[:,6]
Winkel2 = np.copy(Winkel)
Winkel2 = shifting(Winkel2,shift,'deg')
R_inter = interp1d(Winkel, R_raw)


#Freie Energy Formel + Smit-Beljers-Suhl approach
#Evtl Baselgia approach verwenden!!!!-------------------------------
F = M*B*(sin(theta)*sin(thetaB)*cos(phi-phiB)+cos(theta)*cos(thetaB))-(1/2*mu0*M**2-K2s)*sin(theta)**2-K2p*sin(theta)**2*cos(phi-phiU)**2-1/2*K4s*cos(theta)**4-1/8*K4p*(3+cos(4*phi))*sin(theta)**4
#Smit-Beljers-Suhl nur Gültig, wenn theta != 0
#Ansonsten Baselgia verwenden oder gleich nur Baselgia

#Smit-Beljers-Suhl:
#halb_res = gamma**2/(M*sin(theta))**2*(F.diff(theta,2)*F.diff(phi,2) - F.diff(theta,phi)**2)

#Baselgia:
halb_res = gamma**2/M**2*(F.diff(theta,2)*(F.diff(phi,2)/(sin(theta))**2 + cos(theta)/sin(theta)*F.diff(theta))-(F.diff(theta,phi)/sin(theta)-cos(theta)/sin(theta)*F.diff(phi)/sin(theta))**2)

eq_halb_res = Eq(w**2,halb_res)
#Formel der Resonanzfrequenz B_RES, Lösung [1] ist positive Lösung Funktion für Baselgia; Lösung [0] ist positive Lösung Funktion für Smit-Beljers-Suhl
B_RES = solve(eq_halb_res,B)

Bounds = [(0,2*m.pi),(0,3*m.pi)]    #Zum minimieren



def replace():
    P = B_RES.args[1]
    a = P.subs({g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],theta:consts[7],w:consts[8]})
    b = a.subs({K2s:params[0],K2p:params[1],K4s:params[2],K4p:params[3],phiU:params[4]})
    print(params)
    #In b fehlen noch alle Winkel: PhiB, ThetaB, Phi
    return b

def replace_ani_fit(a,b,c):
    #a = Array phi_Range, b = Array EqAngles[:,0], c = Array EqAngles[:,1]
    P1 = B_RES.args[1]
    A1 = P1.subs({g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],theta:consts[7],w:consts[8]})
    #B1 = A1.subs({phiB:c,thetaB:b})
    #funktion, in der alle Ani konstanten Fehlen
    return A1


def fit_funkt(phi_real,K2S,K2P,K4S,K4P,PHI_U):
    global PF1
    fkt = []
    P0 = B_RES.args[1]
    PF = F.subs({g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],theta:consts[7],w:consts[8]}) 
    P1 = P0.subs({g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],theta:consts[7],w:consts[8]})
    print(K2S,K2P,K4S,K4P)  
    for l,i in enumerate(phi_real):
        start_Paras = [Angles_Start[0][l],Angles_Start[1][l]]
        P2 = P1.subs({phi:phi_real[l],K2s:K2S,K2p:K2P,K4s:K4S,K4p:K4P,phiU:PHI_U})
        PF1 = PF.subs({phi:phi_real[l],K2s:K2S,K2p:K2P,K4s:K4S,K4p:K4P,phiU:PHI_U,B:B_RANGE[l]})
        #print(P2)
        minim = minimize(minim_funktion,start_Paras,method = 'TNC')
        P3 = P2.subs({thetaB:minim.x[0],phiB:minim.x[1]})
        fkt.append(float(P3))

    #print(K2S)
    return fkt


def calcB_things(a,b,c):
    #a ist Array vom Winkel PhiB, b ist Array vom Winkel phi, c ist Array für ThetaB
    X = B_origin
    Y = X.subs({phiB:a,phi:b,thetaB:c})
    #print(a,b,c)
    return Y

def create_B_RANGE(a):
    return R_inter(a*180/m.pi)

def minim_funktion(x):
    alt_fkt_temp = PF1.subs({phiB:x[1],thetaB:x[0]})
    return float(alt_fkt_temp)

def zu_minimierende_funktion(x):
    alt_fkt_temp = test_fkt.subs({phiB:x[1],thetaB:x[0]})
    return alt_fkt_temp

def alternative_fkt(B1,phi_RANGE):
    a1 = F.subs({g:consts[1],M:consts[2],mu0:consts[3],mub:consts[4],hbar:consts[5],gamma:consts[6],theta:consts[7],w:consts[8]})
    b1 = a1.subs({K2s:params[0],K2p:params[1],K4s:params[2],K4p:params[3],phiU:params[4],B:B1,phi:phi_RANGE})
    #print(a1)
    return b1

def minimizeer(it):
    #Funktion die zum minimieren dient
    start_Paras = Angles_Start
    #result = minimize(zu_minimierende_funktion,start_Paras,method = 'SLSQP',bounds=Bounds)
    result = minimize(zu_minimierende_funktion,start_Paras,method = 'TNC')
    return result

def get_eq_angles(testetst,Angles_Start):
    global test_fkt
    #print(it_while)
    #print(Angles_Start[1])
    for it,los_anglos in enumerate(phi_RANGE):
        start_Paras = [Angles_Start[0][it],Angles_Start[1][it]]
        test_fkt = testetst[it]
        result = minimize(zu_minimierende_funktion,start_Paras,method = 'TNC',bounds=Bounds)
        #print('ThetaB:',result.x[0],';','PhiB:',result.x[1])
        EqAngles_func[0][it] = result.x[0]
        EqAngles_func[1][it] = result.x[1]
        #print('PhiB:',result.x[1],'\tPhi:',los_anglos,'\tEqAngles[1]:',EqAngles[1][it])
    return EqAngles_func


#Beginn der Simulation

phi_min = min(Winkel)*m.pi/180
phi_max = max(Winkel)*m.pi/180
#phi_min = 0
#phi_max = 359*m.pi/180
phi_step = 2*m.pi/90
phi_RANGE = np.arange(phi_min,phi_max,phi_step,dtype='float')
B_RANGER = np.vectorize(create_B_RANGE)


phi_shifted = np.copy(phi_RANGE)
phi_shifted = shifting(phi_shifted,shift,'rad')
B_RANGE = B_RANGER(phi_RANGE)

#Array für Equilibrium Angles:
EqAngles_func = np.ndarray((2,len(phi_shifted)))

B_origin = replace()
B1 = np.vectorize(calcB_things)
zu_min_fct = np.vectorize(zu_minimierende_funktion)
replace_for_ani_fit = np.vectorize(replace_ani_fit)

#------Neuster Versuch------------
#1. B1 berechnen mit der annahme, dass die Eq-Winkel = Winkel des externen Feldes
#2. Die echten Eq-Winkel aus dem errechneten B1 bestimmen
#3. Berechne neuen Wert B2 aus den neuen echten Eq-Winkeln
#4. Vergleiche B1 mit B2, wenn differenz kleiner als d_limit breche iteration ab, ansonsten mache nächste Iteration startend bei Schritt 2 mit B1=B2


#Start Parameter der Winkel zur Minimalisierung


Angles_Start = np.ndarray((2,len(phi_shifted)))
for _,winkel in enumerate(phi_shifted):
    Angles_Start[0][_] = 5*m.pi/2
    Angles_Start[1][_] = winkel
'''

Angles_Start = [[],[]]

Angles_Start[0] = m.pi/2#theta
Angles_Start[1] = m.pi#phi
'''

alt_fkt = np.vectorize(alternative_fkt) #zum setzen der zu minimierenden Funktion
B_1 = B1(phi_shifted,phi_shifted,m.pi/2)
B1_fkt = alt_fkt(B_1,phi_shifted) #Array wo B1 eingesetzt wurde, aber thetaB und phiB fehlen!

minim = np.vectorize(minimizeer)

params_Fit = Parameters()
params_Fit.add_many( 
            ('K2S',params[0], True, -400000, 400000, None, None),             
            ('K2P',params[1], True, -2000, 50000, None, None),
            ('K4S',params[2], False, -400000, 400000, None, None),             
            ('K4P',params[3], True, -50000, 50000, None, None),
            ('PHI_U',params[4], True, 0, 2*m.pi, None, None)
            )
model = Model(fit_funkt)


it_while = 0
start = time.time()
while True:
    start_it = time.time()
    print('1. Starting iteration',it_while+1)
    B1_fkt = alt_fkt(B_RANGE,phi_shifted)
    EqAngles = get_eq_angles(B1_fkt,Angles_Start) #Berechnet die gleichgewichts Winkel
    #Nun muss B2 aus diesen gleichgew. Winkeln berechnet werden.
    B_2 = B1(EqAngles[1],phi_shifted,EqAngles[0])
    #print('Die Summe von Eq und Phi_range',np.sum(EqAngles[1]-phi_shifted))
    B2_fkt = alt_fkt(B_2,phi_shifted)
    print('2.   ',abs(np.sum(B_1-B_2)))
    ''' 
    if abs(np.sum(B_1-B_2)) < maxBresDelta:
        print('Finished after maxBresDelta was reached!. \nIt took,',it_while,'Iterations')
        fkt_fur_fit = replace_for_ani_fit(phi_shifted,EqAngles[0],EqAngles[1])  #Erstelle Funktion für Fit
        print('Now fitting and recalculating the eq angles!')        
        ani_result = model.fit(B_RANGE, params_Fit, phi_real=phi_shifted) #Führe Fit durch

        #Nun müssen die K´s in B1 gesetzt werden
        gesamtWert_alt = params[0] + params[1] + params[2] + params[3] + params[4]
        params[0] = ani_result.params['K2S'].value
        params[1] = ani_result.params['K2P'].value
        params[2] = ani_result.params['K4S'].value
        params[3] = ani_result.params['K4P'].value
        params[4] = ani_result.params['PHI_U'].value
        params_Fit['K2S'].set(value=ani_result.params['K2S'].value,min=None,max=None)
        params_Fit['K2P'].set(value=ani_result.params['K2P'].value,min=None,max=None)
        params_Fit['K4S'].set(value=ani_result.params['K4S'].value,min=None,max=None)
        params_Fit['K4P'].set(value=ani_result.params['K4P'].value,min=None,max=None)
        params_Fit['PHI_U'].set(value=ani_result.params['PHI_U'].value,min=None,max=None)

        gesamtWert_neu = params[0] + params[1] + params[2] + params[3] + params[4]
        print('Alt vs Neu = ',gesamtWert_alt - gesamtWert_neu)
        if abs(gesamtWert_alt - gesamtWert_neu) == 0:
            break
        B_origin = replace()

        end_it = time.time()
        print('Iteration',it_while,'took :',end_it-start_it,'seconds')
    '''
    if it_while == 2:
        end_it = time.time()
        print('Iteration',it_while,'took :',end_it-start_it,'seconds')
        print('Reached maximum iterations')
        break
    else:
        Angles_Start = np.ndarray((2,len(phi_shifted)))
        Angles_Start = np.copy(EqAngles)

    it_while += 1
    end_it = time.time()
    print('Iteration',it_while,'took :',end_it-start_it,'seconds')

#print(EqAngles[:,1])
end = time.time()
print('ENDE !!!!! It took ',end-start,'seconds to calculate')


#plt.plot(phi_RANGE,B_2, '-b', label='Daten')
#plt.xlabel('Angle [Degrees]')
#plt.ylabel('B-Field [Arb. Units]')
#plt.show()


#fkt_fur_fit = replace_for_ani_fit(phi_shifted,EqAngles[:,0],EqAngles[:,1])

ani_result = model.fit(B_RANGE, params_Fit, phi_real=phi_shifted)
print(ani_result.fit_report())

plt.plot(phi_RANGE,B_RANGE, '-b', label='Interpolated Data')
plt.plot(phi_RANGE,ani_result.best_fit, '-r', label='best fit')
plt.legend()
plt.show()

#Auf Surface dauert es ohne Multiprocessing 143 sek für it_while = 2!

