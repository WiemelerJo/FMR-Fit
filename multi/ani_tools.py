from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from symengine.lib.symengine_wrapper import solve
from sympy import sin, cos, sqrt, Eq, Symbol, Function, Matrix, re
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import shgo, differential_evolution  # activate when using global optimization
from multiprocessing import Pool
from functools import partial
import numpy as np
import matplotlib as mp

mp.use('QT5Agg')
import matplotlib.pyplot as plt
import math as m
import time


def shifting(array, shift, deg):
    if deg == 'deg':
        for l, i in enumerate(array):
            # print('1.',array[l])
            array[l] += shift
            # print('2.',array[l])
        return array
    else:
        shift = shift * m.pi / 180
        for l, i in enumerate(array):
            array[l] += shift
        return array


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
omega = Symbol('omega')
w = Symbol('w')
g = Symbol('g')
hbar = Symbol('hbar')
f = Symbol('f')
F = Function('F')

K2s_WERT = 100000  # 120000#
K2p_WERT = -287  # 1200#
K4s_WERT = 1  # 0
K4p_WERT = 10000  # 2000#
phiU_WERT = 2.15  # 2.12#2.13#
# phiB_WERT = 2*m.pi

maxBresDelta = 0.01

shift = -26

# consts Einträge
consts = np.zeros(8)
consts[0] = 9.4567 * 10 ** 9  # Frequenz
consts[1] = 2.05  # g
consts[2] = 1.1890 * 10 ** 6  # M
consts[3] = 4 * m.pi * 10 ** (-7)  # mu0
consts[4] = 9.27 * 10 ** (-24)  # mub
consts[5] = (6.62606957 * 10 ** (-34)) / (2 * m.pi)  # hbar
consts[6] = consts[1] * consts[4] / consts[5]
consts[7] = 2 * m.pi * consts[0]  # omega bzw. w in der Formel

# params Einträge
params = np.zeros(5)
params[0] = K2s_WERT
params[1] = K2p_WERT
params[2] = K4s_WERT
params[3] = K4p_WERT
params[4] = phiU_WERT

# Laden der gefitteten Resonanzpositionen
D = np.loadtxt("parameter-IP-19-09-2020.dat", dtype='float', skiprows=0)
R_raw = D[:, 3]
Winkel = D[:, 5].flatten()
Winkel2 = np.copy(Winkel)
Winkel2 = shifting(Winkel2, shift, 'deg')
B_inter = interp1d(Winkel, R_raw)  # Interpoliertes B_res, array mit Länge len(Winkel)

# Freie Energy Formel + Baselgia
F = M * B * (sin(theta) * sin(thetaB) * cos(phi - phiB) + cos(theta) * cos(thetaB)) - (
        1 / 2 * mu0 * M ** 2 - K2s) * sin(theta) ** 2 - K2p * sin(theta) ** 2 * cos(
    phi - phiU) ** 2 - 1 / 2 * K4s * cos(theta) ** 4 - 1 / 8 * K4p * (3 + cos(4 * phi)) * sin(theta) ** 4
# F = M*B*(sin(theta)*sin(thetaB)*cos(phi-phiB)+cos(theta)*cos(thetaB))-(1/2*mu0*M**2-K2s)*sin(theta)**2-K2p*sin(theta)**2*cos(phi-phiU)**2-1/2*K4s*cos(theta)**4-1/8*K4p*(3+cos(4*phi))*sin(theta)**4
halb_res = gamma ** 2 / M ** 2 * (
        F.diff(theta, 2) * (F.diff(phi, 2) / (sin(theta)) ** 2 + cos(theta) / sin(theta) * F.diff(theta)) - (
        F.diff(theta, phi) / sin(theta) - cos(theta) / sin(theta) * F.diff(phi) / sin(theta)) ** 2)
eq_halb_res = Eq(omega ** 2, halb_res)
# Formel der Resonanzfrequenz B_RES, Lösung [0] ist positive Lösung der Wurzelfunktion
B_RES = solve(eq_halb_res, B)  # Diese Funktion ist ResField(...) aus Mathematica!
B_RES = B_RES.args[1]

Bounds = [(0.001, 2 * m.pi), (0.001, 2 * m.pi), (0.001, 2 * m.pi)]  # Zum minimieren

params_name = [K2s, K2p, K4s, K4p, phiU]
const_name = [f, g, M, mu0, mub, hbar, gamma, omega]
names = [params_name, const_name]
names = [val for sublist in names for val in sublist]
# names ist flattened[params,consts] : ----->  [K2s, K2p, K4s, K4p, phiU, f, g, M, mu0, muB, hbar, gamma, omega]

rule = [params, consts]
rule = [val for sublist in rule for val in sublist]
# rule ist flattened Werte der Namen

# rules ist dann array wo name mit wert zusammengefügt ist
rules = []
for i in range(len(rule)):
    rules.append([names[i], rule[i]])

rule_consts = []
for i in range(len(const_name)):
    rule_consts.append([const_name[i], consts[i]])

phi_min = min(Winkel) * m.pi / 180
phi_max = max(Winkel) * m.pi / 180
phi_step = m.pi / 91
phi_RANGE = np.arange(phi_min, phi_max, phi_step, dtype='float')  # array of Rad
phi_RANGE_deg = np.multiply(phi_RANGE, 180 / m.pi)

phi_shifted = np.copy(phi_RANGE)
phi_shifted = shifting(phi_shifted, shift, 'rad')
phi_shifted_deg = np.multiply(phi_shifted, 180 / m.pi)

maxWinkel = max(Winkel)
minWinkel = min(Winkel)

for i in range(len(phi_shifted_deg)):
    if phi_shifted_deg[i] > maxWinkel:
        phi_shifted_deg[i] -= maxWinkel
    elif phi_shifted_deg[i] < minWinkel:
        phi_shifted_deg[i] += maxWinkel


def make_phirange(shift_value: float, phi_array: list, deg: bool):
    if deg:
        for l, i in enumerate(phi_array):
            # print('1.',array[l])
            phi_array[l] += shift_value
            # print('2.',array[l])

        for i in range(len(phi_array)):
            if phi_array[i] > maxWinkel:
                phi_array[i] -= maxWinkel
            elif phi_array[i] < minWinkel:
                phi_array[i] += maxWinkel
    else:
        shift_value = shift_value * m.pi / 180
        for l, i in enumerate(phi_array):
            phi_array[l] += shift_value
        for i in range(len(phi_array)):
            if phi_array[i] > phi_max:
                phi_array[i] -= phi_max
            elif phi_array[i] < phi_min:
                phi_array[i] += phi_max
    return phi_array


def solveAngles(B, phiB, fkt):
    # Minimize Angles Theta, Phi, ThetaB
    # start_Paras for the moment these are constant: [m.pi/2, m.pi, m.pi/2]
    # start_Paras = [m.pi/2, m.pi, m.pi/2]
    # --------------------------------------Finding Global Minimum (Mathematica approach)---------------------------------------------
    # Extremly slow, slower than using a global solver in minimize function
    '''MIN = []
                RESULT = []
                for a in np.arange(-1,1,m.pi/8):
                   for b in np.arange(-1,1,m.pi/8):
                     for c in np.arange(-1,1,m.pi/8):
                        start_Paras = [a,b,c]
                        result = minimize(F_fkt, start_Paras, args=(B,phiB,fkt) ,method = 'L-BFGS-B', bounds=Bounds)
                        MIN.append([result.fun,result.x[0],result.x[1],result.x[2]])
                MIN = np.asarray(MIN)
                for index,value in enumerate(MIN[:,0]):
                    if value == min(MIN[:,0]):
                        print(index, "Found!")
                        print(MIN[index])
                        #RESULT.append([MIN[:,1],MIN[:,2],MIN[:,3]])
                        print(RESULT,'Global Minimum')'''

    # After evaluating some minima: Global Minimum search was everytime equivalent to local minima search (sometimes modulo 2*PI). This equivalence must not be true in every case!!!
    # ----------------------------------------------------------------------------------------------------------------------------------

    start_Paras = [m.pi / 2, m.pi, m.pi / 2]
    result = minimize(F_fkt, start_Paras, args=(B, phiB, fkt), method='L-BFGS-B',
                      bounds=Bounds)  # alternatively a global search can be performed using shgo or differential_evolution, look up scipy Docs.
    # result = shgo(F_fkt,Bounds, args=(B,phiB,fkt))  #Works as fast as local search (maybe even faster), but produces weird artifacts in plot!?
    # result = dual_annealing(F_fkt,Bounds, args=(B,phiB,fkt))   # Finds nice solution but takes years to compute one array with Anglestep = m.pi/40
    # result = differential_evolution(F_fkt,Bounds, args=(B,phiB,fkt)) #slower than shgo, but also nice solution
    ##print(result.x[0],result.x[1],result.x[2],"\"Local Minimum\"")
    return result.x[0], result.x[1], result.x[2]


def find_shift():
    print("Hi")


def ResFieldNumInp(rules, phi_val):
    # takes EqAngles and solves B_RES using interpolated data
    B = B_inter(phi_val * 180 / m.pi)
    phiB_val = phi_val

    B_FKT = B_RES.subs({i: l for i, l in rules})  # rules beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rules})  # rules beeing inserted, for minimizing
    Eq_angles = solveAngles(B, phiB_val, FKT)
    result = B_FKT.subs({theta: Eq_angles[0], phi: Eq_angles[1], thetaB: Eq_angles[2], phiB: phiB_val})
    return [float(result), float(Eq_angles[0]), float(Eq_angles[1]), float(Eq_angles[2])]


def F_fkt(x, *args):
    # Function of Free Energy Density, that should be minimized
    # *args is array:[B,phiB]
    # x are the fitted params
    Fit_fkt = args[2].subs({B: args[0], phiB: args[1], theta: x[0], phi: x[1], thetaB: x[2]})
    return float(Fit_fkt)


def Model_Fit_Fkt(phi_real, K2S, K2P, K4S, K4P, PHI_U):
    # is the function used for least square fitting
    # inputs are values for parameter
    # outputs array "fkt", an element in fkt is the value of the function at point x=phi_real, so f(x)
    fkt = []
    a = B_RES.subs({i: l for i, l in rule_consts})  # Parameter, B and EqAngle missing
    for i, l in enumerate(phi_real):
        b = a.subs({phiB: l, K2s: K2S, K2p: K2P, K4s: K4S, K4p: K4P, phiU: PHI_U})
        c = b.subs({theta: Eq_angles[0][i], phi: Eq_angles[1][i], thetaB: Eq_angles[2][i]})
        fkt.append(float(c))
    return fkt


params_Fit = Parameters()
params_Fit.add_many(
    ('K2S', params[0], True, -400000, 400000, None, None),
    ('K2P', params[1], True, -50000, 50000, None, None),
    ('K4S', params[2], False, -400000, 400000, None, None),
    ('K4P', params[3], True, -50000, 50000, None, None),
    ('PHI_U', params[4], True, 0, 2 * m.pi, None, None)
)
model = Model(Model_Fit_Fkt)


def update_rules(params_new):
    # update rules array according to new values inside params_new
    rule_new = [params_new, consts]
    rule_new = [val for sublist in rule_new for val in sublist]
    rules_new = []
    for i in range(len(rule_new)):
        rules_new.append([names[i], rule_new[i]])
    print('Updated rules:', rules_new)
    return rules_new


def iteration(B_Sim_orig, B_Sim2_orig):
    global Eq_angles
    B_Sim = B_Sim_orig
    B_Sim2 = B_Sim2_orig
    it_while = 0
    while np.linalg.norm(B_Sim - B_Sim2) > maxBresDelta:
        print('Error too big! Start iterating')
        if np.array_equal(B_Sim, B_Sim2):
            # check if B_Sim changed
            it_while += 1
            print('Error: Difference between simulated Bres did not change! it_while = ', it_while)
            if it_while == 5:
                break
        B_Sim = np.copy(B_Sim2)
        Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

        ani_result = model.fit(B_Exp, params_Fit, phi_real=phi_RANGE)
        print(ani_result.fit_report())

        temp_paras = []
        for name, param in ani_result.params.items():
            temp_paras.append(float(param.value))
            params_Fit[name].set(value=param.value)
        # temp_paras ist nun array aus Endwerten des Fits: [K2s,K2p,K4s,K4p,PhiU]
        # Aktualisiere nun rules
        for i, l in enumerate(temp_paras):
            params[i] = l

        rules = update_rules(params)
        func = partial(ResFieldNumInp, rules)
        B_Sim2 = pool.map(func, phi_RANGE)
        B_Sim2 = np.asarray(B_Sim2)
        print(np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]))
    return B_Sim2

def create_pre_fit(rules_start):
    pool = Pool()  # Create pool of threads for multiprocessing
    func = partial(ResFieldNumInp,
                   rules_start)  # indirectly add more than one functional arguments without *args or **kwargs
    B_Sim = pool.map(func, phi_RANGE)  # Now this function can be mapped to the pool using an array phi_RANGE
    B_Exp = pool.map(B_inter, phi_RANGE_deg)  # Same as above
    B_Sim = np.asarray(B_Sim)  # convert List to numpy array
    return B_Sim[:,0],B_Exp

def init(fit: bool, rules_start):
    pool = Pool()  # Create pool of threads for multiprocessing
    func = partial(ResFieldNumInp,
                   rules_start)  # indirectly add more than one functional arguments without *args or **kwargs
    B_Sim = pool.map(func, phi_RANGE)  # Now this function can be mapped to the pool using an array phi_RANGE
    B_Exp = pool.map(B_inter, phi_shifted_deg)  # Same as above
    B_Sim = np.asarray(B_Sim)  # convert List to numpy array
    Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

    # Boolean to determine whether you want to fit or just display data
    if fit:
        # Works as follows:
        # This is an implementation of dynamic programming. Solve a numerical problem by splitting it into sub problems.
        # Using iterations, solving of the sub problems will give the solution of the main problem.
        # 0. Evaluate B_res simulation with starting K´s: B_Sim1
        # 1. Fit Model_Fit_FKT to data set in order to get K-Parameters
        #   -Now the EqAngles changed, therefore need to be reevaluated
        # 2. Get new EqAngles using K´s from the Least Square Fit
        # 3. Calculate using new EqAngles and fitted K`s B_Sim2
        # ---------Now the iteration could start--------------
        # 4. calculate error between B_sim1 and B_Sim2, if bigger than maxBresDelta, start iterating
        # 5. Iteration: set B_Sim1 = B_Sim2 and go to step 1. until error is smaller than maxBresDelta

        ani_result = model.fit(B_Exp, params_Fit, phi_real=phi_RANGE)
        print(ani_result.fit_report())

        temp_paras = []
        for name, param in ani_result.params.items():
            temp_paras.append(float(param.value))
            params_Fit[name].set(value=param.value)
        # temp_paras ist nun array aus Endwerten des Fits: [K2s,K2p,K4s,K4p,PhiU]
        # Aktualisiere nun rules
        for i, l in enumerate(temp_paras):
            params[i] = l

        rules = update_rules(params)
        func = partial(ResFieldNumInp, rules)
        B_Sim2 = pool.map(func, phi_RANGE)
        B_Sim2 = np.asarray(B_Sim2)
        print("Magnituden Unterschied von: ", np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]), ", maxBresDelta: ",
              maxBresDelta, " bigger than magnitude? ", np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]) > maxBresDelta)

        B_Sim2 = iteration(B_Sim, B_Sim2)

        plt.plot(phi_RANGE_deg, B_Sim2[:, 0], '-g', label='Simulation2')
        # plt.plot(phi_RANGE_deg, ani_result.best_fit, '-r', label='best fit')
    plt.plot(phi_RANGE_deg, B_Sim[:, 0], '-g', label='Simulation1')
    plt.plot(phi_RANGE_deg, B_Exp, '-b', label='Interpolation')
    plt.scatter(phi_RANGE_deg, B_Exp, label='Experiment')
    plt.xlabel('Angle [Degrees]')
    plt.ylabel('B-Field [Arb. Units]')
    plt.legend()
    plt.show()


'''if __name__ == '__main__':
    init(False, rules)'''

# Ohne Global Minimum search 8,6 sek.
# Shgo: 19.2 sec.
