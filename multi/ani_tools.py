from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from symengine.lib.symengine_wrapper import solve
from sympy import sin, cos, sqrt, Eq, Symbol, Function, Matrix, re, sympify
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import shgo, differential_evolution  # activate when using global optimization
#from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Pool
from functools import partial
import numpy as np
import matplotlib as mp

mp.use('QT5Agg')
import matplotlib.pyplot as plt
import math as m
import time


class Worker(object):
    """docstring for Worker"""
    def __init__(self, arg):
        super(Worker, self).__init__()
        self.arg = arg
        
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

omega = Symbol('omega')
g = Symbol('g')
M = Symbol('M')
K4s = Symbol('K4s')

# ------const Params(always present)--------
mu0 = Symbol('mu0')
mub = Symbol('mub')
B = Symbol('B')
phi = Symbol('phi')
theta = Symbol('theta')
thetaB = Symbol('thetaB')
phiB = Symbol('phiB')
gamma = Symbol('gamma')
# w = Symbol('w')
hbar = Symbol('hbar')
f = Symbol('f')
F = Function('F')


def init_load(filename, FreeE, fit_params, fixed_params, shifts, anglestep, Fit: bool, Plot: bool, *args):
    #global B_RES
    # filename string,  value, path to parameter.dat file used for fitting
    # FreeE string,  of Free Energy Density
    # fit_params dict => Parameter + start value
    # fixed_params dict => Fixed parameters + value
    # shift: float
    # anglestep: integer
    # Fit: bool, to control whether to fit or not
    # Plot: bool, same as Fit

    # The rest will be passed as *arg
    # phi_array: array, of angles Phi of phi in degrees shifted!
    # phi_array_not_shifted: orig array 
    # R_raw: array of raw Bres fitted, taken from a file

    if args:
        phi_array = args[0]
        print(filename, FreeE, fit_params, fixed_params, shifts, anglestep, phi_array)
    else:
        print(filename, FreeE, fit_params, fixed_params, shifts, anglestep,Fit,Plot)
    # ------fixed Params--------
    omega = Symbol('omega')
    g = Symbol('g')
    M = Symbol('M')
    K4s = Symbol('K4s') # -----------------------------------------------------------------change that, K4s does not belong to fixed params (hardcoded!!)!!!!!!-------------------------------------------------------

    # ------const Params(always present)--------
    mu0 = Symbol('mu0')
    mub = Symbol('mub')
    B = Symbol('B')
    phi = Symbol('phi')
    theta = Symbol('theta')
    thetaB = Symbol('thetaB')
    phiB = Symbol('phiB')
    gamma = Symbol('gamma')
    # w = Symbol('w')
    hbar = Symbol('hbar')
    f = Symbol('f')
    F = Function('F')

    maxBresDelta = 0.01  # Also adjustable by gui?

    shift = -shifts

    # Add consts to fixed params dict
    fixed_params[mu0] = 4 * m.pi * 10 ** (-7)
    fixed_params[mub] = 9.27 * 10 ** (-24)
    fixed_params[hbar] = (6.62606957 * 10 ** (-34)) / (2 * m.pi)
    #gamma = g*mub/hbar 
    fixed_params[gamma] = fixed_params['g'] * fixed_params[mub] / fixed_params[hbar]
    #       ------------------------------------------gamma missing? -----------------------------------------------

    # Laden der gefitteten Resonanzpositionen
    D = np.loadtxt(filename, dtype='float', skiprows=0)
    R_raw = D[:, 3]
    Winkel = D[:, 5].flatten()
    B_inter = interp1d(Winkel, R_raw)  # Interpoliertes B_res, array mit Länge len(Winkel)

    #FreeE = 'B*M*(sin(theta)*sin(thetaB)*cos(phi - phiB) + cos(theta)*cos(thetaB)) - K2p*sin(theta)**2*cos(phi - phiu)**2 - K4p*(cos(4*phi) + 3)*sin(theta)**4/8 - K4s*cos(theta)**4/2 - (-K2s + M**2*mu0/2)*sin(theta)**2'
    # Freie Energy Formel + Baselgia
    F = sympify(FreeE)  # let sympy interpret the FreeE string as a function

    # Use Baselgia approach for Bres (Has no singularity at theta = 0!)
    halb_res = gamma ** 2 / M ** 2 * (
            F.diff(theta, 2) * (F.diff(phi, 2) / (sin(theta)) ** 2 + cos(theta) / sin(theta) * F.diff(theta)) - (
            F.diff(theta, phi) / sin(theta) - cos(theta) / sin(theta) * F.diff(phi) / sin(theta)) ** 2)
    eq_halb_res = Eq(omega ** 2, halb_res)

    # Formel der Resonanzfrequenz B_RES, Lösung [1] ist positive Lösung der Wurzelfunktion (glaub ich)
    B_RES = solve(eq_halb_res, B)  # Solves the Baselgia approach for B
    B_RES = B_RES.args[1]  # Choose second solution of solve!

    # Now create rule as dict(). rule is then variable that contains every parameter information (params+const)
    rule = {**fit_params, **fixed_params}
    # ----->  [K2s, K2p, K4s, K4p, phiU, f, g, M, mu0, muB, hbar, gamma, omega]

    # Create angle arrays
    phi_min = min(Winkel) * m.pi / 180  # define smallest value
    phi_max = max(Winkel) * m.pi / 180  # define biggest value

    Winkel_min = min(Winkel)
    Winkel_max = max(Winkel)

    phi_step = anglestep  # define stepwidth (resolution)
    phi_RANGE = np.arange(phi_min, phi_max, phi_step, dtype='float')  # array of radians, with stepwidth = anglestep
    phi_RANGE_deg = np.multiply(phi_RANGE, 180 / m.pi)  # convert to degrees

    # phi_array # shifted array from GUI (make_phirange)

    # Create Paramters dict() for the leastsq fit
    

    if Fit:
        #plot: bool, rules_start: dict, phi_RANGE: list, phi_array: list, B_inter: func, B_RES: sympy_object, F: sympy_object
        phi_array = make_phirange(shift, phi_RANGE_deg, True, Winkel_min, Winkel_max)

        func_args = [i for i in fit_params.keys()] # get the names/keys of the fitted params
        func_args = str(func_args).replace('[','').replace(']','').replace("'",'') # prepare the list for the string function
        func_str = create_str_func(func_args)
        exec(func_str,globals())  # This will create the function give by func_str: "Model_Fit_Fkt"

        # Then create Parameter dict() and model for lmfit
        params_Fit = Parameters()
        for name, value in fit_params.items():
            params_Fit.add(name, value)
        model = Model(Model_Fit_Fkt)

        main_loop(Plot, rule, phi_RANGE, phi_array, B_inter, B_RES, F, model, params_Fit, fixed_params, fit_params, maxBresDelta)
    else:
        result = create_pre_fit(rule, phi_RANGE, phi_RANGE_deg,B_inter, B_RES, F)
        return result

def create_str_func(func_args):
    # this function creates a string, that in the later part will be conveted to an expression
    func = """
def Model_Fit_Fkt(phi_real, .ß123, **kwargs):
    
    
    fkt = []
    B_RES = kwargs['B_RES']
    Eq_angles = kwargs['Eq_angles']
    #print("Now it comes!:")
    #print(locals())
    #print(kwargs)

    Local_dict = locals()
    Local_dict.pop('kwargs')
    Local_dict.pop('B_RES')
    Local_dict.pop('phi_real')
    Local_dict.pop('fkt')
    Local_dict.pop('Eq_angles')
    #print(Local_dict)

    a = B_RES.subs({i: l for i, l in kwargs['fixed_params'].items()})  # Parameter, B and EqAngle missing

    for i, l in enumerate(phi_real):
        b = a.subs({i: l for i, l in Local_dict.items()})
        c = b.subs({phiB: l, theta: Eq_angles[0][i], phi: Eq_angles[1][i], thetaB: Eq_angles[2][i]})
        fkt.append(float(c))
    return fkt""".replace('.ß123',func_args)
    return func



'''def Model_Fit_Fkt(phi_real, K2s, K2p, K4p, phi_u):
    # is the function used for least square fitting
    # inputs are values for parameter
    # outputs array "fkt", an element in fkt is the value of the function at point x=phi_real, so f(x)
    fkt = []
    print('\n\n\n\n',K2s, K2p, K4p, phi_u)
    a = B_RES.subs({i: l for i, l in rule_consts})  # Parameter, B and EqAngle missing
    for i, l in enumerate(phi_real):
        b = a.subs({phiB: l, K2s: K2S, K2p: K2P, K4s: K4S, K4p: K4P, phiU: PHI_U})
        c = b.subs({theta: Eq_angles[0][i], phi: Eq_angles[1][i], thetaB: Eq_angles[2][i]})
        fkt.append(float(c))
    return fkt'''

def iteration(B_Sim_orig, B_Sim2_orig, B_RES, fixed_params, B_inter, F, maxBresDelta, model, params_Fit, B_Exp, phi_RANGE, fit_params, pool):
    B_Sim = B_Sim_orig
    B_Sim2 = B_Sim2_orig
    it_while = 0
    while np.linalg.norm(B_Sim - B_Sim2) > maxBresDelta:
        print('Error too big! Start iterating')

        if np.array_equal(B_Sim, B_Sim2):  # Check if both arrays are equal (if there is a change)
            it_while += 1  # if it hasn't changed add 1 to it_while
            print('Error: Difference between simulated Bres did not change! it_while = ', it_while)
            if it_while == 5:
                break  # stop iterating if B_Sim and B_Sim2 are 5 times equal

        B_Sim = np.copy(B_Sim2)  # To start with the result from previous iteration
        Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

        ani_result = model.fit(B_Exp, params_Fit, phi_real=phi_RANGE, B_RES=B_RES, fixed_params=fixed_params, Eq_angles=Eq_angles)  # Fit
        print(ani_result.fit_report())

        # Refresh the parameter dict()
        for name, param in ani_result.params.items():
            fit_params[name] = param.value
            params_Fit[name].set(value=param.value)

        # Refresh rule
        rule = update_rules(fit_params, fixed_params)

        B_FKT = str(B_RES.subs({i: l for i, l in rule.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
        FKT = F.subs({i: l for i, l in rule.items()})  # rule beeing inserted, for minimizing 

        # Recalculate EqAngles
        func = partial(ResFieldNumInp, B_inter, B_FKT, FKT)
        B_Sim2 = pool.map(func, phi_RANGE)
        B_Sim2 = np.asarray(B_Sim2)

        # print error between old and new result
        print(np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]))
    return B_Sim2,params_Fit

def ResFieldNumInp(B_inter, B_RES, F, phi_val):
    #rule, B_inter, B_RES, F
    phiB_val = phi_val
    FKT = F

    # Here Python with multiprocessing behaves weird. If one would give ResFieldNumInp B_RES without subs, the interpreter can no longe pickle the data
    # (means, he somehow cant convert it to binary data and back, could be a bug in symengine!)
    # once B_RES was subs and converted to string it somehow works. One disadvantage of this mehtod is, 
    # that the string B_RES now has to be interpreted each time the function is called
    # i dont now if this has any performance affections (it may even be faster; nothing noticeable yet).
    B_FKT = sympify(B_RES)

    # takes EqAngles and solves B_RES using interpolated datal
    
    B = B_inter(phi_val * 180 / m.pi)
    #B_FKT = B_RES.subs({i: l for i, l in rule.items()})  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    #FKT = F.subs({i: l for i, l in rule.items()})  # rule beeing inserted, for minimizing
    Eq_angles = solveAngles(B, phiB_val, FKT)
    result = B_FKT.subs({'theta': Eq_angles[0], 'phi': Eq_angles[1], 'thetaB': Eq_angles[2], 'phiB': phiB_val})
    #print(float(result))
    return [float(result), float(Eq_angles[0]), float(Eq_angles[1]), float(Eq_angles[2])]

def F_fkt(x, *args):
    # Function of Free Energy Density, that should be minimized
    # *args is array:[B,phiB]
    # x are the fitted params
    Fit_fkt = args[2].subs({'B': args[0], 'phiB': args[1], 'theta': x[0], 'phi': x[1], 'thetaB': x[2]})
    return float(Fit_fkt)

def update_rules(params_new, fixed_params):
    # update rules array according to new values inside params_new
    rule_new = {**params_new, **fixed_params}
    print('Updated rule:', rule_new)
    return rule_new

def solveAngles(B, phiB, fkt):
    # Minimize Angles Theta, Phi, ThetaB
    # start_paras for the moment these are constant: [m.pi/2, m.pi, m.pi/2]
    # start_paras = [m.pi/2, m.pi, m.pi/2]
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
    Bounds = [(0, 2 * m.pi), (0, 2 * m.pi), (0, 2 * m.pi)]
    start_paras = np.asarray([m.pi / 2, m.pi, m.pi / 2])
    result = minimize(F_fkt, start_paras, args=(B, phiB, fkt), method='L-BFGS-B',
                      bounds=Bounds)  # alternatively a global search can be performed using shgo or differential_evolution, look up scipy Docs.
    # result = shgo(F_fkt,Bounds, args=(B,phiB,fkt))  #Works as fast as local search (maybe even faster), but produces weird artifacts in plot!?
    # result = dual_annealing(F_fkt,Bounds, args=(B,phiB,fkt))   # Finds nice solution but takes years to compute one array with Anglestep = m.pi/40
    # result = differential_evolution(F_fkt,Bounds, args=(B,phiB,fkt)) #slower than shgo, but also nice solution
    ##print(result.x[0],result.x[1],result.x[2],"\"Local Minimum\"")
    return result.x[0], result.x[1], result.x[2]

def create_pre_fit(rules_start, phi_RANGE, phi_RANGE_deg, B_inter, B_RES, F):
    pool = Pool()  # Create pool of threads for multiprocessing

    B_FKT = str(B_RES.subs({i: l for i, l in rules_start.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rules_start.items()})  # rule beeing inserted, for minimizing                

    func = partial(ResFieldNumInp, B_inter, B_FKT, FKT) # indirectly add more than one functional arguments without *args or **kwargs
    B_Sim = pool.map(func, phi_RANGE)  # Now this function can be mapped to the pool using an array phi_RANGE
    B_Exp = pool.map(B_inter, phi_RANGE_deg)  # Same as above
    B_Sim = np.asarray(B_Sim)  # convert List to numpy array
    return B_Sim[:, 0], B_Exp, phi_RANGE_deg

def main_loop(plot: bool, rules_start, phi_RANGE, phi_array, B_inter, B_RES, F, model, params_Fit, fixed_params, fit_params,maxBresDelta):
    #phi_RANGE is array of rad not shifted
    #phi_array is aray of deg shifted

    pool = Pool()  # Create pool of threads for multiprocessing

    B_FKT = str(B_RES.subs({i: l for i, l in rules_start.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rules_start.items()})  # rule beeing inserted, for minimizing  

    func = partial(ResFieldNumInp, B_inter, B_FKT, FKT)  # indirectly add more than one functional arguments without *args or **kwargs
    B_Sim = pool.map(func, phi_RANGE)  # Now this function can be mapped to the pool using an array phi_RANGE
    B_Exp = pool.map(B_inter, phi_array)  # Same as above
    B_Sim = np.asarray(B_Sim)  # convert List to numpy array
    Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

    # Works as follows:
    # This is an implementation of dynamic programming. Solve a numerical problem by splitting it into sub problems (i think).
    # Using iterations, solving of the sub problems will give the solution of the main problem.
    # 0. Evaluate B_res simulation with starting K´s: B_Sim1
    # 1. Fit Model_Fit_FKT to data set in order to get K-Parameters
    #   -Now the EqAngles changed, therefore need to be reevaluated
    # 2. Get new EqAngles using K´s from the Least Square Fit
    # 3. Calculate using new EqAngles and fitted K`s B_Sim2
    # ---------Now the iteration could start--------------
    # 4. calculate error between B_sim1 and B_Sim2, if bigger than maxBresDelta, start iterating
    # 5. Iteration: set B_Sim1 = B_Sim2 and go to step 1. until error is smaller than maxBresDelta

    ani_result = model.fit(B_Exp, params_Fit, phi_real=phi_RANGE, B_RES=B_RES, fixed_params=fixed_params, Eq_angles=Eq_angles)
    print(ani_result.fit_report())

    # Refresh the parameter dict()
    for name, param in ani_result.params.items():
        fit_params[name] = param.value
        params_Fit[name].set(value=param.value)

    # Refresh rule
    rule = update_rules(fit_params, fixed_params)

    B_FKT = str(B_RES.subs({i: l for i, l in rule.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rule.items()})  # rule beeing inserted, for minimizing 

    func = partial(ResFieldNumInp, B_inter, B_FKT, FKT)
    B_Sim2 = pool.map(func, phi_RANGE)
    B_Sim2 = np.asarray(B_Sim2)

    delta = np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]) 
    print("Magnituden Unterschied von: ", delta, ", maxBresDelta: ",
          maxBresDelta, " bigger than magnitude? ", delta > maxBresDelta)

    if delta > maxBresDelta:
        # if delta is too big start iterating
        Sim_result = iteration(B_Sim, B_Sim2, B_RES, fixed_params,B_inter, F,maxBresDelta, model, params_Fit, B_Exp, phi_RANGE, fit_params, pool)
        # Sim_result is now a list 5D List: 
        #    [:, 0] = B_Sim
        #    [:, 1] = Eq_Angle1
        #    [:, 2] = Eq_Angle2
        #    [:, 3] = Eq_Angle3
        #    [:, 5] = best fit params dict
        B_Sim2 = Sim_result[0]
        print(Sim_result[1])

    if plot:
        phi_deg = pool.map(rad_to_deg, phi_RANGE)
        plt.cla()
        plt.clf()
        plt.plot(phi_deg, B_Sim2[:, 0], '-g', label='Simulation')
        # plt.plot(phi_RANGE_deg, ani_result.best_fit, '-r', label='best fit')
        #plt.plot(phi_deg, B_Sim[:, 0], '-g', label='Simulation1')
        plt.plot(phi_deg, B_Exp, '-b', label='Interpolation')
        plt.scatter(phi_deg, B_Exp, label='Experiment')
        plt.xlabel('Angle [Degrees]')
        plt.ylabel('B-Field [Arb. Units]')
        plt.legend()
        plt.show()

def rad_to_deg(val):
    val = val*180/m.pi
    return val

def make_phirange(shift_value: float, phi_array: list, deg: bool, minWinkel: float, maxWinkel: float):
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
        #Not used at the moment
        shift_value = shift_value * m.pi / 180
        for l, i in enumerate(phi_array):
            phi_array[l] += shift_value
        for i in range(len(phi_array)):
            if phi_array[i] > phi_max:
                phi_array[i] -= phi_max
            elif phi_array[i] < phi_min:
                phi_array[i] += phi_max
    return phi_array




if __name__ == '__main__':

    #For debug only
    
    #text = "C:/Users/Jonas/Desktop/test.dat B*M*(sin(theta)*sin(thetaB)*cos(phi - phiB) + cos(theta)*cos(thetaB)) - K2p*sin(theta)**2*cos(phi - phiU)**2 - K4p*(cos(4*phi) + 3)*sin(theta)**4/8 - K4s*cos(theta)**4/2 - (-K2s + M**2*mu0/2)*sin(theta)**2 {'K2p': 863.25, 'K2s': 261345.0, 'K4p': 13720.6, 'phiu': 5.0756} {'omega': 62066561101.381386, 'g': 2.05, 'M': 1530000.0, 'K4s': 0} 0.0 0.03490658503988659 False False"
    init_load('C:/Users/Jonas/Desktop/test.dat','B*M*(sin(theta)*sin(thetaB)*cos(phi - phiB) + cos(theta)*cos(thetaB)) - K2p*sin(theta)**2*cos(phi - phiU)**2 - K4p*(cos(4*phi) + 3)*sin(theta)**4/8 - K4s*cos(theta)**4/2 - (-K2s + M**2*mu0/2)*sin(theta)**2', {'K2p': 863.25, 'K2s': 100000.0, 'K4p': 13720.6, 'phiU': 5.0756}, {'omega': 62066561101.381386, 'g': 2.05, 'M': 1530000.0, 'K4s': 0}, 26.0, 0.03490658503988659,True, True,[332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356,   0,
   2,   4,   6,   8,  10,  12,  14,  16,  18,  20,  22,  24,  26,  28,
  30,  32,  34,  36,  38,  40,  42,  44,  46,  48,  50,  52,  54,  56,
  58,  60,  62,  64,  66,  68,  70,  72,  74,  76,  78,  80,  82,  84,
  86,  88,  90,  92,  94,  96,  98, 100, 102, 104, 106, 108, 110, 112,
 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140,
 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168,
 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196,
 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224,
 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252,
 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, 278, 280,
 282, 284, 286, 288, 290, 292, 294, 296, 298, 300, 302, 304, 306, 308,
 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332])

# Ohne Global Minimum search 8,6 sek.
# Shgo: 19.2 sec.
