from lmfit import Model
from lmfit import Parameters
from lmfit import Parameter
from symengine.lib.symengine_wrapper import solve
from symengine import sin, cos, sqrt, sympify, Symbol, Function, Eq
#from sympy import sin, cos, sqrt, Eq, Symbol, Function, Matrix, re, sympify
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import shgo, differential_evolution, basinhopping  # activate when using global optimization
from pathos.multiprocessing import ProcessingPool as Pool

'''try:
    from ray.util.multiprocessing import Pool
except:
    from multiprocessing import Pool'''
#from multiprocessing import Pool
from functools import partial
import numpy as np
import matplotlib as mp

mp.use('QT5Agg')
import matplotlib.pyplot as plt
import math as m
import os
import datetime
import time

# Create Symbols for symengine to symbolically solve equations

omega = Symbol('omega')
g = Symbol('g')
M = Symbol('M')

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


def get_fit_options_from_file(fname):
    with open(fname) as f:
        header_parts = []
        Title = f.readline()

        if Title.split('# ')[-1] == 'FMR-Fit\n':
            header_parts.append(Title)
            header_parts.append(f.readline())
        # print(header_parts[2])

        Lineshape = int(header_parts[1].split('$.')[1].split(' Lineshape')[0])  # Extract lineshape out of params file
        # If Lineshape = 2 its Lorentz
        # If Lineshape = 3 its Dyson
        try:
            fit_num = int(header_parts[1].split(' ')[11])
        except:
            print('Could not extract numbers of fit functions. Assuming one function!')
            fit_num = 1
    return Lineshape, fit_num

def init_load(filename, FreeE, fit_params, fixed_params, shifts, anglestep, Fit: bool, Plot: bool, ani_ori: bool, fit_select:int, *args):
    #global B_RES
    # filename string,  value, path to parameter.dat file used for fitting
    # FreeE string,  of Free Energy Density
    # fit_params dict => Parameter + start value
    # fixed_params dict => Fixed parameters + value
    # shift: float
    # anglestep: integer
    # Fit: bool, to control whether to fit or not
    # Plot: bool, same as Fit

    print(filename, FreeE, fit_params, fixed_params, shifts, anglestep, Fit, Plot)

    maxBresDelta = 0.001  # Also adjustable by gui?

    shift = shifts

    # Add consts to fixed params dict
    fixed_params[mu0] = 4 * m.pi * 10 ** (-7)
    fixed_params[mub] = 9.27 * 10 ** (-24)
    fixed_params[hbar] = (6.62606957 * 10 ** (-34)) / (2 * m.pi)
    #gamma = g*mub/hbar 
    fixed_params[gamma] = fixed_params['g'] * fixed_params[mub] / fixed_params[hbar]

    # Load fitted Params and determine Lineshape

    Lineshape, fit_num = get_fit_options_from_file(filename)

    D = np.loadtxt(filename, dtype='float', skiprows=0) # Data
    if Lineshape == 2: #Lorentz
        R_raw = D[:, 3 * fit_select]
        Winkel = D[:, 3 * fit_num + 2].flatten()
    elif Lineshape == 3: #Dyson
        print("Dyson")
        R_raw = D[:, 4 * fit_select]
        Winkel = D[:, 4 * fit_num + 2].flatten()
    else: # Rest (not implemented yet)
        print("Please use either Lorentz or Dyson shape!")
    
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

    Winkel_min = min(Winkel) # its in deg; Orignial data
    Winkel_max = max(Winkel) # its in deg; Orignial data
    angle_step = anglestep  # define stepwidth (resolution)
    angle_RANGE_deg = np.arange(Winkel_min,Winkel_max,angle_step*180/m.pi)  # Phi Array in degress NOT shifted
    reference_data = B_inter(angle_RANGE_deg) # Interpolate Data from experiment with array steps set by programm
    # Now the idea is, we have an array (B_inter_data_array) of real world data, which is shifted in +- direction,
    # we therefore need to shift this.
    # This reference_data is not shifted, the shift is introduced in angle_RANGE. The lmfit routine then uses
    # the information of the shifted angles to map/fit the interpolated data.

    #try:

    def correct_angle_RANGE(angle_RANGE,len_ref, count=0, add='-', lauf_var=0, shift=0):
        if add == '-':
            angle_RANGE = np.arange(angle_min + shift, angle_max + shift, m.pi / (1 / (angle_step / m.pi) - lauf_var), dtype='float')
        else:
            angle_RANGE = np.arange(angle_min + shift, angle_max + shift, m.pi / (1 / (angle_step / m.pi) + lauf_var), dtype='float')
        if len(angle_RANGE) != len_ref:
            print(len(angle_RANGE),len_ref,count)
            if count < 10:
                correct_angle_RANGE(angle_RANGE, len_ref, count+1,'-', lauf_var+1)
            elif count < 20:
                if count == 10:
                    lauf_var = 0
                correct_angle_RANGE(angle_RANGE,len_ref,count+1,'+', lauf_var+1)
            else:
                print('angle_RANGE and reference_data have different length, change Anglestep (~ +-1 ) and try again')
        else:
            return angle_RANGE

    is_OOP = ani_ori

    if Fit:
        angle_min = (Winkel_min + shift) * m.pi / 180  # define smallest value; shifted
        angle_max = (Winkel_max + shift) * m.pi / 180  # define biggest value; shifted
       # angle_RANGE = np.arange(angle_min, angle_max, angle_step, dtype='float')  # array of radians, with stepwidth = anglestep
        angle_RANGE = np.arange(Winkel_min + shift, Winkel_max + shift, angle_step*180/m.pi) * m.pi/180
        # There is the possibility, that the length of angle_RANGE can be smaller or bigger than reference_data
        if len(angle_RANGE) != len(reference_data):
            angle_RANGE = correct_angle_RANGE(angle_RANGE, len(reference_data),0 , '-', 0, shift*m.pi/180)

        # Then limit the array so that it doesnt exceed 0-360 degrees
        for i, l in enumerate(angle_RANGE):
            if l >= 359.99 * m.pi / 180:
                angle_RANGE[i] -= 360.0 * m.pi / 180
            elif l <= 0.0:
                angle_RANGE[i] += 360.0 * m.pi / 180
        print('Simulating Bres and fitting anisotropy constants. Your computer may be unresponsive')
        end_pfad = init_folder(filename)


        func_args = [i for i in fit_params.keys()] # get the names/keys of the fitted params
        func_args = str(func_args).replace('[','').replace(']','').replace("'",'') # prepare the list for the string function
        func_str = create_str_func(func_args, is_OOP)
        exec(func_str,globals())  # This will create the function given by func_str: "Model_Fit_Fkt"

        # Then create Parameter dict() and model for lmfit
        params_Fit = Parameters()
        for name, value in fit_params.items():
            params_Fit.add(name, value)
        model = Model(Model_Fit_Fkt)

        # plot: bool, rules_start: dict, angle_RANGE: list, phi_array: list, B_inter: func, B_RES: sympy_object, F: sympy_object
        main_loop(Plot, rule, angle_RANGE, reference_data, B_RES, F, model, params_Fit, fixed_params, fit_params, maxBresDelta, end_pfad, is_OOP)
    else:   # Pre Fit
        angle_RANGE = angle_RANGE_deg * m.pi/180
        if len(angle_RANGE) != len(reference_data):
            angle_RANGE = correct_angle_RANGE(angle_RANGE, len(reference_data))

        print('Creating pre fit')
        result = create_pre_fit(rule, angle_RANGE, angle_RANGE_deg, reference_data, B_RES, F, is_OOP)
        return result
    #except Exception as e:
    #    print('Error in ani_tools.init_load(): ',e)


def init_folder(filename):
    # This function should get the filename string from the main script and then checks if folders are present, if not they will be created
    #pfad is the path to the folder the user saved his parameters to
    pfad = os.path.dirname(filename)
    dat_name = filename.replace(pfad + '/', '').replace('.dat','')
    dir_path = pfad + '/{}_Anisotropy'.format(dat_name)

    if not os.path.exists(dir_path):
        print('Folder Anisotropy does not exist')
        print('Creating /Anisotropy ....')
        os.mkdir(dir_path)
    return dir_path + '/'

def create_str_func(func_args, is_OOP):
    # this function creates a string, that in the later part will be converted to an expression
    # It is necessary for taking varying K parameters into account. With this it is possible to create a function with
    # varying explicit function parameters needed for the Fit
    if is_OOP: # OOP
        func = """
def Model_Fit_Fkt(angle_real, .ß123, **kwargs):


    fkt = []
    B_RES = kwargs['B_RES']
    Eq_angles = kwargs['Eq_angles']
    #print("Now it comes!:")
    #print(locals())
    #print(kwargs)

    Local_dict = locals()
    Local_dict.pop('kwargs')
    Local_dict.pop('B_RES')
    Local_dict.pop('angle_real')
    Local_dict.pop('fkt')
    Local_dict.pop('Eq_angles')
    #print(Local_dict)

    a = B_RES.subs({i: l for i, l in kwargs['fixed_params'].items()})  # Parameter, B and EqAngle missing

    for i, l in enumerate(angle_real):
        b = a.subs({i: l for i, l in Local_dict.items()})
        c = b.subs({thetaB: l, theta: Eq_angles[0][i], phi: Eq_angles[1][i], phiB: Eq_angles[2][i]})
        fkt.append(float(c))
    return fkt""".replace('.ß123', func_args)
    else: # IP
        func = """
def Model_Fit_Fkt(angle_real, .ß123, **kwargs):
    
    
    fkt = []
    B_RES = kwargs['B_RES']
    Eq_angles = kwargs['Eq_angles']
    #print("Now it comes!:")
    #print(locals())
    #print(kwargs)

    Local_dict = locals()
    Local_dict.pop('kwargs')
    Local_dict.pop('B_RES')
    Local_dict.pop('angle_real')
    Local_dict.pop('fkt')
    Local_dict.pop('Eq_angles')
    #print(Local_dict)

    a = B_RES.subs({i: l for i, l in kwargs['fixed_params'].items()})  # Parameter, B and EqAngle missing

    for i, l in enumerate(angle_real):
        b = a.subs({i: l for i, l in Local_dict.items()})
        c = b.subs({phiB: l, theta: Eq_angles[0][i], phi: Eq_angles[1][i], thetaB: Eq_angles[2][i]})
        fkt.append(float(c))
    return fkt""".replace('.ß123',func_args)
    return func

def iteration(B_Sim_orig, B_Sim2_orig, B_RES, fixed_params, reference_data, F, maxBresDelta, model, params_Fit, angle_RANGE, fit_params, pool, is_OOP):
    B_Sim = B_Sim_orig
    B_Sim2 = B_Sim2_orig
    it_while = 0
    while np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]) > maxBresDelta:

        print('Error too big! Start iterating; Iteration: ',it_while)
        B_Sim = np.copy(B_Sim2)  # To start with the result from previous iteration
        Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

        ani_result = model.fit(reference_data, params_Fit, angle_real=angle_RANGE, B_RES=B_RES, fixed_params=fixed_params, Eq_angles=Eq_angles)  # Fit
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
        func = partial(ResFieldNumInp, B_FKT, FKT, is_OOP)
        #B_Sim2 = pool.map(func, angle_RANGE, reference_data)
        B_Sim2 = []
        for index, val in enumerate(angle_RANGE):
            B_Sim2.append(func(val, reference_data[index]))
        B_Sim2 = np.asarray(B_Sim2)

        # print error between old and new result
        print(np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]))

        if it_while > 5:   # End statement if algorithm is iterating too often
            print("Stop Iteration, Iter_Count is too high")
            break
        else:
            it_while += 1

    return B_Sim2,fit_params

def ResFieldNumInp(B_RES, F, is_OOP, phi_val, reference_data):
    phiB_val = phi_val
    FKT = F

    # Here Python with multiprocessing behaves weird. If one would give ResFieldNumInp B_RES without subs, the interpreter can no longe pickle the data
    # (means, he somehow cant convert it to binary data and back, could be a bug in symengine!)
    # once B_RES was subs and converted to string it somehow works. One disadvantage of this mehtod is, 
    # that the string B_RES now has to be interpreted each time the function is called
    # i dont now if this has any performance affections (it may even be faster; nothing noticeable yet).
    B_FKT = sympify(B_RES)

    # takes EqAngles and solves B_RES using interpolated datal
    B = reference_data
    #B_FKT = B_RES.subs({i: l for i, l in rule.items()})  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    #FKT = F.subs({i: l for i, l in rule.items()})  # rule beeing inserted, for minimizing
    Eq_angles = solveAngles(B, phiB_val, FKT, is_OOP)
    if is_OOP: # OOP
        result = B_FKT.subs({'theta': Eq_angles[0], 'phi': Eq_angles[1], 'phiB': Eq_angles[2], 'thetaB': phiB_val})
    else: # IP
        result = B_FKT.subs({'theta': Eq_angles[0], 'phi': Eq_angles[1], 'thetaB': Eq_angles[2], 'phiB': phiB_val})
    #print(float(result))
    return [float(result), float(Eq_angles[0]), float(Eq_angles[1]), float(Eq_angles[2])]

def F_fkt(x, *args):
    # Function of Free Energy Density, that should be minimized
    # *args is array:[B,phiB,fkt,is_OOP] or [B,thetaB,fkt,is_OOP]
    # x are the fitted params
    if args[3]: # OOP
        Fit_fkt = args[2].subs({'B': args[0], 'thetaB': args[1], 'theta': x[0], 'phi': x[1], 'phiB': x[2]})
    else:   # IP
        Fit_fkt = args[2].subs({'B': args[0], 'phiB': args[1], 'theta': x[0], 'phi': x[1], 'thetaB': x[2]})
    return float(Fit_fkt)

def update_rules(params_new, fixed_params):
    # update rules array according to new values inside params_new
    rule_new = {**params_new, **fixed_params}
    print('Updated rule:', rule_new)
    return rule_new

def solveAngles(B, phiB, fkt, is_OOP):
    # Todo: Add Robust Fit Method as Global Solver
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
    if is_OOP: # OOP
        #[Theta, Phi, PhiB]
        start_paras = np.asarray([m.pi / 2, m.pi, m.pi])
    else:
        # [Theta, Phi, ThetaB]
        start_paras = np.asarray([m.pi / 2, m.pi, m.pi / 2])
    result = minimize(F_fkt, start_paras, args=(B, phiB, fkt, is_OOP), method='L-BFGS-B',bounds=Bounds)
    # alternatively a global search can be performed using shgo or differential_evolution, look up scipy Docs.
    # result = shgo(F_fkt,Bounds, args=(B,phiB,fkt))  #Works as fast as local search (maybe even faster), but produces weird artifacts in plot!?
    # result = dual_annealing(F_fkt,Bounds, args=(B,phiB,fkt))   # Finds nice solution but takes years to compute one array with Anglestep = m.pi/40
    #result = differential_evolution(F_fkt,Bounds, args=(B,phiB,fkt)) #slower than shgo, but also nice solution
    ##print(result.x[0],result.x[1],result.x[2],"\"Local Minimum\"")

    # Basinhopping might be a good approach for global minima search!
    #args = (B, phiB, fkt, is_OOP)
    #result = basinhopping(F_fkt,start_paras, T=0.8, stepsize=m.pi/8, niter=5,minimizer_kwargs={'args':args,"method":"L-BFGS-B"})
    return result.x[0], result.x[1], result.x[2]

def create_pre_fit(rules_start, angle_RANGE, angle_RANGE_deg, reference_data, B_RES, F, is_OOP):
    #pool = Pool()  # Create pool of threads for multiprocessing
    B_FKT = str(B_RES.subs({i: l for i, l in rules_start.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rules_start.items()})  # rule beeing inserted, for minimizing
    func = partial(ResFieldNumInp, B_FKT, FKT, is_OOP) # indirectly add more than one functional arguments without *args or **kwargs
    #B_Sim = pool.map(func, angle_RANGE, reference_data)  # Now this function can be mapped to the pool using an array angle_RANGE
    B_Sim = []
    for index, val in enumerate(angle_RANGE):
        B_Sim.append(func(val,reference_data[index]))
    B_Sim = np.asarray(B_Sim)  # convert List to numpy array
    #pool.close()
    return B_Sim[:, 0], reference_data, angle_RANGE_deg

def main_loop(plot: bool, rules_start, angle_RANGE, reference_data, B_RES, F, model, params_Fit, fixed_params, fit_params, maxBresDelta, end_pfad, is_OOP):
    #angle_RANGE is array of rad not shifted
    #phi_array is aray of deg shifted
    start = time.time()
    
    #pool = Pool()  # Create pool of threads for multiprocessing

    B_FKT = str(B_RES.subs({i: l for i, l in rules_start.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rules_start.items()})  # rule beeing inserted, for minimizing

    func = partial(ResFieldNumInp, B_FKT, FKT, is_OOP)  # indirectly add more than one functional arguments without *args or **kwargs
    #B_Sim = pool.map(func, angle_RANGE, reference_data)  # Now this function can be mapped to the pool using an array angle_RANGE
    B_Sim = []
    for index, val in enumerate(angle_RANGE):
        B_Sim.append(func(val, reference_data[index]))
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
    ani_result = model.fit(reference_data, params_Fit, angle_real=angle_RANGE, B_RES=B_RES, fixed_params=fixed_params, Eq_angles=Eq_angles)
    print(ani_result.fit_report())

    # Refresh the parameter dict()
    for name, param in ani_result.params.items():
        fit_params[name] = param.value
        params_Fit[name].set(value=param.value)

    # Refresh rule
    rule = update_rules(fit_params, fixed_params)

    B_FKT = str(B_RES.subs({i: l for i, l in rule.items()}))  # rule beeing inserted, missing : B, theta, thetaB, phi, phiB
    FKT = F.subs({i: l for i, l in rule.items()})  # rule beeing inserted, for minimizing 

    func = partial(ResFieldNumInp, B_FKT, FKT, is_OOP)
   #B_Sim2 = pool.map(func, angle_RANGE, reference_data)

    B_Sim2 = []
    for index, val in enumerate(angle_RANGE):
        B_Sim2.append(func(val, reference_data[index]))
    B_Sim2 = np.asarray(B_Sim2)

    delta = np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]) 
    print("Magnituden Unterschied von: ", delta, ", maxBresDelta: ",
          maxBresDelta, " bigger than magnitude? ", delta > maxBresDelta)

    pool = None
    if delta > maxBresDelta:
        # if delta is too big start iterating
        Sim_result = iteration(B_Sim, B_Sim2, B_RES, fixed_params,reference_data, F,maxBresDelta, model, params_Fit, angle_RANGE, fit_params, pool, is_OOP)
        # Sim_result is now a list 5D List: 
        #    [:, 0] = B_Sim
        #    [:, 1] = Eq_Angle1
        #    [:, 2] = Eq_Angle2
        #    [:, 3] = Eq_Angle3
        #    [1] = best fit params dict
        B_Sim2 = Sim_result[0]

        # make output file
        now = datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        end_pfad_now = end_pfad + str(now)

        with open(end_pfad_now + '_Ani_Constants.dat','w') as f:
            for key, val in Sim_result[1].items():
                print([key, val])
                f.write('{0}: {1}'.format(key,val))
                f.write('\n')
        np.savetxt(end_pfad_now + '_Function.dat', B_Sim2[:,0])
        Eq_angles = [B_Sim2[:,1],B_Sim2[:,2],B_Sim2[:,3]]
        np.savetxt(end_pfad_now + '_Eq_Angles.dat', Eq_angles)
    end = time.time() - start
    print('It took',end)

    if plot:
        #phi_deg = pool.map(rad_to_deg, angle_RANGE)
        phi_deg = []
        for val in angle_RANGE:
            phi_deg.append(rad_to_deg(angle_RANGE))
        phi_deg = np.asarray(phi_deg[0])

        plt.cla()
        plt.clf()
        plt.plot(phi_deg, B_Sim2[:, 0], '-g', label='Simulation')
        # plt.plot(angle_RANGE_deg, ani_result.best_fit, '-r', label='best fit')
        #plt.plot(phi_deg, B_Sim[:, 0], '-g', label='Simulation1')
        plt.plot(phi_deg, reference_data, '-b', label='Interpolation')
        plt.scatter(phi_deg, reference_data, label='Experiment')
        plt.xlabel('Angle [Degrees]')
        plt.ylabel('B-Field [Arb. Units]')
        plt.legend()
        plt.show()

def rad_to_deg(val):
    val = val*180/m.pi
    return val

def deg_to_rad(val):
    val = val*m.pi/180
    return val

def make_phirange(shift_value: float, phi_array: list, circ: bool, minWinkel: float, maxWinkel: float):
    #if deg:
    delta_winkel = maxWinkel - minWinkel
    if delta_winkel >= 350: # Check if angluar dependence is from 0-360 deg -> Circle
        for l, i in enumerate(phi_array):   # shift array by value: shift_value
            phi_array[l] += shift_value

             # Check if phi value is over or under the limit of 0.0 - 360.0 degrees
            if phi_array[l] > maxWinkel:
                phi_array[l] -= delta_winkel
            elif phi_array[l] < minWinkel:
                phi_array[l] += delta_winkel

    else:   # Angular dependences without 360 degrees need to be handled differently
        for l, i in enumerate(phi_array):   # shift array by value: shift_value
            phi_array[l] += shift_value
            if phi_array[l] > 359.99:
                phi_array[l] -= 360.0
            elif phi_array[l] < 0.0:
                phi_array[l] += 360.0
    return phi_array

def make_phi_range():
    if add == '-':
        angle_RANGE = np.arange(angle_min + shift, angle_max + shift, m.pi / (1 / (angle_step / m.pi) - lauf_var),
                                dtype='float')
    else:
        angle_RANGE = np.arange(angle_min + shift, angle_max + shift, m.pi / (1 / (angle_step / m.pi) + lauf_var),
                                dtype='float')
    if len(angle_RANGE) != len_ref:
        print(len(angle_RANGE), len_ref, count)
        if count < 10:
            correct_angle_RANGE(angle_RANGE, len_ref, count + 1, '-', lauf_var + 1)
        elif count < 20:
            if count == 10:
                lauf_var = 0
            correct_angle_RANGE(angle_RANGE, len_ref, count + 1, '+', lauf_var + 1)
        else:
            print('angle_RANGE and reference_data have different length, change Anglestep (~ +-1 ) and try again')
    else:
        return angle_RANGE




if __name__ == '__main__':

    #For debug only
    
    #text = "C:/Users/Jonas/Desktop/test.dat B*M*(sin(theta)*sin(thetaB)*cos(phi - phiB) + cos(theta)*cos(thetaB)) - K2p*sin(theta)**2*cos(phi - phiU)**2 - K4p*(cos(4*phi) + 3)*sin(theta)**4/8 - K4s*cos(theta)**4/2 - (-K2s + M**2*mu0/2)*sin(theta)**2 {'K2p': 863.25, 'K2s': 261345.0, 'K4p': 13720.6, 'phiu': 5.0756} {'omega': 62066561101.381386, 'g': 2.05, 'M': 1530000.0, 'K4s': 0} 0.0 0.03490658503988659 False False"
    init_load('C:/Users/Jonas/Desktop/test.dat','-B*M*(sin(theta)*sin(thetaB)*cos(phi - phiB) + cos(theta)*cos(thetaB)) - K2p*sin(theta)**2*cos(phi - phiU)**2 - K4p*(cos(4*phi) + 3)*sin(theta)**4/8 - K4s*cos(theta)**4/2 - (-K2s + M**2*mu0/2)*sin(theta)**2', {'K2p': 863.25, 'K2s': 100000.0, 'K4p': 13720.6, 'phiU': 5.0756}, {'omega': 62066561101.381386, 'g': 2.05, 'M': 1530000.0, 'K4s': 0}, 37.5, 0.03490658503988659,True, True,[332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356,   0,
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
