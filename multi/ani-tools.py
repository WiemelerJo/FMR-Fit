from lmfit import Model
from lmfit import Parameters
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


class Ani_Fit():
    def __init__(self,fname,Ani_model):
        self.fname = fname  #Parameter Filename
        self.ani_model = Ani_model  #Index of a dropdown menu from GUI for selecting different Ani models
        self.setup()

    def setup(self):
        self.K2p = Symbol('K2p')
        self.K2s = Symbol('K2s')
        self.K4p = Symbol('K4p')
        self.K4s = Symbol('K4s')
        self.M = Symbol('M')
        self.mu0 = Symbol('mu0')
        self.mub = Symbol('mub')
        self.B = Symbol('B')
        self.theta = Symbol('theta')
        self.thetaB = Symbol('thetaB')
        self.phi = Symbol('phi')
        self.phiU = Symbol('phiU')
        self.phiB = Symbol('phiB')
        self.gamma = Symbol('gamma')
        self.omega = Symbol('omega')
        self.w = Symbol('w')
        self.g = Symbol('g')
        self.hbar = Symbol('hbar')
        self.f = Symbol('f')
        self. F = Function('F')

        K2s_WERT = 150000  # 120000#
        K2p_WERT = -1287  # 1200#
        K4s_WERT = 1  # 0
        K4p_WERT = 2000  # 2000#
        phiU_WERT = 2.15  # 2.12#2.13#
        # phiB_WERT = 2*m.pi

        self.maxBresDelta = 0.01
        shift = 0 # for first evaluation set to 0

        # consts Einträge
        self.consts = np.zeros(8)
        self.consts[0] = 9.4567 * 10 ** 9  # Frequenz
        self.consts[1] = 2.05  # g
        self.consts[2] = 1.1890 * 10 ** 6  # M
        self.consts[3] = 4 * m.pi * 10 ** (-7)  # mu0
        self.consts[4] = 9.27 * 10 ** (-24)  # mub
        self.consts[5] = (6.62606957 * 10 ** (-34)) / (2 * m.pi)  # hbar
        self.consts[6] = self.consts[1] * self.consts[4] / self.consts[5]
        self.consts[7] = 2 * m.pi * self.consts[0]  # omega bzw. w in der Formel

        # params Einträge
        self.params = np.zeros(5)
        self.params[0] = K2s_WERT
        self.params[1] = K2p_WERT
        self.params[2] = K4s_WERT
        self.params[3] = K4p_WERT
        self.params[4] = phiU_WERT

        # Laden der gefitteten Resonanzpositionen
        D = np.loadtxt("parameter-IP-19-09-2020.dat", dtype='float', skiprows=0)
        self.R_raw = D[:, 3]
        self.Winkel = D[:, 5].flatten()
        self.Winkel2 = np.copy(self.Winkel)
        self.Winkel2 = self.shifting(self.Winkel2, shift, 'deg')
        self.B_inter = interp1d(self.Winkel, self.R_raw)  # Interpoliertes B_res, array mit Länge len(Winkel)

        # Freie Energy Formel + Baselgia
        self.F = self.M * self.B * (sin(self.theta) * sin(self.thetaB) * cos(self.phi - self.phiB) + cos(self.theta) * cos(self.thetaB)) - (
                1 / 2 * self.mu0 * self.M ** 2 - self.K2s) * sin(self.theta) ** 2 - self.K2p * sin(self.theta) ** 2 * cos(
            self.phi - self.phiU) ** 2 - 1 / 2 * self.K4s * cos(self.theta) ** 4 - 1 / 8 * self.K4p * (3 + cos(4 * self.phi)) * sin(self.theta) ** 4
        # F = M*B*(sin(theta)*sin(thetaB)*cos(phi-phiB)+cos(theta)*cos(thetaB))
        # -(1/2*mu0*M**2-K2s)*sin(theta)**2-K2p*sin(theta)**2*cos(phi-phiU)**2
        # -1/2*K4s*cos(theta)**4-1/8*K4p*(3+cos(4*phi))*sin(theta)**4

        halb_res = self.gamma ** 2 / self.M ** 2 * (
                self.F.diff(self.theta, 2) * (self.F.diff(self.phi, 2) / (sin(self.theta)) ** 2 + cos(self.theta) / sin(self.theta) * self.F.diff(self.theta)) - (
                self.F.diff(self.theta, self.phi) / sin(self.theta) - cos(self.theta) / sin(self.theta) * self.F.diff(self.phi) / sin(self.theta)) ** 2)
        eq_halb_res = Eq(self.omega ** 2, halb_res)
        # Formel der Resonanzfrequenz B_RES, Lösung [0] ist positive Lösung der Wurzelfunktion
        B_RES = solve(eq_halb_res, self.B)  # Diese Funktion ist ResField(...) aus Mathematica!
        self.B_RES = B_RES.args[1]

        self.Bounds = [(0.001, 2 * m.pi), (0.001, 2 * m.pi), (0.001, 2 * m.pi)]  # Zum minimieren

        params_name = [self.K2s, self.K2p, self.K4s, self.K4p, self.phiU]
        const_name = [self.f, self.g, self.M, self.mu0, self.mub, self.hbar, self.gamma, self.omega]
        names = [params_name, const_name]
        self.names = [val for sublist in names for val in sublist]
        # names is flattened[params,consts] : ----->  [K2s, K2p, K4s, K4p, phiU, f, g, M, mu0, muB, hbar, gamma, omega]

        rule = [self.params, self.consts]
        rule = [val for sublist in rule for val in sublist]
        # rule is flattened values of the names

        # rules is array names and values combined
        self.rules = []
        for i in range(len(rule)):
            self.rules.append([self.names[i], rule[i]])

        self.rule_consts = []
        for i in range(len(const_name)):
            self.rule_consts.append([const_name[i], self.consts[i]])

        self.params_Fit = Parameters()
        self.params_Fit.add_many(
            ('K2S', self.params[0], True, -400000, 400000, None, None),
            ('K2P', self.params[1], True, -50000, 50000, None, None),
            ('K4S', self.params[2], False, -400000, 400000, None, None),
            ('K4P', self.params[3], True, -50000, 50000, None, None),
            ('PHI_U', self.params[4], True, 0, 2 * m.pi, None, None)
        )
        self.model = Model(self.Model_Fit_Fkt)


        self.gui_shift(-25) #remove when implemented in GUI!!!----------------------------------------------------------------------------



    def gui_shift(self,d):
        #called from DoubleSpinBox while selecting the shift
        # d is the value of the spinbox returned by the connected signal

        #shift = -25
        shift = d
        phi_min = min(self.Winkel) * m.pi / 180
        phi_max = max(self.Winkel) * m.pi / 180
        phi_step = m.pi / 91
        self.phi_RANGE = np.arange(phi_min, phi_max, phi_step, dtype='float')  # array of Rad
        self.phi_RANGE_deg = np.multiply(self.phi_RANGE, 180 / m.pi)

        phi_shifted = np.copy(self.phi_RANGE)
        self.phi_shifted = self.shifting(phi_shifted, shift, 'rad')
        self.phi_shifted_deg = np.multiply(phi_shifted, 180 / m.pi)

        maxWinkel = max(self.Winkel)
        minWinkel = min(self.Winkel)
        for i in range(len(self.phi_shifted_deg)):
            if self.phi_shifted_deg[i] > maxWinkel:
                self.phi_shifted_deg[i] -= maxWinkel
            elif self.phi_shifted_deg[i] < minWinkel:
                self.phi_shifted_deg[i] += maxWinkel


    def shifting(self, array, shift, deg):
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

    def solveAngles(self,B, phiB, fkt):
        # Minimize Angles Theta, Phi, ThetaB
        # start_Paras for the moment these are constant: [m.pi/2, m.pi, m.pi/2]
        # start_Paras = [m.pi/2, m.pi, m.pi/2]
        # --------------------------------------Finding Global Minimum (Mathematica approach)---------------------------------------------
        #Extremly slow, slower than using a global solver in minimize function
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

    def ResFieldNumInp(self, rules, phi_val):
        # takes EqAngles and solves B_RES using interpolated data
        B = self.B_inter(phi_val * 180 / m.pi)
        phiB_val = phi_val

        B_FKT = self.B_RES.subs({i: l for i, l in rules})  # rules beeing inserted, missing : B, theta, thetaB, phi, phiB
        FKT = self.F.subs({i: l for i, l in rules})  # rules beeing inserted, for minimizing
        self.Eq_angles = self.solveAngles(B, phiB_val, FKT)
        result = B_FKT.subs({self.theta: self.Eq_angles[0], self.phi: self.Eq_angles[1], self.thetaB: self.Eq_angles[2], self.phiB: phiB_val}) #calculate B_RES
        return [float(result), float(self.Eq_angles[0]), float(self.Eq_angles[1]), float(self.Eq_angles[2])]

    def F_fkt(self,x, *args):
        # Function of Free Energy Density, that should be minimized
        # *args is array:[B,phiB]
        # x are the fitted params
        Fit_fkt = args[2].subs({self.B: args[0], self.phiB: args[1], self.theta: x[0], self.phi: x[1], self.thetaB: x[2]})
        return float(Fit_fkt)

    def Model_Fit_Fkt(self,phi_real, K2S, K2P, K4S, K4P, PHI_U):
        # is the function used for least square fitting
        fkt = []
        a = self.B_RES.subs({i: l for i, l in self.rule_consts})  # Parameter, B and EqAngle missing
        for i, l in enumerate(phi_real):
            b = a.subs({self.phiB: l, self.K2s: K2S, self.K2p: K2P, self.K4s: K4S, self.K4p: K4P, self.phiU: PHI_U})
            c = b.subs({self.theta: self.Eq_angles[0][i], self.phi: self.Eq_angles[1][i], self.thetaB: self.Eq_angles[2][i]})
            fkt.append(float(c))
        return fkt

    def update_rules(self,params_new):
        # update rules array according to new values inside params_new
        rule_new = [params_new, self.consts]
        rule_new = [val for sublist in rule_new for val in sublist]
        rules_new = []
        for i in range(len(rule_new)):
            rules_new.append([self.names[i], rule_new[i]])
        print('Updated rules:', rules_new)
        return rules_new

    def iterate(self):
        pool = Pool()  # Create pool of threads for multiprocessing
        func = partial(self.ResFieldNumInp, self.rules)  # indirectly add more than one functional arguments without *args or **kwargs
        B_Sim = pool.map(func, self.phi_RANGE)  # Now this function can be mapped to the pool using an array phi_RANGE
        B_Exp = pool.map(self.B_inter, self.phi_shifted_deg)  # Same as above
        B_Sim = np.asarray(B_Sim)  # convert List to numpy array
        self.Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

        Fit = True  # Boolean to determine whether you want to fit or just display data

        if Fit:
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

            ani_result = self.model.fit(B_Exp, self.params_Fit, phi_real=self.phi_RANGE)
            print(ani_result.fit_report())

            temp_paras = []
            for name, param in ani_result.params.items():
                temp_paras.append(float(param.value))
                self.params_Fit[name].set(value=param.value)
            # temp_paras ist nun array aus Endwerten des Fits: [K2s,K2p,K4s,K4p,PhiU]
            # Aktualisiere nun rules
            for i, l in enumerate(temp_paras):
                self.params[i] = l

            rules = self.update_rules(self.params)  # update rules
            func = partial(self.ResFieldNumInp, rules)
            B_Sim2 = pool.map(func, self.phi_RANGE)
            B_Sim2 = np.asarray(B_Sim2)
            print("Magnituden Unterschied von: ", np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]), ", maxBresDelta: ",
                  self.maxBresDelta, " bigger than magnitude? ", np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]) > self.maxBresDelta)
            it_while = 0
            while np.linalg.norm(B_Sim - B_Sim2) > self.maxBresDelta:
                print('Error too big! Start iterating')
                if np.array_equal(B_Sim,B_Sim2):
                    # check if B_Sim changed
                    it_while += 1
                    print('Error: Difference between simulated Bres did not change! it_while = ', it_while)
                    if it_while == 5:
                        break
                B_Sim = np.copy(B_Sim2)
                self.Eq_angles = [B_Sim[:, 1], B_Sim[:, 2], B_Sim[:, 3]]

                ani_result = self.model.fit(B_Exp, self.params_Fit, phi_real=self.phi_RANGE)
                print(ani_result.fit_report())

                temp_paras = []
                for name, param in ani_result.params.items():
                    temp_paras.append(float(param.value))
                    self.params_Fit[name].set(value=param.value)
                # temp_paras ist nun array aus Endwerten des Fits: [K2s,K2p,K4s,K4p,PhiU]
                # Aktualisiere nun rules
                for i, l in enumerate(temp_paras):
                    self.params[i] = l

                rules = self.update_rules(params)
                func = partial(self.ResFieldNumInp, rules)
                B_Sim2 = pool.map(func, self.phi_RANGE)
                B_Sim2 = np.asarray(B_Sim2)
                print('Error in this iteration = ', np.linalg.norm(B_Sim[:, 0] - B_Sim2[:, 0]))

            plt.plot(self.phi_RANGE, B_Sim2[:, 0], '-g', label='Simulation2')
            plt.plot(self.phi_RANGE, ani_result.best_fit, '-r', label='best fit')
        plt.plot(self.phi_RANGE, B_Sim[:, 0], '-g', label='Simulation1')
        plt.plot(self.phi_RANGE, B_Exp, '-b', label='Interpolation')
        plt.scatter(self.phi_RANGE, B_Exp, label='Experiment')
        plt.xlabel('Angle [Degrees]')
        plt.ylabel('B-Field [Arb. Units]')
        plt.legend()
        plt.show()



if __name__ == '__main__':
    Ani_Fit(123, 2)
    Ani_Fit(123,2).iterate()


# Ohne Global Minimum search 8,6 sek.
# Shgo: 19.2 sec.