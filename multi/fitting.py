import math as m
import matplotlib.pyplot as plt 
import numpy as np
from lmfit import Model
from lmfit import Parameters
#from PyQt5.QtCore import *

def func(B,args):
	return ((4*args[3]*args[1]**2*(3*args[0]*args[1]-4*args[0]*(B-args[2])**2-8*m.sqrt(3)*args[1]*(B-args[2])))/(m.sqrt(3)*(4*(B-args[2])**2+3*args[1]**2)**2))

def linear(B,slope,offset):
    return slope*B+offset

def dyson1(B, alpha1, dB1, R1, A1):
    return ((4*A1*dB1**2*(3*alpha1*dB1-4*alpha1*(B-R1)**2-8*m.sqrt(3)*dB1*(B-R1)))/(m.sqrt(3)*(4*(B-R1)**2+3*dB1**2)**2))
def dyson2(B, alpha2, dB2, R2, A2):
    return ((4*A2*dB2**2*(3*alpha2*dB2-4*alpha2*(B-R2)**2-8*m.sqrt(3)*dB2*(B-R2)))/(m.sqrt(3)*(4*(B-R2)**2+3*dB2**2)**2))
def dyson3(B, alpha3, dB3, R3, A3):
    return ((4*A3*dB3**2*(3*alpha3*dB3-4*alpha3*(B-R3)**2-8*m.sqrt(3)*dB3*(B-R3)))/(m.sqrt(3)*(4*(B-R3)**2+3*dB3**2)**2))
def dyson4(B, alpha4, dB4, R4, A4):
    return ((4*A4*dB4**2*(3*alpha4*dB4-4*alpha4*(B-R4)**2-8*m.sqrt(3)*dB4*(B-R4)))/(m.sqrt(3)*(4*(B-R4)**2+3*dB4**2)**2))
def dyson5(B, alpha5, dB5, R5, A5):
    return ((4*A5*dB5**2*(3*alpha5*dB5-4*alpha5*(B-R5)**2-8*m.sqrt(3)*dB5*(B-R5)))/(m.sqrt(3)*(4*(B-R5)**2+3*dB5**2)**2))
def dyson6(B, alpha6, dB6, R6, A6):
    return ((4*A6*dB6**2*(3*alpha6*dB6-4*alpha6*(B-R6)**2-8*m.sqrt(3)*dB6*(B-R6)))/(m.sqrt(3)*(4*(B-R6)**2+3*dB6**2)**2))
def dyson7(B, alpha7, dB7, R7, A7):
    return ((4*A7*dB7**2*(3*alpha7*dB7-4*alpha7*(B-R7)**2-8*m.sqrt(3)*dB7*(B-R7)))/(m.sqrt(3)*(4*(B-R7)**2+3*dB7**2)**2))
def dyson8(B, alpha8, dB8, R8, A8):
    return ((4*A8*dB8**2*(3*alpha8*dB8-4*alpha8*(B-R8)**2-8*m.sqrt(3)*dB8*(B-R8)))/(m.sqrt(3)*(4*(B-R8)**2+3*dB8**2)**2))
def dyson9(B, alpha9, dB9, R9, A9):
    return ((4*A9*dB9**2*(3*alpha9*dB9-4*alpha9*(B-R9)**2-8*m.sqrt(3)*dB9*(B-R9)))/(m.sqrt(3)*(4*(B-R9)**2+3*dB9**2)**2))
def dyson10(B, alpha10, dB10, R10, A10):
    return ((4*A10*dB10**2*(3*alpha10*dB10-4*alpha10*(B-R10)**2-8*m.sqrt(3)*dB10*(B-R10)))/(m.sqrt(3)*(4*(B-R10)**2+3*dB10**2)**2))

def lorentz1(B, dB1, R1, A1):
    return(-64*A1*(B-R1)/(9*dB1*(1+(B-R1)**2/(dB1*m.sqrt(3)/4)**2)**2))
def lorentz2(B, dB2, R2, A2):
    return(-64*A2*(B-R2)/(9*dB2*(1+(B-R2)**2/(dB2*m.sqrt(3)/4)**2)**2))
def lorentz3(B, dB3, R3, A3):
    return(-64*A3*(B-R3)/(9*dB3*(1+(B-R3)**2/(dB3*m.sqrt(3)/4)**2)**2))
def lorentz4(B, dB4, R4, A4):
    return(-64*A4*(B-R4)/(9*dB4*(1+(B-R4)**2/(dB4*m.sqrt(3)/4)**2)**2))
def lorentz5(B, dB5, R5, A5):
    return(-64*A5*(B-R5)/(9*dB5*(1+(B-R5)**2/(dB5*m.sqrt(3)/4)**2)**2))
def lorentz6(B, dB6, R6, A6):
    return(-64*A6*(B-R6)/(9*dB6*(1+(B-R6)**2/(dB6*m.sqrt(3)/4)**2)**2))
def lorentz7(B, dB7, R7, A7):
    return(-64*A7*(B-R7)/(9*dB7*(1+(B-R7)**2/(dB7*m.sqrt(3)/4)**2)**2))
def lorentz8(B, dB8, R8, A8):
    return(-64*A8*(B-R8)/(9*dB8*(1+(B-R8)**2/(dB8*m.sqrt(3)/4)**2)**2))
def lorentz9(B, dB9, R9, A9):
    return(-64*A9*(B-R9)/(9*dB9*(1+(B-R9)**2/(dB9*m.sqrt(3)/4)**2)**2))
def lorentz10(B, dB10, R10, A10):
    return(-64*A10*(B-R10)/(9*dB10*(1+(B-R10)**2/(dB10*m.sqrt(3)/4)**2)**2))



class Fit(object):
    """docstring for Fit"""
    def __init__(self,  index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max):
        super(Fit, self).__init__()
        self.index_model = index_model
        self.__make_params__(fit_num,init_values,bound_min,bound_max)
        self.set_model(fit_num,index_model)

    def combine(self,fit_num,index_model):
        global combined
        if index_model == 3: # Dyson
            if fit_num == 1:
                def combined(B, slope,offset, alpha1, dB1, R1, A1):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + linear(B,slope,offset)
                    return combi
            if fit_num == 2:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + linear(B,slope,offset)
                    return combi
            if fit_num == 3:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + linear(B,slope,offset)
                    return combi
            if fit_num == 4:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + linear(B,slope,offset)
                    return combi
            if fit_num == 5:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + linear(B,slope,offset)
                    return combi
            if fit_num == 6:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5 ,alpha6, dB6, R6, A6):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + dyson6(B, alpha6, dB6, R6, A6) + linear(B,slope,offset)
                    return combi
            if fit_num == 7:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5 ,alpha6, dB6, R6, A6, alpha7, dB7, R7, A7):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + dyson6(B, alpha6, dB6, R6, A6) + dyson7(B, alpha7, dB7, R7, A7) + linear(B,slope,offset)
                    return combi
            if fit_num == 8:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5 ,alpha6, dB6, R6, A6, alpha7, dB7, R7, A7, alpha8, dB8, R8, A8):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + dyson6(B, alpha6, dB6, R6, A6) + dyson7(B, alpha7, dB7, R7, A7) + dyson8(B, alpha8, dB8, R8, A8) + linear(B,slope,offset)
                    return combi
            if fit_num == 9:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5 ,alpha6, dB6, R6, A6, alpha7, dB7, R7, A7, alpha8, dB8, R8, A8, alpha9, dB9, R9, A9):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + dyson6(B, alpha6, dB6, R6, A6) + dyson7(B, alpha7, dB7, R7, A7) + dyson8(B, alpha8, dB8, R8, A8) + dyson9(B, alpha9, dB9, R9, A9) + linear(B,slope,offset)
                    return combi
            if fit_num == 10:
                def combined(B, slope,offset, alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, alpha3, dB3, R3, A3, alpha4, dB4, R4, A4, alpha5, dB5, R5, A5 ,alpha6, dB6, R6, A6, alpha7, dB7, R7, A7, alpha8, dB8, R8, A8, alpha9, dB9, R9, A9, alpha10, dB10, R10, A10):
                    combi = dyson1(B, alpha1, dB1, R1, A1) + dyson2(B, alpha2, dB2, R2, A2) + dyson3(B, alpha3, dB3, R3, A3) + dyson4(B, alpha4, dB4, R4, A4) + dyson5(B, alpha5, dB5, R5, A5) + dyson6(B, alpha6, dB6, R6, A6) + dyson7(B, alpha7, dB7, R7, A7) + dyson8(B, alpha8, dB8, R8, A8) + dyson9(B, alpha9, dB9, R9, A9) + dyson10(B, alpha10, dB10, R10, A10) + linear(B,slope,offset)
                    return combi

        elif index_model == 2: #Lorentz
            if fit_num == 1:
                def combined(B, slope,offset, dB1, R1, A1):
                    combi = lorentz1(B, dB1, R1, A1) + linear(B,slope,offset)
                    return combi
            if fit_num == 2:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + linear(B,slope,offset)
                    return combi
            if fit_num == 3:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + linear(B,slope,offset)
                    return combi
            if fit_num == 4:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + linear(B,slope,offset)
                    return combi
            if fit_num == 5:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + linear(B,slope,offset)
                    return combi
            if fit_num == 6:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5, dB6, R6, A6):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + lorentz6(B, dB6, R6, A6) + linear(B,slope,offset)
                    return combi
            if fit_num == 7:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5, dB6, R6, A6, dB7, R7, A7):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + lorentz6(B, dB6, R6, A6) + lorentz7(B, dB7, R7, A7) + linear(B,slope,offset)
                    return combi
            if fit_num == 8:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5, dB6, R6, A6, dB7, R7, A7, dB8, R8, A8):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + lorentz6(B, dB6, R6, A6) + lorentz7(B, dB7, R7, A7) + lorentz8(B, dB8, R8, A8) + linear(B,slope,offset)
                    return combi
            if fit_num == 9:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5, dB6, R6, A6, dB7, R7, A7, dB8, R8, A8, dB9, R9, A9):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + lorentz6(B, dB6, R6, A6) + lorentz7(B, dB7, R7, A7) + lorentz8(B, dB8, R8, A8) + lorentz9(B, dB9, R9, A9) + linear(B,slope,offset)
                    return combi
            if fit_num == 10:
                def combined(B, slope,offset, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, dB4, R4, A4, dB5, R5, A5, dB6, R6, A6, dB7, R7, A7, dB8, R8, A8, dB9, R9, A9, dB10, R10, A10):
                    combi = lorentz1(B, dB1, R1, A1) + lorentz2(B, dB2, R2, A2) + lorentz3(B, dB3, R3, A3) + lorentz4(B, dB4, R4, A4) + lorentz5(B, dB5, R5, A5) + lorentz6(B, dB6, R6, A6) + lorentz7(B, dB7, R7, A7) + lorentz8(B, dB8, R8, A8) + lorentz9(B, dB9, R9, A9) + lorentz10(B, dB10, R10, A10) + linear(B,slope,offset)
                    return combi         


    def set_model(self,fit_num, index_model):
        global model
        self.combine(fit_num, index_model)
        model = Model(combined)

    def __make_params__(self,fit_num,init_values,bound_min,bound_max):
        global paramsD
        global paramsL
        paramsD = Parameters()
        paramsL = Parameters()

        paramsD.add_many( 
            ('slope',init_values[0], True, bound_min[0], bound_max[0], None, None),             
            ('offset',init_values[1], True, bound_min[1], bound_max[1], None, None)
            )
        paramsL.add_many(  
            ('slope',init_values[0], True, bound_min[0], bound_max[0], None, None),             
            ('offset',init_values[1], True, bound_min[1], bound_max[1], None, None)
            )     

        for numbers in range(1,fit_num+1):
            paramsD.add_many(
            ('alpha'+str(numbers), init_values[4*(numbers-1)+2], True, bound_min[4*(numbers-1)+2], bound_max[4*(numbers-1)+2], None, None),
            ('dB'+str(numbers),init_values[1+4*(numbers-1)+2], True, bound_min[1+4*(numbers-1)+2], bound_max[1+4*(numbers-1)+2], None, None),
            ('R'+str(numbers),init_values[2+4*(numbers-1)+2], True, bound_min[2+4*(numbers-1)+2], bound_max[2+4*(numbers-1)+2], None, None),
            ('A'+str(numbers),init_values[3+4*(numbers-1)+2], True, bound_min[3+4*(numbers-1)+2], bound_max[3+4*(numbers-1)+2], None, None)
                )
            paramsL.add_many(   
            ('dB'+str(numbers),init_values[3*(numbers-1)+2], True, bound_min[3*(numbers-1)+2], bound_max[3*(numbers-1)+2], None, None),
            ('R'+str(numbers),init_values[1+3*(numbers-1)+2], True, bound_min[1+3*(numbers-1)+2], bound_max[1+3*(numbers-1)+2], None, None),
            ('A'+str(numbers),init_values[2+3*(numbers-1)+2], True, bound_min[2+3*(numbers-1)+2], bound_max[2+3*(numbers-1)+2], None, None),
                )

    def give_params(self,fit_num, parameter_table_names_final,index_model,Adata2,Bdata2, j_min,j):
        temp_paras = []
        if index_model == 2:
            #Lorentz
            for name, param in result.params.items():
                temp_paras.append(float(param.value))
        else:
            #Dyson
            for name, param in result.params.items():
                temp_paras.append(float(param.value))
        return temp_paras

    def give_param_table(self,index_model):
        if index_model == 2:
            return paramsL
        else:
            return paramsD

    def give_model(self,index_model):
        if index_model == 2:
            return model
        else:
            return model

    def give_combined(self):
        return combined

    def fit(self,index_model,Adata2,Bdata2, j_min,j):
        global result
        if index_model == 2:
            result = model.fit(Adata2[j_min:j], paramsL, B=Bdata2[j_min:j])
        else:
            result = model.fit(Adata2[j_min:j], paramsD, B=Bdata2[j_min:j])
        return result

