import math as m
import matplotlib.pyplot as plt 
import numpy as np
from lmfit import Model
from lmfit import Parameters
from func_gen import Gen_Lorentz, Gen_Dyson
#from PyQt5.QtCore import *

class Fit(object):
    """docstring for Fit"""
    def __init__(self,  index_model, fit_num,Adata2,Bdata2, j_min,j,init_values,bound_min,bound_max):
        super(Fit, self).__init__()
        self.index_model = index_model
        self.__make_params__(fit_num,init_values,bound_min,bound_max)
        self.set_model(fit_num,index_model)

    def set_model(self,fit_num, index_model):
        if index_model == 2:
            fit_func = Gen_Lorentz(fit_num)
        elif index_model == 3:
            fit_func = Gen_Dyson(fit_num)
        else:
            print('Not supported!')
        exec(fit_func.get_str(), globals())
        self.model = Model(model_fit_func)

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
            if self.index_model == 3:
                paramsD.add_many(
                ('alpha'+str(numbers), init_values[4*(numbers-1)+2], True, bound_min[4*(numbers-1)+2], bound_max[4*(numbers-1)+2], None, None),
                ('dB'+str(numbers),init_values[1+4*(numbers-1)+2], True, bound_min[1+4*(numbers-1)+2], bound_max[1+4*(numbers-1)+2], None, None),
                ('R'+str(numbers),init_values[2+4*(numbers-1)+2], True, bound_min[2+4*(numbers-1)+2], bound_max[2+4*(numbers-1)+2], None, None),
                ('A'+str(numbers),init_values[3+4*(numbers-1)+2], True, bound_min[3+4*(numbers-1)+2], bound_max[3+4*(numbers-1)+2], None, None)
                    )
            elif self.index_model == 2:
                paramsL.add_many(
                ('dB'+str(numbers),init_values[3*(numbers-1)+2], True, bound_min[3*(numbers-1)+2], bound_max[3*(numbers-1)+2], None, None),
                ('R'+str(numbers),init_values[1+3*(numbers-1)+2], True, bound_min[1+3*(numbers-1)+2], bound_max[1+3*(numbers-1)+2], None, None),
                ('A'+str(numbers),init_values[2+3*(numbers-1)+2], True, bound_min[2+3*(numbers-1)+2], bound_max[2+3*(numbers-1)+2], None, None),
                    )
            else:
                print('Not supported')

        #for numbers in range(2,len()):

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

    def give_model(self):
        return self.model

    def fit(self,index_model,Adata2,Bdata2, j_min,j):
        global result
        if index_model == 2:
            result = self.model.fit(Adata2[j_min:j], paramsL, B=Bdata2[j_min:j])
        else:
            result = self.model.fit(Adata2[j_min:j], paramsD, B=Bdata2[j_min:j])
        return result

