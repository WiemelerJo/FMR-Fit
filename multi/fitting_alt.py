from lmfit import Model
from lmfit import Parameters
import math as m


paramsD = Parameters()
paramsL = Parameters()
'''
-----------------------------Bessere version um Funktionen zu erstellen, wird jedoch nicht von lmfit unterstützt!----------------------------------
def dyson(B,*args):
    print(B,args)
    print(args[2])
    return_stmt = ((4*args[3]*args[1]**2*(3*args[0]*args[1]-4*args[0]*(B-args[2])**2-8*m.sqrt(3)*args[1]*(B-args[2])))/(m.sqrt(3)*(4*(B-args[2])**2+3*args[1]**2)**2))
    if fit_num == 1:
        return(return_stmt+args[4]*B+args[5]) #args[4] und [5] sind hier slope und offset
    else:   
        for num in range(2,fit_num+1):
            return_stmt += ((4*args[3+(num-1)*4]*args[1+(num-1)*4]**2*(3*args[0+(num-1)*4]*args[1+(num-1)*4]-4*args[0+(num-1)*4]*(B-args[2+(num-1)*4])**2-8*m.sqrt(3)*args[1+(num-1)*4]*(B-args[2+(num-1)*4])))/(m.sqrt(3)*(4*(B-args[2+(num-1)*4])**2+3*args[1+(num-1)*4]**2)**2))
        return(return_stmt+args[0+(fit_num)*4]*B+args[1+(fit_num)*4])

def lorentz(B,*args):
    return_stmt = -(args[2]*(2*(B-args[0])))/(m.pi*args[1]**3*((B-args[0]**2)/(args[1]**2)+1)**2)
    if fit_num == 1:
        return(return_stmt+args[3]*B+args[4]) #args[3] und [4] sind hier slope und offset
    else:   
        for num in range(2,fit_num+1):
            return_stmt += (-(args[2+(num-1)*3]*(2*(B-args[0+(num-1)*3])))/(m.pi*args[1+(num-1)*3]**3*((B-args[0+(num-1)*3]**2)/(args[1+(num-1)*3]**2)+1)**2))
        return(return_stmt+args[0+(fit_num)*3]*B+args[1+(fit_num)*3])
'''

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






def lorentz1(B, dB1, R1, A1, offset, slope):
    return(-(A1*(2*(B-dB1)))/(m.pi*R1**3*((B-dB1**2)/(R1**2)+1)**2))+slope*B+offset
def lorentz2(B, dB1, R1, A1, dB2, R2, A2, offset, slope):
    return(-(A1*(2*(B-dB1)))/(m.pi*R1**3*((B-dB1**2)/(R1**2)+1)**2))+slope*B+offset+(-(A2*(2*(B-dB2)))/(m.pi*R2**3*((B-dB2**2)/(R2**2)+1)**2))+slope*B+offset
def lorentz3(B, dB1, R1, A1, dB2, R2, A2, dB3, R3, A3, offset, slope):
    return(-(A1*(2*(B-dB1)))/(m.pi*R1**3*((B-dB1**2)/(R1**2)+1)**2))+slope*B+offset+(-(A2*(2*(B-dB2)))/(m.pi*R2**3*((B-dB2**2)/(R2**2)+1)**2))+(-(A3*(2*(B-dB3)))/(m.pi*R3**3*((B-dB3**2)/(R3**2)+1)**2))+slope*B+offset


class Fit(object):
    """docstring for Fit"""
    def __init__(self, index_model, fit_num,Adata2,Bdata2, j_min,j):
        #super(Fit, self).__init__()
        print('Gestartet')

    def __make_params__(self,fit_num):
        global paramsD
        global paramsL
        for numbers in range(1,fit_num+1):
            paramsD.add_many(
            ('alpha'+str(numbers), 0.001, True, 0, 1, None, None),
            ('dB'+str(numbers),0.034, True, 0, 1500, None, None),
            ('R'+str(numbers),0.1, True, 0, 1.2, None, None),
            ('A'+str(numbers),7, True, 0, 150, None, None),
                )
            paramsL.add_many(   
            ('dB'+str(numbers),0.034, True, 0, 1500, None, None),
            ('R'+str(numbers),0.1, True, 0, 1.2, None, None),
            ('A'+str(numbers),7, True, 0, 150, None, None),
                )

        paramsD.add_many(   
            ('offset',2, True, -5, 15, None, None),
            ('slope',0.5, True, -9, 9, None, None)
            )
        paramsL.add_many(   
            ('offset',2, True, -5, 15, None, None),
            ('slope',0.5, True, -9, 9, None, None)
            )


    def fit(self,i, index_model, fit_num,Adata2,Bdata2, j_min,j):
        model = Model(self.combine(index_model, fit_num))
        if index_model == 2:
            result = model.fit(Adata2[j_min:j], paramsL, B=Bdata2[j_min:j])
        elif index_model == 3:
            result = model.fit(Adata2[j_min:j], paramsD, B=Bdata2[j_min:j])
        else:
            result = model.fit(Adata2[j_min:j], paramsL, B=Bdata2[j_min:j])
        return result

    def run(self,i, index_model, fit_num,Adata2,Bdata2, j_min,j):
        self.__make_params__(fit_num)
        #self.fit(i, index_model, fit_num,Adata2,Bdata2, j_min,j)
        print(self.fit(i, index_model, fit_num,Adata2,Bdata2, j_min,j).fit_report())
        print('Wörk') 

    def combine(self, index_model, fit_num):
        if index_model == 1:
            print("Gaussian not supported yet, switching to Lorentz")
            combi = lorentz1
        elif index_model == 2:
            if fit_num == 1:
                combi = lorentz1 
            elif fit_num == 2:
                combi = lorentz2
            else:
                combi = lorentz3
        elif index_model == 3:
            if fit_num == 1:
                combi = dyson1
            elif fit_num == 2:
                combi = dyson2
            else:
                combi = dyson3
        else:
            combi = lorentz1
        return combi








'''

def dyson_func(B, alpha, dB, R, A, slope, offset):
    return ((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))+slope*B+offset

def dyson(B, alpha, dB, R, A):
    return ((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))
def dyson2(B, alpha2, dB2, R2, A2):
    return ((4*A2*dB2**2*(3*alpha2*dB2-4*alpha2*(B-R2)**2-8*m.sqrt(3)*dB2*(B-R2)))/(m.sqrt(3)*(4*(B-R2)**2+3*dB2**2)**2))
def dyson3(B, alpha3, dB3, R3, A3):
    return ((4*A3*dB3**2*(3*alpha3*dB3-4*alpha3*(B-R3)**2-8*m.sqrt(3)*dB3*(B-R3)))/(m.sqrt(3)*(4*(B-R3)**2+3*dB3**2)**2))


def linear(B, slope, offset):
    return slope*B+offset


def lorentz_func(B,dB, R, A, slope, offset):
    return(-(A*(2*(B-dB)))/(m.pi*R**3*((B-dB**2)/(R**2)+1)**2))+slope*B+offset

def lorentz(B,dB, R, A):
    return(-(A*(2*(B-dB)))/(m.pi*R**3*((B-dB**2)/(R**2)+1)**2))
def lorentz2(B,dB2, R2, A2):
    return(-(A2*(2*(B-dB2)))/(m.pi*R2**3*((B-dB2**2)/(R2**2)+1)**2))
def lorentz3(B,dB3, R3, A3):
    return(-(A3*(2*(B-dB3)))/(m.pi*R3**3*((B-dB3**2)/(R**2)+1)**2))


    def combine(self,index_model,fit_num):
        if index_model == 1:
            print("Gaussian not supported yet, switching to Lorentz")
            combi = lorentz + linear
        elif index_model == 2:
            if fit_num == 1:
                combi = lorentz + linear
            elif fit_num == 2:
                combi = lorentz2 + linear
            else:
                combi = lorentz3 + linear
        elif index_model == 3:
            if fit_num == 1:
                combi = dyson + linear
            elif fit_num == 2:
                combi = dyson2 + linear
            else:
                combi = dyson3 + linear
        else:
            combi = Lorentz + linear
        return combi



if fit_num ==1:
	paramsD.add_many(
	('alpha', 0.001, True, 0, 1, None, None),
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None)
	)
elif fit_num ==2:
	paramsD.add_many(
	('alpha', 0.001, True, 0, 1, None, None),
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None),
	('alpha2', 0.001, True, 0, 1, None, None),
	('dB2',0.04, True, 0, 1500, None, None),
	('R2',0.08, True, 0, 1.2, None, None),
	('A2',5, True, 0, 150, None, None)
	)
elif fit_num ==3:
	paramsD.add_many(
	('alpha', 0.001, True, 0, 1, None, None),
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None),
	('alpha2', 0.001, True, 0, 1, None, None),
	('dB2',0.05, True, 0, 1500, None, None),
	('R2',0.06, True, 0, 1.2, None, None),
	('A2',3, True, 0, 150, None, None),
	('alpha3', 0.001, True, 0, 1, None, None),
	('dB3',0.025, True, 0, 1500, None, None),
	('R3',0.06, True, 0, 1.2, None, None),
	('A3',5, True, 0, 150, None, None)
	)

elif fit_num ==1:
	paramsL.add_many(
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None)
	)
elif fit_num ==2:
	paramsL.add_many(
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None),
	('dB2',0.04, True, 0, 1500, None, None),
	('R2',0.08, True, 0, 1.2, None, None),
	('A2',5, True, 0, 150, None, None)
	)
elif fit_num ==3:
	paramsL.add_many(
	('dB',0.034, True, 0, 1500, None, None),
	('R',0.1, True, 0, 1.2, None, None),
	('A',7, True, 0, 150, None, None),
	('offset',2, True, -5, 15, None, None),
	('slope',0.5, True, -9, 9, None, None),
	('dB2',0.05, True, 0, 1500, None, None),
	('R2',0.06, True, 0, 1.2, None, None),
	('A2',3, True, 0, 150, None, None),
	('dB3',0.025, True, 0, 1500, None, None),
	('R3',0.06, True, 0, 1.2, None, None),
	('A3',5, True, 0, 150, None, None)
	)
'''