import math as m

def Lorentz(B, dB, R, A):
    return(-64*A*(B-R)/(9*dB*(1+(B-R)**2/(dB*m.sqrt(3)/4)**2)**2))

def Dyson(B, alpha, dB, R, A):
    return((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))

def functions_value(B, slope, offset, index_model, params,fit_num):
    #params should now be the missing repeating parameters
    #then differentiate between Lorentz and Dyson
    final = 0

    if index_model == 2:
        #Lorentz-Model
        #generate the variables
        for i in range(1,fit_num+1):
            final += Lorentz(B,params[3*(i-1)],params[1+3*(i-1)],params[2+3*(i-1)])
        return final + B * slope + offset 
    else:
        #Dyson-Model
        #generate the variables

        for i in range(1,fit_num+1):
            final += Dyson(B,params[4*(i-1)],params[1+4*(i-1)],params[2+4*(i-1)],params[3+4*(i-1)])
        return final + B * slope + offset

def single_func(B, slope, offset, index_model, params, fit_num):
    try:
        if index_model == 2:
            #Lorentz-Model
            final = Lorentz(B,params[3*(fit_num-1)],params[1+3*(fit_num-1)],params[2+3*(fit_num-1)])
            return final + B * slope + offset 
        else:
            #Dyson-Model
            final = Dyson(B,params[4*(fit_num-1)],params[1+4*(fit_num-1)],params[2+4*(fit_num-1)],params[3+4*(fit_num-1)])
            return final + B * slope + offset
    except Exception as e:
        print("Error in single_func:",e)
        print("Maybe wrong number of fcuntions")