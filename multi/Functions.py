import math as m

def Lorentz(B, dB, R, A):
    return(-64*A*(B-R)/(9*dB*(1+(B-R)**2)/(dB*m.sqrt(3)/4)**2)**2)

def Dyson(B, alpha, dB, R, A):
    return((4*A*dB**2*(3*alpha*dB-4*alpha*(B-R)**2-8*m.sqrt(3)*dB*(B-R)))/(m.sqrt(3)*(4*(B-R)**2+3*dB**2)**2))

def functions_value(B, slope, offset, index_model, *args):
    #args should now be the missing repeating parameters
    #then differentiate between Lorentz and Dyson
    arg_length = len(args) # count how many arguments where passed
    final = 0

    if index_model == 2:
        #Lorentz-Model
        #generate the variables
        arg_num = arg_length/3 #amount of Functions

        for i in range(1,arg_num+1):
            final += Lorentz(B,args[3*(i_num-1)],args[1+3*(i_num-1)],args[2+3*(i_num-1)])
        return final + B * slope + offset

    else:
        #Dyson-Model
        #generate the variables
        arg_num = int(arg_length/4) #amount of Functions

        for i in range(1,arg_num+1):
            final += Lorentz(B,args[4*(i_num-1)],args[1+4*(i_num-1)],args[2+4*(i_num-1)],args[3+4*(i_num-1)])
        return final + B * slope + offset
