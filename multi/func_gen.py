# Skript for generating the Functions as a String
# Idea is, these Classes should create and return a string which in the end via excec(....., globals()) can be turned
# into a python function and in the end into a model for fitting using lmfit

class Gen_Lorentz():
    def __init__(self,fit_num):
        self.num = fit_num

    def get_str(self):
        str_func_args = str() # Set empty str
        str_body = str() # Set empty str
        str_end = "slope*B+offset" # Set linear fkt
        for i in range(1,self.num+1):   # create function/s
            str_func_args += ", dB{0}, R{0}, A{0}".format(i)
            str_body +=  "-64*A{0}*(B-R{0})/(9*dB{0}*(1+(B-R{0})**2/(dB{0}*m.sqrt(3)/4)**2)**2)".format(i)
        func = "def model_fit_func(B,slope,offset{0}):\n\treturn({1} + {2})".format(str_func_args,str_body,str_end)
        return func

class Gen_Dyson():
    def __init__(self,fit_num):
        self.num = fit_num

    def get_str(self):
        str_func_args = str()
        str_body = str()
        str_end = "slope*B+offset"
        for i in range(1,self.num+1):
            str_func_args += ", alpha{0}, dB{0}, R{0}, A{0}".format(i)
            if i == 1:
                str_body += "(4*A{0}*dB{0}**2*(3*alpha{0}*dB{0}-4*alpha{0}*(B-R{0})**2-8*m.sqrt(3)*dB{0}*(B-R{0})))/(m.sqrt(3)*(4*(B-R{0})**2+3*dB{0}**2)**2)".format(i)
            else:
                str_body += "+(4*A{0}*dB{0}**2*(3*alpha{0}*dB{0}-4*alpha{0}*(B-R{0})**2-8*m.sqrt(3)*dB{0}*(B-R{0})))/(m.sqrt(3)*(4*(B-R{0})**2+3*dB{0}**2)**2)".format(
                    i)
        func = "def model_fit_func(B, slope,offset{0}):\n\treturn({1} + {2})".format(str_func_args,str_body,str_end)
        return func

#-------Example on how to get the function string and excecute it to a pyhon function
#obj = Gen_Dyson(10)    # Save Object with 10 Functions to obj
#model = obj.get_str()  # call gen_str() to generate str
#exec(model,globals())   # excecute to python function in global namespace