#Skript for generating the Functions as a String
head = str()
function = str()

'''for i in range(1,11):
    head_new = " dB{0}, R{0}, A{0},".format(i)
    function_new = " + lorentz{0}(B, dB{0}, R{0}, A{0})".format(i)
    function = function + function_new
    head = head + head_new
    string = "if fit_num == {0}:\n\tdef combined(B, slope,offset,{1}):\n\t\tcombi = {2}\n\t\treturn combi".format(i,head,function)
    print(string)
    '''

for i in range(1,11):
    string = "def lorentz{0}(B, dB{0}, R{0}, A{0}):\n\treturn(-64*A{0}*(B-R{0})/(9*dB{0}*(1+(B-R{0})**2/(dB{0}*m.sqrt(3)/4)**2)**2))".format(i)
    print(string)

'''def lorentz1(B, dB1, R1, A1):
    return (-64*A1*(B-R1)/(9*dB1*(1+(B-R1)**2/(dB1*m.sqrt(3)/4)**2)**2))'''