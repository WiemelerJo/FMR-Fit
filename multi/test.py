from lmfit import Model
from lmfit import Parameters

parameter = ['a','b','c','d']
parametercorr = []
fitnum = 3
zähler = 1
paramsD = Parameters()

for index in range(1,fitnum+1):
    for i in parameter:
        parametercorr.append(i+str(index))
        print(i+str(index),zähler)
        zähler += 1
print(parametercorr)
'''
start = [0.03,0.08,7]
start_2 = [0.03,0.08,7]

for l in range(1,10+1):
    temp1 = start[0]
    temp2 = start[1]
    temp3 = start[2]
    start_2.append(temp1+0.01*l)
    start_2.append(temp2+0.01*l)
    start_2.append(temp3+0.5*l)
print(start_2)
'''

default_values_D = [0.0001,0.03,0.08,7, 0.0002,0.03,0.08,7, 0.0003,0.03,0.08,7, 0.0004,0.03,0.08,7, 0.0005,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7]


for numbers in range(1,fitnum+1):
    paramsD.add_many(
    ('alpha'+str(numbers), default_values_D[4*(numbers-1)], True, 0, 1, None, None),
    ('dB'+str(numbers),default_values_D[1+4*(numbers-1)], True, 0, 1500, None, None),
    ('R'+str(numbers),default_values_D[2+4*(numbers-1)], True, 0, 1.2, None, None),
    ('A'+str(numbers),default_values_D[3+4*(numbers-1)], True, 0, 150, None, None))
print(paramsD)
