import symengine
from sympy import sin,cos,sqrt,Eq,Symbol,Function,Matrix,re

A1 = Symbol('A1')
dB1 = Symbol('dB1')
alpha1 = Symbol('alpha1')
R1 = Symbol('R1')
B = Symbol('B')

Dyson = (4*A1*dB1**2*(3*alpha1*dB1-4*alpha1*(B-R1)**2-8*sqrt(3)*dB1*(B-R1)))/(sqrt(3)*(4*(B-R1)**2+3*dB1**2)**2)

print(
    Dyson.diff(R1)

    )