from Hydraulic_Functions import *
import numpy as np
from scipy.optimize import fsolve

# Define fixed values for mdot1, mdot2, H, A, and F
mdot1 = 0.5
mdot2 = 0.6
H = 1000
A = 10
F = 5
cp = 4179
def equation(T1o):
    T1in = 20
    T2in = 60
    cp = 4179
    T2o = T2in - (mdot1/mdot2)*(T1o - T1in)
    y = H*A*(((T2in - T1o) - (T2o - T1in))/np.log((T2in - T1o)/(T2o - T1in))) - mdot1*cp*(T1o - T1in)
    return y 

# Initial guess for T1o
initial_guess_T1o = 25

result = fsolve(equation, initial_guess_T1o)
print("Root for T1o:", result[0])
print(result - 20)
T2o = 60 - (mdot1/mdot2)*(result - 20)
print(60-T2o)
print((equation(result)- mdot1*cp*(result - 20))/(mdot1*cp))
