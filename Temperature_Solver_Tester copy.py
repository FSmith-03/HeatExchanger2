import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution
from temperature_solvers_old import *
import matplotlib.pyplot as plt
def temperature_solverplotter(mdot1, mdot2, H, A, F, N_i):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 57
    T1out = 21
    best_delta_q12 = 1e6
    best_delta_q23 = 1e6
    step = 40/N_i
    qdot3_list = []
    T_List = np.arange(T1in, T2in, step)
    for T2out in T_List:
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        qdot1 = cp*mdot1*(T1out-T1in)
        qdot2 = cp*mdot2*(T2in-T2out)
        qdot3 = H*A*F*LMTD
        delta_q12 = abs(qdot1 - qdot2)
        delta_q23 = abs(qdot2 - qdot3)
        qdot3_list.append(qdot3)
        if delta_q23 < best_delta_q23:
            best_delta_q23 = delta_q23
        else:
            return qdot3_list, qdot3

def temperature_differenceplotter(mdot1, mdot2, H, A, F, N_i):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 57
    T1out = 21
    best_delta_q12 = 1e6
    best_delta_q23 = 1e6
    step = 40/N_i
    T2out_list = []
    Tdifference_list = []
    T1difference_list = []
    T_List = np.arange(T1in, T2in, step)
    T1_list = []
    for T2out in T_List:
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        T1_list.append(T1out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        qdot1 = cp*mdot1*(T1out-T1in)
        qdot2 = cp*mdot2*(T2in-T2out)
        qdot3 = H*A*F*LMTD
        delta_q12 = abs(qdot1 - qdot2)
        delta_q23 = abs(qdot2 - qdot3)
        delta = abs(T2out - T2in)
        delta1 = abs(T1in - T1out)
        Tdifference_list.append(delta)
        T1difference_list.append(delta1)
        if delta_q23 < best_delta_q23:
            best_delta_q23 = delta_q23
            T2out_best = T2out
            delta_best = delta

    return Tdifference_list, T_List, T1difference_list, T1_list, T2out_best, delta_best
            








L, N, Y, N_b, passes = (0.24, 14, 0.012, 6, 2)
A = N*np.pi*0.006*L
F = 1
mdot2, error2 = pressure_intersection_2(N, L, mdot_upper=0.8)
v_tube = v_tube_finder(mdot2, N)
reynolds_tube_value = reynolds_tube_finder(v_tube)

# Cold stream
#Check for shell error from invalid geometries
A_sh = shell_area_finder(Y, N_b, L)
mdot1, error1 = pressure_intersection_1(N, N_b, Y, L, mdot_upper=0.8)

v_shell = v_shell_finder(mdot1, A_sh)
v_nozzle_1 = v_nozzle_finder(mdot1)
d_sh = shell_chic_length_finder(A_sh)
reynolds_shell_value = reynolds_shell_finder(v_shell, d_sh)

# Thermal Analysis
H = H_finder(reynolds_tube_value, reynolds_shell_value, L, N_b, Y, mdot1, passes, N)
#print(mdot1, mdot2, H, A, F)
#T1out, T2out, LMTD = temperature_solvernew2(mdot1, mdot2, H, A, F, 25000, 100)
"""
qdot3_list, qdot3_best = temperature_solverplotter(mdot1, mdot2, H, A, F, 1000)
qdot3_best_list = []
for n in qdot3_list:
    qdot3_best_list.append(qdot3_best)
iteration_list = np.arange(1, len(qdot3_list)+1, 1)
plt.scatter(iteration_list, qdot3_list, s=3, color = 'black', label = 'Qdot Values')
plt.plot(iteration_list, qdot3_best_list, color = "red", label = "Best Qdot value")
plt.xlabel('Iteration Number')
plt.ylabel('Qdot')
plt.title('Qdot solver response')
plt.legend(fontsize = 'large')
plt.show()
"""
Tdifference_List, T_List, T1difference_list, T1_list, T2out, delta = temperature_differenceplotter(mdot1, mdot2, H, A, F, 1000)
plt.scatter(T2out, delta, s=7, color = 'purple', label = 'T2out solution')
plt.plot(T_List, Tdifference_List, color = "red", label = "T2out evolution")
plt.plot(T1_list, T1difference_list, color = "cyan", label = "T1out evolution")
plt.xlabel('Temperature C')
plt.ylabel('Temperarture Difference C')
plt.title('LMTD Solver')
plt.legend(fontsize = 'large')
plt.show()
