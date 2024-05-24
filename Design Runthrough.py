import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution
from temperature_solvers_old import *

def objective_function(L, N, Y, N_b, passes):
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
    H = H_finder(reynolds_tube_value, reynolds_shell_value, L, N_b, Y, mdot1)
    #print(mdot1, mdot2, H, A, F)
    #T1out, T2out, LMTD = temperature_solvernew2(mdot1, mdot2, H, A, F, 25000, 100)
    T1out, T2out, qdot = temperature_solvernew2(mdot1, mdot2, H, A, F, 1000, passes)
    return qdot*passes*2

print(objective_function(0.26, 3, 0.012, 6, 2))

