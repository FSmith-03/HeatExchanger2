import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution
from temperature_solvers_old import *

#N = 13
#Nb = 9
#Y = 0.014
#mdot1 = 0.5
#mdot2 = 0.44
a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
cp = 4179
F = 1
optimalps = [8.04303169, 10.91415054,  0.01623694,  0.55209775,  0.45994879]
N, Nb, Y, mdot1, mdot2 = optimalps
A = N*np.pi*0.006*0.35
L = 0.35
#Hot stream
v_tube = v_tube_finder(mdot2, N)
reynolds_tube = reynolds_tube_finder(v_tube)
v_nozzle_2 = v_nozzle_finder(mdot2)
pressure_loss_tube = presure_loss_tube_finder(v_tube, reynolds_tube, L)
sigma = sigma_finder(N)
#print(sigma)
kc, ke = kc_ke_finder(sigma)
pressure_loss_ends = pressure_loss_ends_finder(v_tube, kc, ke)
pressure_loss_nozzle_2 = pressure_loss_nozzle_finder(mdot2)

pressure_loss_2 = (pressure_loss_ends + pressure_loss_nozzle_2 + pressure_loss_tube)/1e5
pressure_rise_2, state2 = pressure_checker(pressure_loss_2, mdot2, 2)
if state2 == False:
    print("Potential issue, pressure loss in hot stream is more than pressure rise.")

#Cold stream
A_sh = shell_area_finder(Y, Nb, L)
v_shell = v_shell_finder(mdot1, A_sh)
v_nozzle_1 = v_nozzle_finder(mdot1)
d_sh = shell_chic_length_finder(A_sh)
reynolds_shell = reynolds_shell_finder(v_shell, d_sh)
pressure_loss_shell = pressure_loss_shell_finder(reynolds_shell, N, v_shell, 0.2)
pressure_loss_nozzle_1 = pressure_loss_nozzle_finder(v_nozzle_1)
pressure_loss_1 = (pressure_loss_shell + pressure_loss_nozzle_1)/1e5
pressure_rise_1, state1 = pressure_checker(pressure_loss_1, mdot1, 1)
if state1 == False:
    print("Potential issue, pressure loss in cold stream is more than pressure rise.")

#Thermal Analysis
H = H_finder(reynolds_tube, reynolds_shell)
"""
T1out, T2out, LMTD1 = temperature_solver1(mdot1, mdot2, H, A, F)
print("1", LMTD1)
T1out, T2out, LMTD2 = temperature_solver2(mdot1, mdot2, H, A, F)
print("2", LMTD2)
T1out, T2out, LMTD3 = temperature_solver3GPT(mdot1, mdot2, H, A, F)
print("3", LMTD3)
"""
T1out, T2out, LMTD4 = temperature_solvernew2(mdot1, mdot2, H, A, F, 1000)
print("4", T1out, T2out, LMTD4)