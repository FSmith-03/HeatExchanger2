import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute


def brute_force(params, bounds):
     # Extract parameters
    Nb, Y, mdot1, mdot2 = params
    a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
    cp = 4179
    F = 1
    N = 4
    while N < 11:
        A = N*np.pi*0.006*0.35
    
        v_tube = v_tube_finder(mdot2, N)
        reynolds_tube = reynolds_tube_finder(v_tube)
        v_nozzle_2 = v_nozzle_finder(mdot2)
        pressure_loss_tube = presure_loss_tube_finder(v_tube, reynolds_tube)
        sigma = sigma_finder(N)
        kc, ke = kc_ke_finder(sigma)
        pressure_loss_ends = pressure_loss_ends_finder(v_tube, kc, ke)
        pressure_loss_nozzle_2 = pressure_loss_nozzle_finder(mdot2)

        pressure_loss_2 = (pressure_loss_ends + pressure_loss_nozzle_2 + pressure_loss_tube) / 1e5
        pressure_rise_2, state2 = pressure_checker(pressure_loss_2, mdot2, 2)

        # Cold stream
        A_sh = shell_area_finder(Y, Nb)
        v_shell = v_shell_finder(mdot1, A_sh)
        v_nozzle_1 = v_nozzle_finder(mdot1)
        d_sh = shell_chic_length_finder(A_sh)
        reynolds_shell = reynolds_shell_finder(v_shell, d_sh)
        pressure_loss_shell = pressure_loss_shell_finder(reynolds_shell, N, v_shell, 0.2)
        pressure_loss_nozzle_1 = pressure_loss_nozzle_finder(v_nozzle_1)

        pressure_loss_1 = (pressure_loss_shell + pressure_loss_nozzle_1) / 1e5
        pressure_rise_1, state1 = pressure_checker(pressure_loss_1, mdot1, 1)
    
        # Thermal Analysis
        H = H_finder(reynolds_tube, reynolds_shell)
        T1out, T2out, LMTD = temperature_solver2(mdot1, mdot2, H, A, F)
