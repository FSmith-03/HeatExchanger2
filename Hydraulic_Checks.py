from Hydraulic_Functions import *


#Check functions against calculated values for different parameters
N, Nb, Y, mdot1, mdot2 = [10, 9, 0.014, 0.5, 0.5]
a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
cp = 4179
F = 1
A = N*np.pi*0.006*0.35


v_t = v_tube_finder(mdot2, N)
if abs(v_t - 1.786) > 0.1:
    print("Tube velocity error")
else:
    pass

Re_t = reynolds_tube_finder(v_t)
if abs(Re_t - 16298) > 10:
    print("Tube reynolds error")
else:
    pass

v_n2 = v_nozzle_finder(mdot2)
if abs(v_n2 - 1.607) > 0.1:
    print("Nozzle velocity error (2)")
else:
    pass

p_t = presure_loss_tube_finder(v_t, Re_t)
if abs(p_t - 2533.2) > 10:
    print("Tube pressure loss error")
else:
    pass

p_e = pressure_loss_ends_finder(v_t, 0.47, 0.83)
if abs(p_e - 2052.8) > 10:
    print("Ends pressure loss error")
else:
    pass

p_noz2 = pressure_loss_nozzle_finder(v_n2)
if abs(p_noz2 - 2556.9) > 10:
    print("Nozzle pressure loss error (2)")
else:
    pass

sigma = sigma_finder(N)
if abs(sigma - 0.0879) >0.05:
    print("Sigma error")
else:
    pass

A_sh = shell_area_finder(0.014, 9)
if abs(A_sh - 9.6e-4) > 0.0001:
    print("Shell area error")
else:
    pass

dshp = shell_chic_length_finder(A_sh)
if abs(dshp - 0.0191) > 0.001:
    print("Characteristic length error")
else:
    pass

v_sh = v_shell_finder(mdot1, A_sh)
if abs(v_sh - 0.526) > 0.01:
    print("Shell velocity error")
else:
    pass

Re_sh = reynolds_shell_finder(v_sh, dshp)
if abs(Re_sh - 15279) > 10:
    print("Shell reynolds error")
else:
    pass

p_sh = pressure_loss_shell_finder(Re_sh, N, v_sh, a)
if abs(p_sh - 516.6) > 2:
    print("Shell pressure loss error")
else:
    pass

v_n1 = v_nozzle_finder(mdot1)
if abs(v_n1 - 1.607) > 0.01:
    print("Nozzle velocity error 1")
else:
    pass

p_n1 = pressure_loss_nozzle_finder(v_n1)
if abs(p_n1 - 2556.9) > 5:
    print("Nozzle pressure loss error 1")
else:
    pass