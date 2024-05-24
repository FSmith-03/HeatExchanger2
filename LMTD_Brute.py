import numpy as np
import Hydraulic_Functions as hf
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute
import matplotlib.pyplot as plt
import itertools

copper_tube_mass = 0.20 #kg/m
acrylic_pipe_mass = 0.65 #kg/m
nozzle_mass = 0.025 #kg/nozzles
abs_sheet_mass = 2.39 #kg/m^2 Baffles & Splitters
abs_sheet_thickness = 1.5/1000 #m thick
photopolymer_resin_mass = 1150 #kg/m^2 Tubesheets and endplates
silicone_ring_mass = 0.8/1000 #kg/ring
o36_ring_mass = 5.3/1000 #kg/ring

T1in = 20
T2in = 60
cp = 4179
d_i = 0.006
#L = 0.35
d_sh = 0.064
#N = 10
a = 0.2

def total_mass(N, L, nozzles, baffles, baffle_area, passes, end_plate_vol, sil_rings, o36_rings):
      copper_tube_mass_total = copper_tube_mass * N * L * passes
      acrylic_pipe_mass_total = acrylic_pipe_mass * L
      nozzle_mass_total = nozzle_mass * nozzles
      abs_sheet_mass_total = (abs_sheet_mass * baffles * baffle_area * abs_sheet_thickness) + (passes-1)*L*(d_sh)*abs_sheet_thickness #check for splitters
      photopolymer_resin_mass_total = photopolymer_resin_mass * end_plate_vol
      silicone_ring_mass_total = silicone_ring_mass * sil_rings
      o36_ring_mass_total = o36_ring_mass * o36_rings
      mass_tot = copper_tube_mass_total + acrylic_pipe_mass_total + nozzle_mass_total + abs_sheet_mass_total + photopolymer_resin_mass_total + silicone_ring_mass_total + o36_ring_mass_total 
      if mass_tot <1.2:
         return(True)
      else:
         print("mass over")
         #print(mass_tot)
         print(copper_tube_mass_total, acrylic_pipe_mass_total, nozzle_mass_total,abs_sheet_mass_total,photopolymer_resin_mass_total,silicone_ring_mass_total,o36_ring_mass_total)
         return(False)

def tube_length(N, L, passes):
    copper_length = N*L*passes
    if copper_length <3.5:
        return(True)
    else:
        print("length over")
        return(False)

def frange(start, stop, step):
    while start < stop:
        yield round(start, 10)  # rounding to handle floating-point precision
        start += step


def brute_force_maximizer(objective_function, variable_ranges, step_sizes):
    max_value = float('-inf')
    max_combination = None
    print("Variable Ranges:", variable_ranges)
    print("Step Sizes:", step_sizes)
    LMTD_List = []
    combination_list = []
    # Generate all possible combinations of variable values
    combinations = itertools.product(*[frange(min_val, max_val, step) for (min_val, max_val), step in zip(variable_ranges, step_sizes)])
    
    # Iterate over all combinations
    for combination in combinations:
        # Evaluate the objective function with the current combination of variable values
        value = objective_function(*combination)
        
        # Update maximum value and combination if necessary
        LMTD_List.append(value)
        if value > max_value:
            print("Replacing")
            combination_list.append([value ,combination])
            max_value = value
            max_combination = combination
            
    return max_value, max_combination, LMTD_List, combination_list

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
    H = H_finder(reynolds_tube_value, reynolds_shell_value, L, N_b, Y, mdot1, passes, N)
    #print(mdot1, mdot2, H, A, F)
    #T1out, T2out, LMTD = temperature_solvernew2(mdot1, mdot2, H, A, F, 25000, 100)
    T1out, T2out, qdot = temperature_solvernew2(mdot1, mdot2, H, A, F, 1000, passes)


    

    #mass check
    mass_under = total_mass(N, L, 4, N_b, np.pi * d_sh /4, passes, np.pi * d_sh /4 *1.5/1000, 8, 8 )

    #check copper length
    tube_length_value = tube_length(N, L, passes)

    if mass_under == True and tube_length_value == True:
      return qdot
    else:
       return 0
    

variable_ranges = [(0.21, 0.35), (4, 8), (0.012, 0.016), (3, 7), (4, 5)]  # Example variable ranges
step_sizes = [0.01, 1, 0.001, 1, 1]

max_value, max_combination, LMTD_List, combination_list = brute_force_maximizer(objective_function, variable_ranges, step_sizes)
print("Maximum Qdot value:", max_value)
print("Maximizing combination of variables:", max_combination)
print("Next best sets of combinations:", combination_list[-5:])
x_List = np.arange(1, len(LMTD_List)+1, 1)

plt.scatter(x_List, LMTD_List, s=5, color = 'cyan', label = 'LMTD Trend')
plt.xlabel('Iteration Number')
plt.ylabel('LMTD')
plt.title('Improvement in LMTD with iteration count')
plt.legend(fontsize = 'large')
plt.show()