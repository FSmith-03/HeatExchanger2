import numpy as np
import Hydraulic_Functions as hf
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute
import matplotlib.pyplot as plt
import itertools

def brute_force(bounds, step_i, graphing = False):
     # Extract initial parameters
    Nb_0 = bounds[1][0]
    Y_0 = bounds[2][0]
    Nb_1 = bounds[1][1]
    Y_1 = bounds[2][1]

    #Extract iteration step size
    Y_step = (Y_1-Y_0)/step_i
    #Extract starting values:
    Nb = Nb_0
    Y = Y_0
    a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
    cp = 4179
    F = 1
    N = 4
    counter = 0
    timeout = 1e10
    output = [0, 0, 0]
    design_params = [100, 100, 100, 100, 100]
    x = []
    y = []
    x_opt = []
    y_opt = []
    params_opt = []
    while N < 11:
        while Nb < Nb_1:
            while Y < Y_1:
                A = N*np.pi*0.006*0.35
                mdot2, error2 = pressure_intersection_2(N, mdot_upper=0.8)
                v_tube = v_tube_finder(mdot2, N)
                reynolds_tube = reynolds_tube_finder(v_tube)

                # Cold stream
                #Check for shell error from invalid geometries
                A_sh = shell_area_finder(Y, Nb)
                mdot1, error1 = pressure_intersection_1(N, Nb, Y, mdot_upper=0.8)
                
                v_shell = v_shell_finder(mdot1, A_sh)
                v_nozzle_1 = v_nozzle_finder(mdot1)
                d_sh = shell_chic_length_finder(A_sh)
                reynolds_shell = reynolds_shell_finder(v_shell, d_sh)

                # Thermal Analysis
                H = H_finder(reynolds_tube, reynolds_shell)
                #print(mdot1, mdot2, H, A, F)
                T1out, T2out, LMTD = temperature_solvernew1(mdot1, mdot2, H, A, F, 25000, 100)
                print(T1out, T2out, LMTD)
                #print(mdot2)
                #Append the LMTD only if there is no pressure issue and the LMTD is better than the previous estimate
                if LMTD > output[2]:
                    params = [N, Nb, Y, mdot1, mdot2]
                    #print(params)
                    output = [T1out, T2out, LMTD]
                    point = [counter, LMTD]
                    x_opt.append(counter)
                    y_opt.append(LMTD)
                    params_opt.append(params)
                    counter +=1
                #Temperature solver throws out -1000 if there is an issue with the input values.
                check = LMTD == -1000
                if check == False:
                    x.append(counter)
                    y.append(LMTD)
                print(counter)
                if counter == timeout:
                    raise TimeoutError
                Y += Y_step
            Y = Y_0
            Nb += 1
        Y = Y_0
        Nb = Nb_0
        N += 1 
    print("Optimal design parameters are:", params)
    print("Optimal temperatures are:", output)
    if graphing == True:
        plt.scatter(x, y, s=10, color = 'cyan', label = 'LMTD Trend')
        plt.plot(x_opt, y_opt, color = 'red', label = 'LMTD Optimal')
        plt.xlabel('Iteration Number')
        plt.ylabel('LMTD')
        plt.title('Improvement in LMTD with iteration count')
        plt.legend(fontsize = 'large')
        plt.show()

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
            
            max_value = value
            max_combination = combination
            
    return max_value, max_combination, LMTD_List

def objective_function(L, N, Y, N_b, passes):
    A = N*np.pi*0.006*L
    F = 1
    mdot2, error2 = pressure_intersection_2(N, mdot_upper=0.8)
    v_tube = v_tube_finder(mdot2, N)
    reynolds_tube_value = reynolds_tube_finder(v_tube)

    # Cold stream
    #Check for shell error from invalid geometries
    A_sh = shell_area_finder(Y, N_b)
    mdot1, error1 = pressure_intersection_1(N, N_b, Y, mdot_upper=0.8)
    
    v_shell = v_shell_finder(mdot1, A_sh)
    v_nozzle_1 = v_nozzle_finder(mdot1)
    d_sh = shell_chic_length_finder(A_sh)
    reynolds_shell_value = reynolds_shell_finder(v_shell, d_sh)

    # Thermal Analysis
    H = H_finder(reynolds_tube_value, reynolds_shell_value)
    #print(mdot1, mdot2, H, A, F)
    T1out, T2out, LMTD = temperature_solvernew1(mdot1, mdot2, H, A, F, 25000, 100)


    

    #mass check
    mass_under = total_mass(N, L, 4, N_b, np.pi * d_sh /4, passes, np.pi * d_sh /4 *1.5/1000, 8, 8 )

    U_pipe = hf.H_finder(reynolds_tube_value, reynolds_shell_value)
    A_pipe = N * d_i* L * np.pi
    #print(NTU(U_pipe, A_pipe, mdot1, mdot2))

    #check copper length
    tube_length_value = tube_length(N, L, passes)

    if mass_under == True and tube_length_value == True:
      return LMTD
    else:
       return 0
    

variable_ranges = [(0.15, 0.25), (2, 13), (0.012, 0.013), (1, 10), (1, 3)]  # Example variable ranges
step_sizes = [0.01, 1, 0.001, 1, 1]

max_value, max_combination, LMTD_List = brute_force_maximizer(objective_function, variable_ranges, step_sizes)
print("Maximum value:", max_value)
print("Maximizing combination of variables:", max_combination)

x_List = np.arange(1, len(LMTD_List)+1, 1)

plt.scatter(x_List, LMTD_List, s=5, color = 'cyan', label = 'LMTD Trend')
plt.xlabel('Iteration Number')
plt.ylabel('LMTD')
plt.title('Improvement in LMTD with iteration count')
plt.legend(fontsize = 'large')
plt.show()