import numpy as np
import Hydraulic_Functions as hf
from Hydraulic_Functions import *
import itertools

T1in = 20
T2in = 60
cp = 4179
d_i = 0.006
#L = 0.35
d_sh = 0.064
#N = 10
a = 0.2

def c_min(m1, m2):
 return(min(m1*cp,m2*cp))

def c_max(m1, m2):
 return(max(m1*cp,m2*cp))

def RC(m1,m2):
 ratio = c_min(m1,m2)/c_max(m1,m2)
 return (ratio)

def NTU(U, A, m1, m2):
 c_hot = 4179
 c_cold = 4179
 NTU_out = U*A/c_min(m1,m2)
 return (NTU_out)

#U_pipe = 9878.7
#A_pipe = N * d_i* L * np.pi


def effect_NTU_counterflow(m1, m2, U, A): #for constant cp
    NTU_value = NTU(U, A, m1, m2)
    numer = (1 - np.exp(-NTU_value * (1 - RC(m1, m2))))
    denom = (1 - RC(m1, m2) * np.exp(-NTU_value * (1 - RC(m1, m2))))
    effectiveness = (1 - np.exp(-NTU_value * (1 - RC(m1, m2)))) / (1 - RC(m1, m2) * np.exp(-NTU_value * (1 - RC(m1, m2)))) #for RC is not equal to 1 i.e. different cps
    #effectiveness = NTU_value/(1+NTU_value) #pure counter for RC = 1
    #effectiveness = 2/(2+np.sqrt(2)*(1/np.tan(NTU_value/np.sqrt(2)))) # 1–2 shell-and-tube exchanger
    return (effectiveness)

def effect_NTU_parallelflow(m1, m2, U, A): #for constant cp
    NTU_value = NTU(U, A, m1, m2)
    effectiveness = (1- np.exp(-NTU_value(1+RC(m1,m2))))/(1+RC(m1,m2))
    #effectiveness = 0.5*(1-np.exp(-NTU_value)) #pure parallel rc=1
    return (effectiveness)

def effect_NTU_1_shell_2_pass (m1, m2, U, A): #effectivness of each shell pass
  NTU_value = NTU(U, A, m1, m2)
  RC_value = RC(m1, m2)
  Gamma = NTU_value * np.sqrt(1+RC_value**2)
  effectiveness = 2/((1+RC_value)+ np.sqrt(1+RC_value**2) * 1/(np.tan(Gamma/2)))
  return (effectiveness) 

def effect_NTU_shellandpass(m1, m2, U, A, passes): # effectivness N shell passes, 2N tube passes
  NTU_value = NTU(U, A, m1, m2)
  RC_value = RC(m1, m2)
  effect_1sp = effect_NTU_1_shell_2_pass (m1, m2, U, A)
  if (1 - effect_1sp) != 0:  # Check if the denominator is not zero
    effectiveness_multi = (((1 - effect_1sp * RC_value) / (1 - effect_1sp)) ** passes - 1) * (((1 - effect_1sp * RC_value) / (1 - effect_1sp)) ** passes - RC_value) ** (-1)
  else:
    effectiveness_multi = 0
  
  return (effectiveness_multi) 

"""
mdot1 = 0.5
mdot2 = 0.47
Y = 0.014 # Pitch for square arrangement
N_b = 9
N_passes = 2


shell_area_value = hf.shell_area_finder(Y, N_b)
v_shell_value = hf.v_shell_finder(mdot1, shell_area_value)
shell_chic_length_value = hf.shell_chic_length_finder(shell_area_value)
reynolds_shell_value = hf.reynolds_shell_finder(v_shell_value, shell_chic_length_value)
v_tube_value = hf.v_tube_finder(mdot2, N)
reynolds_tube_value = hf.reynolds_tube_finder(v_tube_value)


U_pipe = hf.H_finder(reynolds_tube_value, reynolds_shell_value)
A_pipe = N * d_i* L * np.pi
Qdotmax = mdot2*cp*(T2in-T1in)
#print(U_pipe, A_pipe)
Qdotactual_counter = effect_NTU_1_shell_2_pass(mdot1, mdot2, U_pipe, A_pipe) * Qdotmax


T2out_NTU_counter = T2in - Qdotactual_counter/(mdot2*cp)
T1out_NTU_counter = Qdotactual_counter/(mdot1*cp) + T1in

print(T2out_NTU_counter)
print(T1out_NTU_counter)

print(effect_NTU_counterflow(mdot1, mdot2, U_pipe, A_pipe))
print(effect_NTU_shellandpass(mdot1, mdot2, U_pipe, A_pipe, N_passes))
"""


copper_tube_mass = 0.20 #kg/m
acrylic_pipe_mass = 0.65 #kg/m
nozzle_mass = 0.025 #kg/nozzles
abs_sheet_mass = 2.39 #kg/m^2 Baffles & Splitters
abs_sheet_thickness = 1.5/1000 #m thick
photopolymer_resin_mass = 1150 #kg/m^2 Tubesheets and endplates
silicone_ring_mass = 0.8/1000 #kg/ring
o36_ring_mass = 5.3/1000 #kg/ring

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
         #print("mass over")
         #print(mass_tot)
         #print(copper_tube_mass_total, acrylic_pipe_mass_total, nozzle_mass_total,abs_sheet_mass_total,photopolymer_resin_mass_total,silicone_ring_mass_total,o36_ring_mass_total)
         print(mass_tot)
         return(False)

def tube_length(N, L, passes):
    copper_length = N*L*passes
    if copper_length <3.5:
        return(True)
    else:
        #print("length over")
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

    # Generate all possible combinations of variable values
    combinations = itertools.product(*[frange(min_val, max_val, step) for (min_val, max_val), step in zip(variable_ranges, step_sizes)])
    
    # Iterate over all combinations
    for combination in combinations:
        # Evaluate the objective function with the current combination of variable values
        value = objective_function(*combination)
        
        # Update maximum value and combination if necessary
        if value > max_value:
            #print(N_b)
            print("Replacing")
            
            max_value = value
            max_combination = combination
            
    return max_value, max_combination

def objective_function(L, N, Y, N_b, passes):
    mdot1, error1 = pressure_intersection_1(N, N_b, Y, L)
    mdot2, error2 = pressure_intersection_2(N, L)
    shell_area_value = shell_area_finder_NTU(Y, N_b, L)
    v_shell_value = hf.v_shell_finder(mdot1, shell_area_value)
    shell_chic_length_value = hf.shell_chic_length_finder(shell_area_value)
    reynolds_shell_value = hf.reynolds_shell_finder(v_shell_value, shell_chic_length_value)
    v_tube_value = hf.v_tube_finder(mdot2, N)
    reynolds_tube_value = hf.reynolds_tube_finder(v_tube_value)
    V_nozzle_value = v_nozzle_finder(mdot2)


    

    #mass check
    mass_under = total_mass(N, L, 4, N_b, np.pi * d_sh /4, passes, np.pi * d_sh /4 *1.5/1000, 8, 8 )

    U_pipe = hf.H_finder(reynolds_tube_value, reynolds_shell_value)
    A_pipe = N * d_i* L * np.pi * passes
    #print(NTU(U_pipe, A_pipe, mdot1, mdot2))

    #check copper length
    tube_length_value = tube_length(N, L, passes)

    if mass_under == True and tube_length_value == True:
      print(U_pipe) 
      return effect_NTU_shellandpass(mdot1, mdot2, U_pipe, A_pipe, passes)
    else:
       return 0
    

variable_ranges = [(0.25, 0.35), (5, 6), (0.012, 0.013), (1, 7), (1, 3)]  # Example variable ranges
step_sizes = [0.001, 1, 0.001, 1, 1]

max_value, max_combination = brute_force_maximizer(objective_function, variable_ranges, step_sizes)
print("Maximum value:", max_value)
print("Maximizing combination of variables:", max_combination)







