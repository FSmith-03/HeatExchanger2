import scipy
from ntu import *
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute

# Define the function to maximise (effectiveness)
def objective_function_NTU(params):
    # Extract parameters
    N, Nb, Y, mdot1, mdot2 = params
    a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
    cp = 4179
    F = 1.0



    shell_area_value = shell_area_finder(Y, N_b)
    v_shell_value = v_shell_finder(mdot1, shell_area_value)
    shell_chic_length_value = shell_chic_length_finder(shell_area_value)
    reynolds_shell_value = reynolds_shell_finder(v_shell_value, shell_chic_length_value)
    v_tube_value = v_tube_finder(mdot2, N)
    reynolds_tube_value = reynolds_tube_finder(v_tube_value)


    U_pipe = H_finder(reynolds_tube_value, reynolds_shell_value)
    A_pipe = N * d_i* L * np.pi

    

    # Hot stream
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
    
    
    # Check if both pressure rises are greater than the pressure drops
    if state1 and state2:
        #print(pressure_loss_1, pressure_rise_1)
        #print(pressure_loss_2, pressure_rise_2)
        return -effect_NTU_counterflow(mdot1, mdot2, U_pipe, A_pipe)    
    else:
        print(pressure_loss_1, pressure_rise_1)
        print(pressure_loss_2, pressure_rise_2)
        return 1e6
        

# Define the bounds for each parameter
bounds = [(1,10), (5,19), (0.005, 0.020), (0.5, 0.7), (0.3, 0.7)]

# Perform global optimization using differential evolution
result = differential_evolution(objective_function_NTU, bounds)

#Extract optimal parameters
optimal_params = result.x
optimal_effective = -result.fun

#result_brute = brute(objective_function_NTU, ranges=bounds, Ns=10, finish=None)
#min_value = objective_function_NTU(result_brute)

# Extract optimal parameters
#optimal_params = result_brute
#optimal_effective = -min_value


# Print results
print("Optimal Parameters:", optimal_params)
print("Optimal Effectiveness:", optimal_effective)

Qdotactual_counter = optimal_effective * Qdotmax
T2out_NTU_counter = T2in - Qdotactual_counter/(mdot2*cp)
T1out_NTU_counter = Qdotactual_counter/(mdot1*cp) + T1in

print(T2out_NTU_counter)
print(T1out_NTU_counter)