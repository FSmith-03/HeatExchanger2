import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute
from temperature_solvers_old import *
# Define the function to minimize (negative LMTD)
def objective_function(params):
    # Extract parameters
    N, Nb, Y, mdot1, mdot2 = params
    a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
    cp = 4179
    F = 1
    A = N*np.pi*0.006*0.35

    # Calculate LMTD (example function)
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
    
    # Thermal Analysis
    H = H_finder(reynolds_tube, reynolds_shell)
    T1out, T2out, LMTD = temperature_solver2(mdot1, mdot2, H, A, F)
    # Check if both pressure rises are greater than the pressure drops
    if state1 and state2:
        # Negative LMTD to convert maximization problem to minimization
        #print(pressure_loss_1, pressure_rise_1)
        #print(params)
        return -LMTD
    else:
        # Return a large value to discourage selecting this solution
        print(params)
        return 1e6

# Define the bounds for each parameter
bounds = [(1,10), (5,19), (0.005, 0.020), (0.1, 1), (0.1, 1)]

# Perform global optimization using differential evolution
result = differential_evolution(objective_function, bounds)

# Extract optimal parameters
optimal_params = result.x
optimal_LMTD = -result.fun  # Convert back to positive LMTD

#result_brute = brute(objective_function, ranges=bounds, Ns=10, finish=None)
#min_value = objective_function(result_brute)

# Extract optimal parameters
#optimal_params = result_brute
#optimal_LMTD = min_value

# Print results
print("Optimal Parameters:", optimal_params)
print("Optimal LMTD:", optimal_LMTD)
