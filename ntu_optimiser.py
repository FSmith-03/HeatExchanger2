import scipy
from ntu import *
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution

# Define the function to maximise (effectiveness)
def objective_function(params):
    # Extract parameters
    N, Nb, Y, mdot1, mdot2 = params
    a = 0.2         #a=0.2 for triangular pitch and 0.34 for square pitch
    cp = 4179
    F = 1.0

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
    pressure_rise_1, state1 = pressure_checker(pressure_loss_1, mdot2, 2)

    # Thermal Analysis
    H = H_finder(reynolds_tube, reynolds_shell)
    
    
    # Check if both pressure rises are greater than the pressure drops
    if state1 and state2:
        return effect_NTU_counterflow
    

# Define the bounds for each parameter
bounds = [(8,20), (5,15), (0.005, 0.020), (0.3, 0.6), (0.3, 0.6)]

# Perform global optimization using differential evolution
result = differential_evolution(objective_function, bounds)

# Extract optimal parameters
optimal_params = result.x
optimal_effective = -result.fun  # Convert back to positive LMTD

# Print results
print("Optimal Parameters:", optimal_params)
print("Optimal LMTD:", optimal_LMTD)