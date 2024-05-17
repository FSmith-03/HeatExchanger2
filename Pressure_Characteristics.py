from Hydraulic_Functions import *
import numpy as np
import matplotlib.pyplot as plt
N = 6
Y = 0.014
Nb = 6
def pressure_loss_2_plot(mdot2):
    A = N*np.pi*0.006*0.35

    #Hot stream
    v_tube = v_tube_finder(mdot2, N)
    reynolds_tube = reynolds_tube_finder(v_tube)
    v_nozzle_2 = v_nozzle_finder(mdot2)
    pressure_loss_tube = presure_loss_tube_finder(v_tube, reynolds_tube)
    sigma = sigma_finder(N)
    #print(sigma)
    kc, ke = kc_ke_finder(sigma)
    pressure_loss_ends = pressure_loss_ends_finder(v_tube, kc, ke)
    pressure_loss_nozzle_2 = pressure_loss_nozzle_finder(mdot2)

    pressure_loss_2 = (pressure_loss_ends + pressure_loss_nozzle_2 + pressure_loss_tube)/1e5
    return pressure_loss_2

def pressure_loss_1_plot(mdot1):
    #Cold stream
    A_sh = shell_area_finder(Y, Nb)
    v_shell = v_shell_finder(mdot1, A_sh)
    v_nozzle_1 = v_nozzle_finder(mdot1)
    d_sh = shell_chic_length_finder(A_sh)
    reynolds_shell = reynolds_shell_finder(v_shell, d_sh)
    pressure_loss_shell = pressure_loss_shell_finder(reynolds_shell, N, v_shell, 0.2)
    pressure_loss_nozzle_1 = pressure_loss_nozzle_finder(v_nozzle_1)
    pressure_loss_1 = (pressure_loss_shell + pressure_loss_nozzle_1)/1e5
    return pressure_loss_1

def pressure_rise(mdot, stream):
    x_cold_data = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
    y_cold_data = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776])
    x_hot_data = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694])
    y_hot_data = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615])
    degree = 5
    if stream == 1:
        # Create a Vandermonde matrix of the independent variable x
        X = np.vander(x_cold_data, degree + 1, increasing=True)
        # Perform least squares polynomial regression
        coefficients, _, _, _ = np.linalg.lstsq(X, y_cold_data, rcond=None)
        reference_pressure_rise = np.polyval(coefficients[::-1], mdot)
        return reference_pressure_rise
    elif stream == 2:
        # Create a Vandermonde matrix of the independent variable x
        X = np.vander(x_hot_data, degree + 1, increasing=True)
        # Perform least squares polynomial regression
        coefficients, _, _, _ = np.linalg.lstsq(X, y_hot_data, rcond=None)
        reference_pressure_rise = np.polyval(coefficients[::-1], mdot)
        return reference_pressure_rise
    else:
        print("Please use a value of 1 or 2 for the stream value.")

mdot1 = np.arange(0.2, 0.8, 0.05)
mdot2 = np.arange(0.2, 0.6, 0.05)
y_loss_1 = pressure_loss_1_plot(mdot1)
y_loss_2 = pressure_loss_2_plot(mdot2)
y_rise_1 = pressure_rise(mdot1, 1)
y_rise_2 = pressure_rise(mdot2, 2)

plt.plot(mdot1, y_loss_1, color = 'cyan', label = 'Cold Pressure Loss')
plt.plot(mdot2, y_loss_2, color = 'red', label = 'Hot Pressure Loss')
plt.plot(mdot1, y_rise_1, color = 'cyan', label = 'Cold Pressure Rise')
plt.plot(mdot2, y_rise_2, color = 'red', label = 'Hot Pressure Rise')
plt.xlabel('Mdot')
plt.ylabel('Pressure Loss')
plt.title('Improvement in LMTD with iteration count')
plt.legend(fontsize = 'large')
plt.show()