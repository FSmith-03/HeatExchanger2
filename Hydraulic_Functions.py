import numpy as np

#function to find tube velocity
def v_tube_finder(mdot2, N):
    mdot_tube = mdot2/N
    d_i = 0.006
    v_tube = mdot_tube/(990.1*np.pi*0.25*(d_i)**2)
    return v_tube

#function to find tube reynolds number
def reynolds_tube_finder(v_tube):
    rho = 990.1
    mu = 6.51e-4
    d_i = 0.006
    reynolds_tube = (rho*v_tube*d_i)/mu
    return reynolds_tube

#function to find nozzle velocity
def v_nozzle_finder(mdot2):
    rho = 990.1
    d_n = 0.02
    v_nozzle = mdot2/(rho*np.pi*0.25*d_n**2)
    return v_nozzle

def presure_loss_tube_finder(v_tube, reynolds_tube):
    rho = 990.1
    d_i = 0.006
    L = 0.35
    f = (1.82*np.log10(reynolds_tube)-1.64)**-2
    pressure_loss_tube = 0.5*rho*(v_tube)**2*f*(L/d_i)
    return pressure_loss_tube

def sigma_finder(N):
    d_i = 0.006
    d_p = 0.064
    sigma = N*d_i**2/d_p**2
    return sigma

def kc_ke_finder(sigma_ref):
    kc = [0.46, 0.425, 0.39, 0.34, 0.3, 0.26, 0.21, 0.18, 0.14, 0.1]
    ke = [0.8, 0.63, 0.46, 0.33, 0.2, 0.1, 0.03, -0.03, -0.08, -0.1]
    sigma = np.arange(0.1, 1.1, 0.1)
    degree_kc = 5       # degree for polynomial regression
    if abs(sigma_ref)>1:
        print("Absolute value of sigma is greater than one, please use a value of sigma less than 1")
    else:
        # Create a Vandermonde matrix of the independent variable sigma
        Xkc = np.vander(sigma, degree_kc + 1, increasing=True)
        # Perform least squares polynomial regression
        coefficients, _, _, _ = np.linalg.lstsq(Xkc, kc, rcond=None)
        # Predict kc value for the reference_sigma
        reference_kc = np.polyval(coefficients[::-1], sigma_ref)

        # Create a Vandermonde matrix of the independent variable sigma
        Xke = np.vander(sigma, degree_kc + 1, increasing=True)
        # Perform least squares polynomial regression
        coefficients, _, _, _ = np.linalg.lstsq(Xke, ke, rcond=None)
        # Predict kc value for the reference_sigma
        reference_ke = np.polyval(coefficients[::-1], sigma_ref)

        return [reference_kc, reference_ke]

#print(kc_ke_finder(1.2))
def pressure_loss_ends_finder(v_tube, kc, ke):
    rho = 990.1
    pressure_loss_ends = 0.5*rho*(v_tube)**2*(kc+ke)
    return pressure_loss_ends


def shell_area_finder(Y, N_b):
    L = 0.35
    d_sh = 0.064
    d_o = 0.008
    B = L/(N_b+1)
    shell_area = (d_sh/Y)*(Y-d_o)*B
    if shell_area > 0:
        return shell_area
    else:
        return "shell_error"
    

def v_shell_finder(mdot1, shell_area):
    rho = 990.1
    v_shell = mdot1/(rho*shell_area)
    return v_shell

def shell_chic_length_finder(shell_area):
    d_sh = 0.064
    pipe_area = np.pi*0.25*d_sh**2
    shell_chic_length = d_sh*(shell_area/pipe_area)
    return shell_chic_length

def reynolds_shell_finder(v_shell, shell_chic_length):
    rho = 990.1
    mu = 6.51e-4
    reynolds_shell = (rho*v_shell*shell_chic_length)/mu
    return reynolds_shell

def pressure_loss_shell_finder(reynolds_shell, N, v_shell, a):
    rho = 990.1
    pressure_loss_shell = 4*a*(reynolds_shell**(-0.15))*N*rho*v_shell**2
    return pressure_loss_shell

def pressure_loss_nozzle_finder(v_nozzle):
    rho = 990.1
    pressure_loss_nozzle = rho*v_nozzle**2
    return pressure_loss_nozzle

#Note input to pressure checker needs to be in bars
def pressure_checker(pressure_loss, mdot, stream):
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
    elif stream == 2:
        # Create a Vandermonde matrix of the independent variable x
        X = np.vander(x_hot_data, degree + 1, increasing=True)
        # Perform least squares polynomial regression
        coefficients, _, _, _ = np.linalg.lstsq(X, y_hot_data, rcond=None)
        reference_pressure_rise = np.polyval(coefficients[::-1], mdot)
    else:
        print("Please use a value of 1 or 2 for the stream value.")
    
    #Check to ensure the pressure rise in the compressor is greater than the pressure drop in the HX

    if pressure_loss > reference_pressure_rise:
        return reference_pressure_rise, False
    else:
        return reference_pressure_rise, True

#print(pressure_checker(0.04360, 0.45, 2))

def H_finder(reynolds_tube, reynolds_shell):
    Pr = 4.31
    c = 0.15    #change for triangular/square pitch
    kw = 0.632
    k_tube = 386
    d_i = 0.006
    d_o = 0.008
    Nu_i = 0.023*(reynolds_tube)**0.8*4.31**0.3     #care for change in mass flow rate changing the reynolds number
    Nu_o = c*reynolds_shell**0.6*Pr**0.3
    h_i = (Nu_i*kw)/d_i
    h_o = (Nu_o*kw)/d_o
    Hinv = (1/(h_i))+((d_i*np.log(d_o/d_i))/(2*k_tube))+(1/(h_o))*(d_i/d_o)
    H = 1/Hinv
    return H

def temperature_solvernew1(mdot1, mdot2, H, A, F, N_i, error):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 57
    T1out = 21
    step = (T2in-T1in)/N_i
    timeout = 1000000
    counter = 0
    while T2out > T1in:
        check = True
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        if ((T2in-T1out)/(T2out-T1in)) == 1:
            check = False
            T2out = T2out - step
        if check:
            LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
            qdot1 = cp*mdot1*(T1out-T1in)
            qdot2 = cp*mdot2*(T2in-T2out)
            qdot3 = H*A*F*LMTD
            answer_new = [qdot1, qdot2, qdot3, T1out, T2out, LMTD]
            delta_qdot1_3_new = abs(answer_new[0] - answer_new[2])
            delta_qdot2_3_new = abs(answer_new[1] - answer_new[2])
            if delta_qdot1_3_new < error and delta_qdot2_3_new < error:
                return T1out, T2out, LMTD
            T2out = T2out - step
        counter += 1
        if counter == timeout:
            return -1000, -1000, -1000
    return -1000, -1000, -1000

def efficiency_finder(LMTD, H, A, F, mdot2):
    cp = 4179
    Qdot = LMTD*H*A*F
    epsilon = Qdot/(mdot2*cp*(60-20))
    return epsilon

#Print statements to test some functions
#print(temperature_solver1(0.5, 0.47, 4179, 9878.7, 0.0858, 1))
#print(temperature_solver2(0.5, 0.47, 4179, 9878.7, 0.0858, 1))
#print(efficiency_finder(28.2, 9878.7, 0.0858, 1, 0.47))
#print(H_finder(11283, 15281))
#print(presure_loss_tube_finder(v_tube_finder(0.45, 13), reynolds_tube_finder(v_tube_finder(0.45, 13)) ))