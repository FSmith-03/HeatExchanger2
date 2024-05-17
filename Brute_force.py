import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution, brute
import matplotlib.pyplot as plt

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

bounds = [(1,10), (5,19), (0.012, 0.020)] 
print(brute_force(bounds, 5, graphing=True))
#print(shell_area_finder(0.005, 5))


