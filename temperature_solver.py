from Hydraulic_Functions import *
import numpy as np
from scipy.optimize import differential_evolution, brute, fsolve



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
            print(qdot1, qdot2, qdot3)
            answer_new = [qdot1, qdot2, qdot3, T1out, T2out, LMTD]
            delta_qdot1_3_new = abs(answer_new[0] - answer_new[2])
            delta_qdot2_3_new = abs(answer_new[1] - answer_new[2])
            print(delta_qdot1_3_new, delta_qdot2_3_new)
            if delta_qdot1_3_new < error and delta_qdot2_3_new < error:
                return T1out, T2out, LMTD
            T2out = T2out - step
        counter += 1
        if counter == timeout:
            return 1000, 1000, 1000
    return 1000, 1000, 1000
A = 13*np.pi*0.006*0.35
print(temperature_solvernew1(0.1, 0.1, 1882.4979991622515, 0.02638937829015426, 1, 2500, 10))