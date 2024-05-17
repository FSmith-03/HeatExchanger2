import numpy as np

def temperature_solvernew(mdot1, mdot2, H, A, F, N_i, error):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59
    T1out = 21
    step = (T2in-T1in)/N_i
    answer_list = []
    while T2out > T1in:
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        qdot1 = cp*mdot1*(T1out-T1in)
        qdot2 = cp*mdot2*(T2in-T2out)
        qdot3 = H*A*F*LMTD
        answer_new = [qdot1, qdot2, qdot3, T1out, T2out, LMTD]
        delta_qdot1_3_new = abs(answer_new[0] - answer_new[2])
        delta_qdot2_3_new = abs(answer_new[1] - answer_new[2])
        if delta_qdot1_3_new < error and delta_qdot2_3_new < error:
            return T1out, T2out, LMTD
            break
        T2out -= step

def temperature_solver1(mdot1, mdot2, H, A, F):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59
    T1out = 21
    counter = 0
    N = 250
    for n in range(N):
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        T2outnew = T2in - (H*A*F/(mdot2*cp))*LMTD
        if abs(T2outnew-T2out) < 0.1:
            return [T1out, T2out, LMTD]
        T2out = T2outnew
        #print(T2out)
        counter += 1
        if counter == N-1:
            print("No solution found")

def temperature_solver1GPT(mdot1, mdot2, H, A, F):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59
    T1out = 21
    counter = 0
    N = 101
    for n in range(N):
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        T2outnew = T2in - (H*A*F/(mdot2*cp))*LMTD
        if abs(T2outnew-T2out) < 0.01:
            return [T1out, T2out, LMTD]
        T2out = T2outnew
        #print(T2out)
        counter += 1
        if counter == N-1:
            print("No solution found")
            return [None, None, None]  # Return default values indicating no solution found
        
def temperature_solver2(mdot1, mdot2, H, A, F):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59      #Starting values for iteration
    T1out = 21
    counter = 0
    N = 10000
    step = 0.01
    for n in range(N):
        T1out = T1in + (mdot2/mdot1)*(T2in - T2out)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        T2outnew = T2in - (H*A*F/(mdot2*cp))*LMTD
        #Iteration by reducing the T2out estimate sequentially by a step until the estimate reaches a minimum
        if T2outnew > T2out:
            return [T1out, T2out, LMTD]
        T2out = T2out - step
        #print(T2out)
        counter += 1
        if counter == N-1:
            print("No solution found increase number of iterations")

def temperature_solver3(mdot1, mdot2, H, A, F):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59      #Starting values for iteration
    T1out = 21
    counter = 0
    step = 0.1
    N = (T2in-T1in)/step
    LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
    qdot1 = cp*mdot1*(T1out-T1in)
    qdot2 = cp*mdot2*(T2in-T2out)
    qdotl = H*A*LMTD*F
    while abs(qdotl-qdot1) > 0.1:
        T1out += step
        qdot1 = cp*mdot1*(T1out-T1in)
        T2out = T2in-(mdot1/mdot2)*(T1out-T1in)
        LMTD = ((T2in - T1out)-(T2out-T1in))/(np.log((T2in-T1out)/(T2out-T1in)))
        qdotl =  H*A*LMTD*F
    return T1out, T2out, LMTD


def temperature_solver3GPT(mdot1, mdot2, H, A, F):
    T1in = 20
    T2in = 60
    cp = 4179
    T2out = 59      # Starting values for iteration
    T1out = 21
    counter = 0
    step = 0.1
    max_iterations = 1000  # Setting a maximum number of iterations
    N = (T2in-T1in)/step
    while counter < max_iterations:
        counter += 1
        T1out_new = T1out + step
        qdot1 = cp * mdot1 * (T1out - T1in)
        T2out = T2in - (mdot1 / mdot2) * (T1out - T1in)
        LMTD = ((T2in - T1out) - (T2out - T1in)) / np.log((T2in - T1out) / (T2out - T1in))
        qdotl = H * A * LMTD * F
        if abs(qdotl - qdot1) <= 0.1:
            break  # Exiting the loop if convergence is achieved
        T1out = T1out_new  # Update T1out for the next iteration

    if counter == max_iterations:
        print("Max iterations reached without convergence.")

    return T1out, T2out, LMTD