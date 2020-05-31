import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import pandas as pd

df = pd.read_csv("india.csv",header = 0)
S = list(df['S'])
E = list(df['E'])
I = list(df['I'])
Q = list(df['Q'])
C = list(df['C'])
R = list(df['R'])
days = list(df['days'])
ns = len(days) #no of time points

def SEIQCR(t, r, p):
    #parameters
    A = 0.0182
    q = 0.01
    d1 = 0.0453
    d2 = 0.06744
    d3 = 0.00227
    mu1 = 0.203488
    mu2 = 0.3027
    mu3 = 0.0125
    delta = 0.0073
    delta1 = 0.196
    p1 = 0.85

    S, E, I, Q, C, R = r
    alpha, beta, delta2 = p

    #equations
    fs = A*S - (delta*S) - (beta*S*I) - (alpha*C*S) - q*S
    fe = (beta*I*S) - (delta*E) - (delta1*E) + (alpha*C*S)
    fi = delta1*E - (delta*I) - (delta2*I) - (mu1*I) - (d1*I)
    fq = (p1*delta2*I) - (delta*Q) - (mu2*Q) - (d2*Q)
    fc = (delta2*I*(1-p1)) - (delta*C) - (mu3*C) - (d3*C)
    fr = mu1*I + mu2*Q + mu3*C - delta*R
    return fs, fe, fi, fq, fc, fr

def simulate(p):
    est = np.zeros((6, ns))
    for i in range(ns):
        sol = solve_ivp(SEIQCR, (0, i), (60461803,2703,94, 127, 2948,1), method = 'LSODA', args = [p])
        est[0][i] = sol.y[0][-1]
        est[1][i] = sol.y[1][-1]
        est[2][i] = sol.y[2][-1]
        est[3][i] = sol.y[3][-1]
        est[4][i] = sol.y[4][-1]
        est[5][i] = sol.y[5][-1]
    return est

def objective(p):
    est = simulate(p)
    obj = 0.0
    for i in range(ns):
        obj += (((S[i]-est[0][i])/est[0][i])**2) + (((E[i]-est[1][i])/est[1][i])**2) + (((I[i]-est[2][i])/est[2][i])**2) + (((Q[i]-est[3][i])/est[3][i])**2) + (((C[i]-est[4][i])/est[4][i])**2) + (((R[i]-est[5][i])/est[5][i])**2)
    return obj

#parameter initial guess`
alpha = 5e-8
beta = 8e-8
delta2 = 0.7
p = [alpha, beta, delta2]

# sol = solve_ivp(SEIQCR, (0, 2), (60461803,2703,94, 127, 2948,1), method = 'LSODA',args = [p])
# print(sol.y)

#show objective
# print('Objective function is ' + str(objective(p)) + '\n');
#minimize
bnds = ((0,1),(0,1),(0,1))
solution = minimize(objective, p,bounds = bnds, options={'disp': True})
print(solution.x)
