import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt

A = 0.0182
q = 0.01
d1 = 0.0453
d2 = 0.06744
d3 = 0.00227
mu1 = 0.203488
mu2 = 0.3027
mu3 = 0.0125
delta = 0.073
delta1 = 0.196
p1 = 0.85
alpha = 2.81534153e-8
beta = 3.21361552e-18
delta2 = 7.00000000e-01

def SEIQCR(t, r):
    S, E, I, Q, C, R = r
    fs = A*S - (delta*S) - (beta*S*I) - (alpha*C*S) - q*S
    fe = (beta*I*S) - (delta*E) - (delta1*E) + (alpha*C*S)
    fi = delta1*E - (delta*I) - (delta2*I) - (mu1*I) - (d1*I)
    fq = (p1*delta2*I) - (delta*Q) - (mu2*Q) - (d2*Q)
    fc = (delta2*I*(1-p1)) - (delta*C) - (mu3*C) - (d3*C)
    fr = mu1*I + mu2*Q + mu3*C - delta*R
    return fs, fe, fi, fq, fc, fr



# reading data
df = pd.read_csv("india.csv",header = 0)
S = list(df['S'])
E = list(df['E'])
I = list(df['I'])
Q = list(df['Q'])
C = df['C']
R = list(df['R'])
days = list(df['days'])
sol = solve_ivp(SEIQCR, (0, 365), (1379985689, 9056, 93, 63, 9470, 14), method = 'BDF')
# plt.plot(sol.t, sol.y[0], label = "Susceptible")
# plt.plot(days, S, label = "Susceptible")
# plt.plot(sol.t, sol.y[1], label = "Exposed")
# plt.plot(days, E, label = "Exposed")
# plt.plot(sol.t, sol.y[2], label = "Infected")
# plt.plot(days, I, label = "Infected")
# plt.plot(sol.t, sol.y[3], label = "Quarantined")
# plt.plot(days, Q, label = "Quarantined")
# plt.plot(sol.t, sol.y[4], label = "Carrier")
# plt.plot(days, C, label = 'Carrier')
plt.plot(sol.t, sol.y[5], label = "Recovered")
plt.plot(days, R, label = "Recovered")
# print(sol.y[5])
# plt.ylabel("Population")
# plt.xlabel("Time(days)")
# plt.title("Italy COVID-19 rate when q = 0.01")
plt.grid()
plt.legend()
plt.show()
