import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt

A = 0.007896
q = 0.01
d1 = 0.24203
d2 = 0.55586
d3 = 0.07849
mu1 = 0.4554
mu2 = 1.21382
mu3 = 0.1325
delta = 0.0000213
delta1 = 0.196
alpha = 7.98025247e-8
beta = 3.99332316e-18
delta2 = 9.99962040e-01
p1 = 0.85

def SEIQCR(t, r):
    S, E, I, Q, C, R = r
    fs = A*S - (delta*S) - (beta*S*I) - (alpha*C*S) - q*S
    fe = (beta*I*S) - (delta*E) - (delta1*E) + (alpha*C*S)
    fi = delta1*E - (delta*I) - (delta2*I) - (mu1*I) - (d1*I)
    fq = (p1*delta2*I) - (delta*Q) - (mu2*Q) - (d2*Q)
    fc = (delta2*I*(1-p1)) - (delta*C) - (mu3*C) - (d3*C)
    fr = mu1*I + mu2*Q + mu3*C - delta*R
    return fs, fe, fi, fq, fc, fr

# R = np.linspace(0.1,3,500)
# R = R.tolist();
#
# E = []
# E1 = []
# I = []
# I1 = []
# R1 = []
# for r in R:
#     E.append(0)
#     I.append(0)
#     E1.append((r-1)*(d+q)/beta)
#     I1.append(epsilon*(d+q)*(r-1)/((d+d1+delta)*beta))
# print((E1[499]-E1[0])/2.9)
# print((I1[499]-I1[0])/2.9)
# plt.plot(R, E1, label = "Stable")
# plt.plot(R, E, label = "Unstable")
# plt.xlabel("Reproduction Rate")
# plt.ylabel("Exposed Individuals")
# plt.title("Transcritical Bifurcation for Exposed Individuals")
# plt.grid()
# plt.legend()
# plt.show()

# reading data
df = pd.read_excel("Italy Variable Data.xlsx",header = 0)
S = list(df['S'])
E = list(df['E'])
I = list(df['I'])
Q = list(df['Q'])
C = df['C']
R = list(df['R'])
days = list(df['days'])
C.rolling(50).sum()
C = list(C)
sol = solve_ivp(SEIQCR, (0, 365), (60461803,2703,94, 127, 2948,1), method = 'BDF')
plt.plot(sol.t, sol.y[0], label = "Susceptible")
plt.plot(days, S, label = "Susceptible")
# plt.plot(sol.t, sol.y[1], label = "Exposed")
# plt.plot(days, E, label = "Exposed")
# plt.plot(sol.t, sol.y[2], label = "Infected")
# plt.plot(days, I, label = "Infected")
# plt.plot(sol.t, sol.y[3], label = "Quarantined")
# plt.plot(days, Q, label = "Quarantined")
# plt.plot(sol.t, sol.y[4], label = "Carrier")
# plt.plot(days, C, label = 'Carrier')
# plt.plot(sol.t, sol.y[5], label = "Recovered")
# plt.plot(days, R, label = "Recovered")
# print(sol.y[5])
# plt.ylabel("Population")
# plt.xlabel("Time(days)")
# plt.title("Italy COVID-19 rate when q = 0.01")
plt.grid()
plt.legend()
plt.show()
