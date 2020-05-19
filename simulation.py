import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
A = 0.007896
q = 0.01
mu1 = 0.4554
mu2 = 1.21382
mu3 = 0.1325
d1 = 0.24203
d2 = 0.55586
d3 = 0.07849
delta = 0.000213
delta1 = 0.196

def SCIR(t, r):
    S, E, I, Q, C, R = r
    fs = A - (delta*S) - (beta*S*I) - (alpha*C*S) - q*S
    fe = (beta*I*S) - (delta*E) - (delta1*E) + (p2*alpha*C*S)
    fi = delta1*E - (delta*I) - (delta2*I) - (mu1*I) - (d1*I)
    fq = (p1*delta2*I) - (delta*Q) - (mu2*Q) - (d2*Q)
    fc = (delta2*I*(1-p1)) - (delta*C) + ((1-p2)*alpha*C*S) - (mu3*C) - (d3*C)
    fr = mu1*I + mu2*Q + mu3*C - delta*R
    return fs, fe, fi, fq, fc, fr

R = np.linspace(0.1,3,500)
R = R.tolist();

E = []
E1 = []
I = []
I1 = []
R1 = []
for r in R:
    E.append(0)
    I.append(0)
    E1.append((r-1)*(d+q)/beta)
    I1.append(epsilon*(d+q)*(r-1)/((d+d1+delta)*beta))
print((E1[499]-E1[0])/2.9)
print((I1[499]-I1[0])/2.9)
# plt.plot(R, E1, label = "Stable")
# plt.plot(R, E, label = "Unstable")
# plt.xlabel("Reproduction Rate")
# plt.ylabel("Exposed Individuals")
# plt.title("Transcritical Bifurcation for Exposed Individuals")
# plt.grid()
# plt.legend()
# plt.show()


# R = (A*beta)/((d + q)*(d + epsilon))
# print(R)
# sol = solve_ivp(SCIR, (0, 365), (136956,10,1.8985,0.3260), method='DOP853', dense_output = 'True')
# plt.plot(sol.t, sol.y[0], label = "Susceptible")
# plt.plot(sol.t, sol.y[1], label = "Exposed")
# plt.plot(sol.t, sol.y[2], label = "Infected")
# plt.plot(sol.t, sol.y[3], label = "Recovered")
# plt.ylabel("Population (* 10000)")
# plt.xlabel("Time(days)")
# plt.title("India COVID-19 rate when q = 0.001")
# plt.grid()
# plt.legend()
# plt.show()
