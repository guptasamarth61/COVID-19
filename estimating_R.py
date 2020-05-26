import numpy as np
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

d1_star = delta+delta2+mu1+d1
d2_star = delta+mu2+d2
d3_star = delta+mu3+d3
R = np.linspace(0.1,3,1000)
R = R.tolist()
E = []
I = []
Q = []
C = []
S = []
E1 = []
I1 = []
C1 = []
S1 = []
Q1 = []
for r in R:
    E.append(d1_star*A*delta1*(r-1)/((delta+delta1)*r*delta1))
    E1.append(0)
    Q1.append(0)
    I1.append(0)
    C1.append(0)
    S1.append(A/(delta+q))
    Q.append(p1*delta2*(r-1)*A*delta1/(d2_star*r*(delta+delta1)))
    I.append(A*delta1*(r-1)/((delta+delta1)*r))
    C.append((1-p1)*delta2*(r-1)*A*delta1/(r*d3_star*(delta+delta1)))
    S.append(A/(r*(delta+q)))

# plt.plot(R, E, label = 'Exposed endemic eq')
# plt.plot(R, E1, label = 'Exposed disease-free eq')
# plt.plot(R, S, label = 'Susceptible endemic eq')
# plt.plot(R, S1, label = 'Susceptible disease-free eq')
# plt.plot(R, I, label = 'Infected endemic eq')
# plt.plot(R, I1, label = 'Infected disease-free eq')
# plt.plot(R, Q, label = 'Quarantined endemic eq')
# plt.plot(R, Q1, label = 'Quarantined disease-free eq')
plt.plot(R, C, label = 'Carrier endemic eq')
plt.plot(R, C1, label = 'Carrier disease-free eq')
plt.legend()
plt.grid()
plt.show()
# R = (beta*A*delta1)/((delta+q)*(delta+delta1)*d1_star) + (alpha*A*delta1*delta2*(1-p1))/((delta+q)*(delta+delta1)*d1_star*d3_star)
# print(R)

# E = d1*I/delta1
# Q = p1*delta2*I/d2
# C = (1-p1)*delta2*I/d3
