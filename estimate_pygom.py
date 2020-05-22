from pygom import Transition, TransitionType, DeterministicOde, SquareLoss
from scipy.integrate import odeint
from scipy.optimize import minimize
import numpy as np
import pandas as pd

df = pd.read_excel("Italy Variable Data.xlsx",header = 0)
S = list(df['S'])
E = list(df['E'])
I = list(df['I'])
Q = list(df['Q'])
C = list(df['C'])
R = list(df['R'])
y = []
t = []
for i in range(1, len(S)):
    y.append([S[i], E[i], I[i], Q[i], R[i]])
    t.append(i)

states = ['S', 'E', 'I', 'Q', 'C', 'R']
params = ['A', 'q', 'mu1', 'mu2', 'mu3', 'd1', 'd2','d3', 'delta', 'delta1', 'delta2', 'p1', 'alpha', 'beta']
odeList = [ Transition(origin='S',equation='A - (delta*S) - (beta*S*I) - (alpha*C*S) - q*S',transition_type=TransitionType.ODE),
           Transition(origin='E',equation='(beta*I*S) - (delta*E) - (delta1*E) + (alpha*C*S)',transition_type=TransitionType.ODE),
           Transition(origin='I',equation='delta1*E - (delta*I) - (delta2*I) - (mu1*I) - (d1*I)',transition_type=TransitionType.ODE),
           Transition(origin='Q',equation='(p1*delta2*I) - (delta*Q) - (mu2*Q) - (d2*Q)',transition_type=TransitionType.ODE),
           Transition(origin='C',equation='(delta2*I*(1-p1)) - (delta*C) - (mu3*C) - (d3*C)',transition_type=TransitionType.ODE),
           Transition(origin='R',equation='mu1*I + mu2*Q + mu3*C - delta*R',transition_type=TransitionType.ODE)]
model = DeterministicOde (states, params, ode =odeList)
init_state = [60461803,2703,94, 127, 2948,1]
param_eval =  [('A', 0.007896),('q',0.01), ('mu1', 0.4554), ('mu2', 1.21382), ('mu3', 0.1325), ('d1',0.24203), ('d2', 0.55586),('d3', 0.07849), ('delta', 0.000213), ('delta1',0.196), ('delta2',0.996), ('alpha',0.000000196), ('beta', 0.000034196), ('p1',0.96)]
model.intial_values = (init_state, [0])
model.parameters = param_eval
# sol = odeint(model.ode, init_state, t[1:])
# print(sol)
# print(sol.size)
theta = [0.5, 0.9, 0.05, 0.05]
bounds = [(0,1),(0,1),(0,1),(0,1)]

objFH = SquareLoss(theta=theta, ode=model, x0=init_state, t0=[0], t=t, y=y, state_name=['S', 'E', 'I', 'Q', 'R'], target_param = ['delta2', 'p1', 'alpha', 'beta'])
res = minimize(fun=objFH.cost, jac=objFH.sensitivity, x0=theta, bounds = bounds, options={'disp': True})
