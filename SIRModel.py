import numpy as np 
from matplotlib import pyplot as plt 

a = 0.000285
r = 0.15
S_0 = 5000
I_0 = 1
R_0 = 0
dt = 0.01
endPoint = 40
n = (int) (40/dt)

t = np.arange(0,40,0.01)
S = np.zeros(n)
I = np.zeros(n)
R = np.zeros(n)

def sDerivative(s, i):
    dsdt = -a*s*i
    return dsdt

def iDerivative(s, i):
    dsdt = a*s*i-r*i
    return dsdt

def rDerivative(s, i):
    dsdt = r*i
    return dsdt

for i in range(n):
    S[i] = S_0
    I[i] = I_0
    R[i] = R_0

    S_1_e = S_0 + sDerivative(S_0, I_0)*dt
    I_1_e = I_0 + iDerivative(S_0, I_0)*dt
    R_1_e = R_0 + rDerivative(S_0, I_0)*dt

    S_1 = S_0 + 0.5*(sDerivative(S_0,I_0) + sDerivative(S_1_e,I_1_e))*dt
    I_1 = I_0 + 0.5*(iDerivative(S_0,I_0) + iDerivative(S_1_e,I_1_e))*dt
    R_1 = R_0 + 0.5*(rDerivative(S_0,I_0) + rDerivative(S_1_e,I_1_e))*dt

    S_0 = S_1
    R_0 = R_1
    I_0 = I_1

fig = plt.figure()
plt.plot(t,S,'r')
plt.plot(t,I,'b')
plt.plot(t,R,'g')
plt.show()




