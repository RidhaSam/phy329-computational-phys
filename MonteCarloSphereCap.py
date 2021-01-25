import numpy as np
import random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 100000
r = 2
d = 1.5
pointsIn = 0
trueVol = (np.pi*d**2)/3*(3*r-d)
cubeVol = 8*r**3

for i in range(n):
    x = rd.uniform(-r,r)
    y = rd.uniform(-r,r)
    z = rd.uniform(-r,r)
    dist = x**2 + y**2 + z**2
    if (dist < r**2 and (r-d) < z < r):
        pointsIn += 1
    else:
        pass

print(trueVol)
print(pointsIn/n*cubeVol)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(xPoints2, yPoints2, zPoints2, c='b')
plt.show()

