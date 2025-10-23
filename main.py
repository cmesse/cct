import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import frenet


# the basecurve
C = frenet.BasecurveCCT( 60, 60, 0.25, 68, 4 )

n = len(C.t)
x = np.zeros( n )
y = np.zeros( n )
z = np.zeros( n )

for k in range( n ) :
    r = C.r( C.t[k])
    x[k] = r[0]
    y[k] = r[1]
    z[k] = r[2]

fig = plt.figure()


ax = fig.add_subplot(111, projection='3d')

ax.plot(x,y,z,'-b')
ax.set_aspect('equal')

A = frenet.CrossSection(4,4,1 )
G = frenet.Geometry(C, A, air_radius=20.0, tape_res=1.0, air_res = 20.0)
G.save("/tmp/test.geo")

plt.show()
