import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import frenet


# the basecurve
C = frenet.BasecurveCCT( 60, 60, 1, 68, 4 )
A = frenet.CrossSection(3,10,2)


G = frenet.Geometry(C,A)
G.save("/tmp/test.geo")