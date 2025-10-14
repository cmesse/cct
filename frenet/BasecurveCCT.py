import numpy as np
from frenet.Basecurve import Basecurve


class BasecurveCCT( Basecurve ) :

    def __init__(self, R1: float, R2: float, pitch: float, angle: float, nturns: int ):

        Basecurve.__init__( self )

        # first radius
        self.R1 = R1

        # second radius
        self.R2 = R2

        # pitch divided by 1/(2*pi)
        self.q = pitch / (2.0 * np.pi )

        # tan of tilt angle
        self.tan_alpha = np.tan( angle * np.pi/180 )

        self.t, self.s = self.make_equidistant(0, nturns * np.pi * 2, nturns * self.num_points_per_turn )



    def r( self, t: float):
        x = np.zeros([3,1])

        # Russenschuck (3.34)
        x[0] = self.R1 * np.cos(t)
        x[1] = self.R2 * np.sin(t)
        x[2] = self.R2 * ( np.sin(t) * self.tan_alpha + self.q * t )
        return x

    def v( self, t):
        t = np.atleast_1d(t)
        dxdt = np.zeros((len(t), 3))

        dxdt[:, 0] = -self.R1 * np.sin(t)
        dxdt[:, 1] = self.R2 * np.cos(t)
        dxdt[:, 2] = self.R2 * ( np.cos(t) * self.tan_alpha + self.q )
        if dxdt.shape[0] == 1:
            return dxdt[0]
        return dxdt

    def a( self, t ):
        t = np.atleast_1d(t)
        d2xdt2 = np.zeros((len(t), 3))


        d2xdt2[:, 0] = -self.R1 * np.cos(t)
        d2xdt2[:, 1] = -self.R2 * np.sin(t)
        d2xdt2[:, 2] =  self.R2 * ( -np.sin(t) * self.tan_alpha )
        if d2xdt2.shape[0] == 1:
            return d2xdt2[0]
        return d2xdt2

    def b(self, t):
        t = np.atleast_1d(t)
        d3xdt3 = np.zeros((len(t), 3))

        d3xdt3[:, 0] = self.R1 * np.sin(t)
        d3xdt3[:, 1] = -self.R2 * np.cos(t)
        d3xdt3[:, 2] = self.R2 * (-np.cos(t) * self.tan_alpha)
        if d3xdt3.shape[0] == 1:
            return d3xdt3[0]
        return d3xdt3
