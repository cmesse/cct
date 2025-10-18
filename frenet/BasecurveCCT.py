import math

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

        self.tmin = 0
        self.tmax = nturns * np.pi * 2 + np.pi
        n = nturns * self.num_points_per_turn

        self.t = np.linspace( self.tmin, self.tmax, n )

        #self.t, self.s = self.make_equidistant(self.tmin, self.tmax, nturns * self.num_points_per_turn )

        self._ta = self.tmin + np.pi
        self._tb = self.tmax - np.pi

        self.cx0 = np.zeros(8)
        self.cz0 = np.zeros(8)
        self.cx1 = np.zeros(8)
        self.cz1 = np.zeros(8)

        self._is_initialized = False

        self._init_polys()

    def r( self, t: float):
        x = np.zeros(3)

        if t < self._ta and self._is_initialized :
            x[0] =  np.polyval( self.cx0, t )
            f = x[0]/self.R1
            x[1] = self.R2 * math.sqrt(1-f*f)
            x[2] = np.polyval( self.cz0, t )
        elif t > self._tb and self._is_initialized:
            x[0] = np.polyval(self.cx1, t)
            f = x[0] / self.R1
            x[1] = self.R2 * math.sqrt(1 - f * f)
            x[2] = np.polyval(self.cz1, t)
        else:
            # Russenschuck (3.34)
            x[0] = self.R1 * np.cos(t)
            x[1] = self.R2 * np.sin(t)
            x[2] = self.R2 * ( np.sin(t) * self.tan_alpha + self.q * t )

        return x

    def v( self, t: float):
        dxdt = np.zeros(3)
        dt = 1e-6

        if t < self._ta and self._is_initialized :

            dxdt[0] =  np.dot(self.cx0, self._deriv1( t ) )

            x0 = np.polyval( self.cx0, t-0.5*dt )
            y0 = self.R2 * math.sqrt(1-x0*x0/(self.R1*self.R1))

            x1 = np.polyval( self.cx0, t+0.5*dt )
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            dxdt[1] = (y1-y0)/dt
            dxdt[2] = np.dot(self.cz0, self._deriv1( t ) )

        elif t > self._tb and self._is_initialized :
            dxdt[0] =  np.dot(self.cx1, self._deriv1( t ) )
            x0 = np.polyval( self.cx1, t-0.5*dt )
            y0 = self.R2 * math.sqrt(1-x0*x0/(self.R1*self.R1))
            x1 = np.polyval( self.cx1, t+0.5*dt )
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            dxdt[1] = (y1-y0)/dt
            dxdt[2] = np.dot(self.cz1, self._deriv1( t ) )
        else:
            dxdt[0] = -self.R1 * np.sin(t)
            dxdt[1] = self.R2 * np.cos(t)
            dxdt[2] = self.R2 * ( np.cos(t) * self.tan_alpha + self.q )
        return dxdt

    def a( self, t: float ):

        d2xdt2 = np.zeros(3)
        dt = 1e-6

        if t < self._ta and self._is_initialized :
            d2xdt2[0] = np.dot(self.cx0, self._deriv2(t))



            x0 = np.polyval(self.cx0, t - dt)
            x1 = np.polyval(self.cx0, t)
            x2 = np.polyval(self.cx0, t + dt)

            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))

            dy0 = (y1-y0)/dt
            dy1 = (y2-y1)/dt
            d2xdt2[1] = (dy1-dy0)/dt
            d2xdt2[2] = np.dot(self.cz0, self._deriv2(t))
        elif t > self._tb and self._is_initialized:
            d2xdt2[0] = np.dot(self.cx1, self._deriv2(t))

            x0 = np.polyval(self.cx1, t - dt)
            x1 = np.polyval(self.cx1, t)
            x2 = np.polyval(self.cx1, t + dt)

            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt
            dy1 = (y2 - y1) / dt
            d2xdt2[1] = -(dy1 - dy0) / dt
            d2xdt2[2] = np.dot(self.cz1, self._deriv2(t))
        else:
            d2xdt2[0] = -self.R1 * np.cos(t)
            d2xdt2[1] = -self.R2 * np.sin(t)
            d2xdt2[2] =  self.R2 * ( -np.sin(t) * self.tan_alpha )

        return d2xdt2

    def b(self, t: float ):
        d3xdt3 = np.zeros(3)
        dt = 1e-6
        if t < self._ta and self._is_initialized :
            x0 = np.polyval(self.cx0, t - 1.5*dt)
            x1 = np.polyval(self.cx0, t - 0.5*dt)
            x2 = np.polyval(self.cx0, t + 0.5*dt)
            x3 = np.polyval(self.cx0, t + 1.5*dt)
            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))
            y3 = self.R2 * math.sqrt(1 - x3 * x3 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt # @ - dt
            dy1 = (y2 - y1) / dt # @ 0
            dy2 = (y3 - y2) / dt # @ + dt

            ddy0 = (dy1-dy0)  / dt # @ -0.5 dt
            ddy1 = ( dy2-dy1) / dt # @ 0.5 dt

            d3xdt3[0] = np.dot(self.cx0, self._deriv3(t))
            d3xdt3[1] = (ddy1-ddy0)/dt
            d3xdt3[2] = np.dot(self.cz0, self._deriv3(t))
        elif t > self._tb and self._is_initialized :
            x0 = np.polyval(self.cx1, t - 1.5*dt)
            x1 = np.polyval(self.cx1, t - 0.5*dt)
            x2 = np.polyval(self.cx1, t + 0.5*dt)
            x3 = np.polyval(self.cx1, t + 1.5*dt)
            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))
            y3 = self.R2 * math.sqrt(1 - x3 * x3 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt # @ - dt
            dy1 = (y2 - y1) / dt # @ 0
            dy2 = (y3 - y2) / dt # @ + dt

            ddy0 = (dy1-dy0)  / dt # @ -0.5 dt
            ddy1 = ( dy2-dy1) / dt # @ 0.5 dt

            d3xdt3[0] = np.dot(self.cx1, self._deriv3(t))
            d3xdt3[1] = (ddy1-ddy0)/dt
            d3xdt3[2] = np.dot(self.cz1, self._deriv3(t))
        else:
            d3xdt3[0] = self.R1 * np.sin(t)
            d3xdt3[1] = -self.R2 * np.cos(t)
            d3xdt3[2] = self.R2 * (-np.cos(t) * self.tan_alpha)
        return d3xdt3

    def _init_polys(self):

        f0 = np.zeros(3)
        f0[0] = 0
        f0[1] = self.R2
        f0[2] = 0
        df0   = 0.5*self.v(self.tmin)
        ddf0  = np.zeros(3)
        dddf0 = np.zeros(3)

        f1 = self.r(self._ta)
        df1 = self.v(self._ta)
        ddf1 = self.a(self._ta)
        dddf1 = self.b(self._ta)

        df = f1-f0

        self.cx0 = self._compute_values(0, self.tmin, f0, df0, ddf0, dddf0, self._ta, f1, df1, ddf1, dddf1 )
        self.cz0 = self._compute_values(2, self.tmin, f0, df0, ddf0, dddf0, self._ta, f1, df1, ddf1, dddf1)

        dz = self.r(self.tmin+0.5*np.pi)[2]-self.r(self.tmin)[2]
        print("r0", self.r(self.tmin+0.5*np.pi))
        print("r1", self.r(self.tmax-0.5*np.pi))

        r = self.r(self.tmax - 0.5*np.pi )
        f0[0] = 0
        f0[1] = -self.R2
        f0[2] = r[2] + dz

        df0 = -0.5 * self.v(self.tmax)

        f1 = self.r(self._tb)
        df1 = self.v(self._tb)
        ddf1 = self.a(self._tb)
        dddf1 = self.b(self._tb)

        self.cx1 = self._compute_values(0, self.tmax, f0, df0, ddf0, dddf0, self._tb, f1, df1, ddf1, dddf1)
        self.cz1 = self._compute_values(2, self.tmax, f0, df0, ddf0, dddf0, self._tb, f1, df1, ddf1, dddf1)

        self._is_initialized  = True

    def _compute_values(self, k: int, t0: float , f0: float, df0: float, ddf0: float, dddf0: float, t1: float, f1: float, df1: float, ddf1: float, dddf1: float ):

        V = np.zeros([8,8])
        V[0][0] = t0 * t0 * t0 * t0 * t0 * t0 * t0
        V[0][1] = t0 * t0 * t0 * t0 * t0 * t0
        V[0][2] = t0 * t0 * t0 * t0 * t0
        V[0][3] = t0 * t0 * t0 * t0
        V[0][4] = t0 * t0 * t0
        V[0][5] = t0 * t0
        V[0][6] = t0
        V[0][7] = 1
        V[1][0] = 7 * t0 * t0 * t0 * t0 * t0 * t0
        V[1][1] = 6 * t0 * t0 * t0 * t0 * t0
        V[1][2] = 5 * t0 * t0 * t0 * t0
        V[1][3] = 4 * t0 * t0 * t0
        V[1][4] = 3 * t0 * t0
        V[1][5] = 2 * t0
        V[1][6] = 1
        V[1][7] = 0
        V[2][0] = 42 * t0 * t0 * t0 * t0 * t0
        V[2][1] = 30 * t0 * t0 * t0 * t0
        V[2][2] = 20 * t0 * t0 * t0
        V[2][3] = 12 * t0 * t0
        V[2][4] = 6 * t0
        V[2][5] = 2
        V[2][6] = 0
        V[2][7] = 0
        V[3][0] = 210 * t0 * t0 * t0 * t0
        V[3][1] = 120 * t0 * t0 * t0
        V[3][2] = 60 * t0 * t0
        V[3][3] = 24 * t0
        V[3][4] = 6
        V[3][5] = 0
        V[3][6] = 0
        V[3][7] = 0
        V[4][0] = t1 * t1 * t1 * t1 * t1 * t1 * t1
        V[4][1] = t1 * t1 * t1 * t1 * t1 * t1
        V[4][2] = t1 * t1 * t1 * t1 * t1
        V[4][3] = t1 * t1 * t1 * t1
        V[4][4] = t1 * t1 * t1
        V[4][5] = t1 * t1
        V[4][6] = t1
        V[4][7] = 1
        V[5][0] = 7 * t1 * t1 * t1 * t1 * t1 * t1
        V[5][1] = 6 * t1 * t1 * t1 * t1 * t1
        V[5][2] = 5 * t1 * t1 * t1 * t1
        V[5][3] = 4 * t1 * t1 * t1
        V[5][4] = 3 * t1 * t1
        V[5][5] = 2 * t1
        V[5][6] = 1
        V[5][7] = 0
        V[6][0] = 42 * t1 * t1 * t1 * t1 * t1
        V[6][1] = 30 * t1 * t1 * t1 * t1
        V[6][2] = 20 * t1 * t1 * t1
        V[6][3] = 12 * t1 * t1
        V[6][4] = 6 * t1
        V[6][5] = 2
        V[6][6] = 0
        V[6][7] = 0
        V[7][0] = 210 * t1 * t1 * t1 * t1
        V[7][1] = 120 * t1 * t1 * t1
        V[7][2] = 60 * t1 * t1
        V[7][3] = 24 * t1
        V[7][4] = 6
        V[7][5] = 0
        V[7][6] = 0
        V[7][7] = 0

        f = np.zeros(8)
        f[0] = f0[k]
        f[1] = df0[k]
        f[2] = ddf0[k]
        f[3] = dddf0[k]
        f[4] = f1[k]
        f[5] = df1[k]
        f[6] = ddf1[k]
        f[7] = dddf1[k]

        c = np.linalg.solve( V, f )
        return c

    def _deriv1(self, t: float ):
        f = np.zeros(8)
        f[0] = 7 * t * t * t * t * t * t
        f[1] = 6 * t * t * t * t * t
        f[2] = 5 * t * t * t * t
        f[3] = 4 * t * t * t
        f[4] = 3 * t * t
        f[5] = 2 * t
        f[6] = 1
        return f

    def _deriv2(self, t: float ):
        f = np.zeros(8)
        f[0] = 42 * t * t * t * t * t
        f[1] = 30 * t * t * t * t
        f[2] = 20 * t * t * t
        f[3] = 12 * t * t
        f[4] = 6 * t
        f[5] = 2
        return f
    def _deriv3(self, t: float ):
        f = np.zeros(8)
        f[0] = 210 * t * t * t * t
        f[1] = 120 * t * t * t
        f[2] = 60 * t * t
        f[3] = 24 * t
        f[4] = 6
        return f