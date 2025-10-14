import numpy as np


class Basecurve :

    def __init__(self):
        self._nintpoints = 7
        self._intpoints, self._weights = np.polynomial.legendre.leggauss(self._nintpoints)
        self.t = None

    # the actual Basecurve
    def r( self, t: float ):
        raise NotImplementedError()

    # the velocity
    def v( self, t: float ):
        raise NotImplementedError()

    # the accelleration
    def a( self, t: float ):
        raise NotImplementedError()

    # the jerk
    def b( self, t: float ):
        raise NotImplementedError()

    def segment_length(self, ta: float, tb: float ):

        l = 0
        for k in range( self._nintpoints ) :
            t = 0.5*((1-self._intpoints[k])*ta + ( 1 + self._intpoints[k])*tb)
            v = self.v(t)
            l += self._weights[k] * np.linalg.norm(v)

        l *= (tb-ta)*0.5

        return l

    def make_equidistant(self, ta: float, tb: float, n: int ):

        t = np.linspace( ta, tb, n )

        # compute the full length of the Basecurve
        l = 0
        for k in range(1,n):
            l += self.segment_length(t[k-1],t[k])

        # expected distance
        dl = l / (n-1)

        for k in range(1,n):
            # initial guesses
            dt = abs(t[k] - t[k-1])
            t0 = t[k-1]
            t1 = t0 + 0.01*dt


            f1 = (self.segment_length(t0, t1) - dl) / dl
            f2 = f1
            t2 = t1
            while f1*f2 > 0 :
                t2 = t2 + dt
                f2 = (self.segment_length(t0, t2)-dl)/dl
            f3 = 1

            c = 0
            while abs(f3) > 1e-7:
                # intersection point
                if c < 20 :
                    t3 = t1 - f1*(t2-t1)/(f2-f1)
                else:
                    t3 = 0.5*(t1+t2)
                f3 = (self.segment_length(t0, t3)-dl)/dl

                if f1*f3 < 0 :
                    t2 = t3
                    f2 = f3
                else:
                    t1 = t3
                    f2 = f3
                c = c + 1

            t[k] = t3

        return t

    def transform(self, t: float ):
        v = self.v( t )
        a = self.a( t )

        # tangent
        T = v / np.linalg.norm(v)

        # binomial vector
        vxa = np.linalg.cross(v,a)
        B = vxa/np.linalg.norm(vxa)

        # normal vector
        N = np.linalg.cross(B,T)
        N /= np.linalg.norm(N)

        # transformation matrix
        R = np.zeros([3,3])
        R[:,0] = N
        R[:,1] = B
        R[:,2] = T

        return R

    def kappa_tau(self, t: float ):
        v = self.v(t)
        a = self.a(t)
        vxa = np.linalg.cross(v, a)

        nvxa = np.linalg.norm(vxa)

        kappa = nvxa/(v**3)
        tau = nvxa * self.b(t)/(nvxa**2)

        return kappa, tau
