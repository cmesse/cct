import numpy as np
from frenet.Curve import Curve
from frenet.Basecurve import Basecurve
from frenet.CrossSection import CrossSection
from frenet.Point import Point
from frenet.Surface import *

class Tape:

    def __init__(self, basecurve: Basecurve, cross_section: CrossSection, index: int ):
        self.index = index
        self.id = index + 1
        self.basecurve = basecurve
        self.cross_section = cross_section

        self.points_left  = []
        self.points_right = []

        # todo: to be set into parameter object
        self.resolution = 5
        self.make_ends = False
        self.delta_z = 200

        self._make_points()

        self.curves = []
        self.curveloops = []
        self.surfaces = []

        self._make_curves()

        self._make_surface()





        self.z0 = np.nan
        self.z1 = np.nan

    def _make_points(self):
        n = len(self.basecurve.t)

        for k in range(n):
            t = self.basecurve.t[k]
            m = self.basecurve.r(t)
            T = self.basecurve.transform(t)

            # left point
            p0 = np.array([
                self.cross_section.leftpoints[self.index][0],  # width direction (n)
                self.cross_section.leftpoints[self.index][1],  # thickness direction (b)
                0.0
            ])
            p = m.flatten() + T @ p0

            self.points_left.append(Point(p[0], p[1], p[2], self.resolution))

            # right point
            q0 = np.array([
                self.cross_section.rightpoints[self.index][0],
                self.cross_section.rightpoints[self.index][1],
                0.0
            ])
            q = m.flatten() + T @ q0

            self.points_right.append(Point(q[0], q[1], q[2], self.resolution))

        if self.make_ends:
            # get coordinates of first point


            # get derivative
            v = self.basecurve.v( self.basecurve.t[0])

            # get help parameter
            xi0 = -self.delta_z / v[2]

            # compute length of first segment
            s = self.basecurve.segment_length( self.basecurve.t[0], self.basecurve.t[1])

            # approximate number of steps
            n = round( self.delta_z / s )

            # values to add
            xi = np.linspace( 0, xi0, n )

            # crate left points
            x0 = self.points_left[0].x
            y0 = self.points_left[0].y
            z0 = self.points_left[0].z
            for k in range(1,n) :
                P = Point( x0 + xi[k] * v[0], y0+ xi[k] * v[1], z0 + xi[k]*v[2], self.resolution)
                self.points_left.insert( 0, P )

            # crate right points
            x0 = self.points_right[0].x
            y0 = self.points_right[0].y
            z0 = self.points_right[0].z
            for k in range(1, n):
                P = Point(x0 + xi[k] * v[0], y0 + xi[k] * v[1], z0 + xi[k] * v[2], self.resolution)
                self.points_right.insert(0, P)

            # end extension
            v = self.basecurve.v(self.basecurve.t[-1])

            # get help parameter
            xi0 = self.delta_z / v[2]

            # compute length of first segment
            s = self.basecurve.segment_length( self.basecurve.t[-2], self.basecurve.t[-1])

            # approximate number of steps
            n = round( self.delta_z / s )

            # values to add
            xi = np.linspace( 0, xi0, n )

            x0 = self.points_left[-1].x
            y0 = self.points_left[-1].y
            z0 = self.points_left[-1].z
            for k in range(1, n):
                P = Point(x0 + xi[k] * v[0], y0 + xi[k] * v[1], z0 + xi[k] * v[2], self.resolution)
                self.points_left.append(P)

            x0 = self.points_right[-1].x
            y0 = self.points_right[-1].y
            z0 = self.points_right[-1].z
            for k in range(1, n):
                P = Point(x0 + xi[k] * v[0], y0 + xi[k] * v[1], z0 + xi[k] * v[2], self.resolution)
                self.points_right.append(P)
    def _make_curves(self):

        F = Curve("Line")
        F.points.append(self.points_left[0])
        F.points.append(self.points_right[0])
        self.curves.append(F)

        R = Curve("Spline")
        R.points = self.points_right
        self.curves.append(R)

        B = Curve("Line")
        B.points.append(self.points_right[-1])
        B.points.append(self.points_left[-1])
        self.curves.append(B)

        L = Curve("Spline")
        for p in self.points_left :
            L.points.insert(0, p)

        self.curves.append(L)

    def _make_surface(self):
        L = CurveLoop()
        L.curves = self.curves
        self.curveloops.append( L )
        S = Surface()
        S.loops.append( L )
        self.surfaces.append(S)