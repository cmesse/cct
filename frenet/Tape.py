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
        self.make_ends = True
        self.delta_z = 200

        self._make_points()

        self.curves = []
        self.loops = []
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
            t = self.basecurve.t[0]
            m = self.basecurve.r(t)
            self.z0 = m[2, 0] - self.delta_z
            v0 = self.basecurve.v(t)

            P0 = self.points_left[0]

            xi = (self.z0 - P0.z) / v0[2]
            P = Point(P0.x + xi * v0[0], P0.y + xi * v0[1], self.z0, self.resolution)

            Q0 = self.points_right[0]

            xi = (self.z0 - Q0.z) / v0[2]
            Q = Point(Q0.x + xi * v0[0], Q0.y + xi * v0[1], self.z0, self.resolution)

            self.points_left.insert(0, P)
            self.points_right.insert(0, Q)

            t = self.basecurve.t[-1]
            m = self.basecurve.r(t)
            self.z1 = m[2, 0] + self.delta_z
            v1 = self.basecurve.v(t)

            P1 = self.points_left[-1]

            xi = (self.z1 - P1.z) / v1[2]
            P = Point(P1.x + xi * v1[0], P1.y + xi * v1[1], self.z1, self.resolution)
            self.points_left.append(P)

            Q1 = self.points_right[-1]

            xi = (self.z1 - Q1.z) / v1[2]
            Q = Point(Q1.x + xi * v1[0], Q1.y + xi * v1[1], self.z1, self.resolution)
            self.points_right.append(Q)


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
        L = Loop()
        L.curves = self.curves
        self.loops.append( L )
        S = Surface()
        S.loops.append( L )
        self.surfaces.append(S)