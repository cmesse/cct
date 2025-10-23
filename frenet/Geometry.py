from frenet.AsciiFile import AsciiFile
from frenet.Basecurve import Basecurve
from frenet.CrossSection import CrossSection
from frenet.Tape import Tape
from frenet.TapeBlock import TapeBlock
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop
from frenet.Volume import Volume, SurfaceLoop
import numpy as np

class Geometry :

    def __init__(self, basecurve: Basecurve, cross_section: CrossSection, air_radius: float = 50.0, tape_res: float = 5.0, air_res: float = 10.0 ):
        self.basecurve = basecurve
        self.cross_section = cross_section
        self.air_radius = air_radius
        self.air_res = air_res
        self.tape_res = tape_res
        self.tapes = []
        self.tape_blocks = []
        self.points = []
        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self.surfaceloops = []
        self.volumes = []

        # Air domain components
        self.air_points = []
        self.air_curves = []
        self.air_surfaces = []

        self._create_tapes()
        self._create_tape_blocks()

    def _create_tapes(self):

        n = self.cross_section.numtapes
        for k in range(n):
            self.tapes.append(Tape(self.basecurve,self.cross_section,k,self.tape_res))

    def _create_tape_blocks(self):
        if self.cross_section.numtapes > 1 :
            for k in range(1,self.cross_section.numtapes):
                self.tape_blocks.append(TapeBlock(self.tapes[k-1],self.tapes[k]))
        

    def _collect_points(self):
        self.points = []
        c = 0
        for t in self.tapes :
            for k in t.points_left :
                c += 1
                k.id = c
                self.points.append( k )
            for k in t.points_right :
                c += 1
                k.id = c
                self.points.append( k )

    def _collect_geometry(self):
        cid = 0
        clid = 0
        sid = 0
        slid = 0
        vid = 0
        for t in self.tapes:
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.curveloops :
                clid+=1
                l.id = clid
                self.curveloops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)

        for t in self.tape_blocks :
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.curveloops :
                clid+=1
                l.id = clid
                self.curveloops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)
            for l in t.surfaceloops :
                slid += 1
                l.id = slid
                self.surfaceloops.append( l )
            for v in t.volumes :
                vid += 1
                v.id = vid
                self.volumes.append( v )

    def _compute_bounding_box(self):
        """Compute bounding box of all coil points"""
        # Collect all points first
        self._collect_points()

        if len(self.points) == 0:
            return None

        # Initialize with first point
        x_min = x_max = self.points[0].x
        y_min = y_max = self.points[0].y
        z_min = z_max = self.points[0].z

        # Find min/max coordinates
        for p in self.points:
            x_min = min(x_min, p.x)
            x_max = max(x_max, p.x)
            y_min = min(y_min, p.y)
            y_max = max(y_max, p.y)
            z_min = min(z_min, p.z)
            z_max = max(z_max, p.z)

        return {
            'x_min': x_min, 'x_max': x_max,
            'y_min': y_min, 'y_max': y_max,
            'z_min': z_min, 'z_max': z_max
        }

    def _create_air_domain(self):
        """Create air domain box around the coil with coherent boundary surfaces"""
        bbox = self._compute_bounding_box()
        if bbox is None:
            return

        # Add margin in x and y directions, terminals are flush in z
        x_min = bbox['x_min'] - self.air_radius
        x_max = bbox['x_max'] + self.air_radius
        y_min = bbox['y_min'] - self.air_radius
        y_max = bbox['y_max'] + self.air_radius
        z_min = bbox['z_min']  # Flush with front terminal
        z_max = bbox['z_max']  # Flush with back terminal

        # Create 8 corner points of the air box
        # Front face (z_min)
        p1 = Point(x_min, y_min, z_min, self.air_res)
        p2 = Point(x_max, y_min, z_min, self.air_res)
        p3 = Point(x_max, y_max, z_min, self.air_res)
        p4 = Point(x_min, y_max, z_min, self.air_res)

        # Back face (z_max)
        p5 = Point(x_min, y_min, z_max, self.air_res)
        p6 = Point(x_max, y_min, z_max, self.air_res)
        p7 = Point(x_max, y_max, z_max, self.air_res)
        p8 = Point(x_min, y_max, z_max, self.air_res)

        p9 = Point(0.0, 0.0, z_min, self.air_res)
        p10 = Point(0.0, 0.0, z_max, self.air_res)

        self.air_points = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]

        # Create curves for the box edges
        # Front face outer boundary curves
        c_front_bottom = Curve("Circle")
        c_front_bottom.points = [p1, p9, p2]
        c_front_right = Curve("Circle")
        c_front_right.points = [p2, p9, p3]
        c_front_top = Curve("Circle")
        c_front_top.points = [p3, p9, p4]
        c_front_left = Curve("Circle")
        c_front_left.points = [p4, p9, p1]

        # Back face outer boundary curves
        c_back_bottom = Curve("Circle")
        c_back_bottom.points = [p5, p10, p6]
        c_back_right = Curve("Circle")
        c_back_right.points = [p6, p10, p7]
        c_back_top = Curve("Circle")
        c_back_top.points = [p7, p10, p8]
        c_back_left = Curve("Circle")
        c_back_left.points = [p8, p10, p5]

        # Connecting edges between front and back
        c_conn1 = Curve("Line")
        c_conn1.points = [p1, p5]
        c_conn2 = Curve("Line")
        c_conn2.points = [p2, p6]
        c_conn3 = Curve("Line")
        c_conn3.points = [p3, p7]
        c_conn4 = Curve("Line")
        c_conn4.points = [p4, p8]

        self.air_curves = [
            c_front_bottom, c_front_right, c_front_top, c_front_left,  # 0-3
            c_back_bottom, c_back_right, c_back_top, c_back_left,      # 4-7
            c_conn1, c_conn2, c_conn3, c_conn4                         # 8-11
        ]

        # Front face: outer boundary with holes for coil
        cl_front_outer = CurveLoop()
        cl_front_outer.curves = [c_front_bottom, c_front_right, c_front_top, c_front_left]
        cl_front_outer.signs = [1, 1, 1, 1]

        s_front = Surface()
        s_front.is_plane = True
        s_front.loops.append(cl_front_outer)
        for t in self.tape_blocks:
            s_front.loops.append(t.front)

        # Back face: outer boundary with holes for coil
        cl_back_outer = CurveLoop()
        cl_back_outer.curves = [c_back_bottom, c_back_right, c_back_top, c_back_left]
        cl_back_outer.signs = [1, 1, 1, 1]

        s_back = Surface()
        s_back.is_plane = True
        s_back.loops.append(cl_back_outer)
        for t in self.tape_blocks:
            s_back.loops.append(t.back)

        # Bottom face (y_min)
        cl_bottom = CurveLoop()
        cl_bottom.curves = [c_front_bottom, c_conn2, c_back_bottom, c_conn1]
        cl_bottom.signs = [1, 1, -1, -1]
        s_bottom = Surface()
        s_bottom.is_plane = False
        s_bottom.loops.append(cl_bottom)

        # Right face (x_max)
        cl_right = CurveLoop()
        cl_right.curves = [c_front_right, c_conn3, c_back_right, c_conn2]
        cl_right.signs = [1, 1, -1, -1]
        s_right = Surface()
        s_right.is_plane = False
        s_right.loops.append(cl_right)

        # Top face (y_max)
        cl_top = CurveLoop()
        cl_top.curves = [c_front_top, c_conn4, c_back_top, c_conn3]
        cl_top.signs = [1, 1, -1, -1]
        s_top = Surface()
        s_top.is_plane = False
        s_top.loops.append(cl_top)

        # Left face (x_min)
        cl_left = CurveLoop()
        cl_left.curves = [c_front_left, c_conn1, c_back_left, c_conn4]
        cl_left.signs = [1, 1, -1, -1]
        s_left = Surface()
        s_left.is_plane = False
        s_left.loops.append(cl_left)

        self.air_surfaces = [s_front, s_back, s_bottom, s_right, s_top, s_left]
        self.air_curveloops = [cl_front_outer, cl_back_outer, cl_bottom, cl_right, cl_top, cl_left]

    def _add_air_domain_to_geometry(self):
        """Add air domain components to the geometry lists"""
        # Assign IDs to air points (starting after coil points)
        pid = len(self.points)
        for p in self.air_points:
            pid += 1
            p.id = pid
            self.points.append(p)

        # Assign IDs to air curves (starting after coil curves)
        cid = len(self.curves)
        for c in self.air_curves:
            cid += 1
            c.id = cid
            self.curves.append(c)

        # Assign IDs to air curve loops (starting after coil curve loops)
        clid = len(self.curveloops)
        for cl in self.air_curveloops:
            clid += 1
            cl.id = clid
            self.curveloops.append(cl)

        # Assign IDs to air surfaces (starting after coil surfaces)
        sid = len(self.surfaces)
        for s in self.air_surfaces:
            sid += 1
            s.id = sid
            self.surfaces.append(s)

        # Create surface loop for the air box (outer boundary)
        air_box_loop = SurfaceLoop()
        air_box_loop.surfaces = self.air_surfaces
        air_box_loop.signs = [1, -1, 1, 1, 1, 1]  # Proper orientation for outer boundary

        # Assign ID to air box surface loop
        slid = len(self.surfaceloops)
        slid += 1
        air_box_loop.id = slid
        self.surfaceloops.append(air_box_loop)

        # Add coil side surfaces to the air boundary
        for t in self.tape_blocks:
            air_box_loop.surfaces.append(t.left)
            air_box_loop.signs.append(-1)
            air_box_loop.surfaces.append(t.right)
            air_box_loop.signs.append(-1)
        air_box_loop.surfaces.append(self.tape_blocks[0].bottom)
        air_box_loop.signs.append(-1)
        air_box_loop.surfaces.append(self.tape_blocks[-1].top)
        air_box_loop.signs.append(-1)

        # Create air volume with coil as interior hole
        # Air volume = Air box (outer) - Coil (inner hole)
        # We reference the coil surface loops (already in self.surfaceloops) as interior boundaries
        air_volume = Volume()
        air_volume.loops.append(air_box_loop)

        # Add all coil surface loops as interior boundaries (holes)
        # These are the surface loops from TapeBlock volumes
        #for sl in self.surfaceloops[:-1]:  # All except air_cylinder_loop (which is last)
        #    air_volume.loops.append(sl)

        # Assign ID to air volume
        vid = len(self.volumes)
        vid += 1
        air_volume.id = vid
        self.volumes.append(air_volume)

    def save(self,path: str):
        self._collect_points()
        self._collect_geometry()

        self._create_air_domain()
        self._add_air_domain_to_geometry()

        F = AsciiFile()
        for p in self.points :
            F.Buffer.append( p.write())

        for c in self.curves :
            F.Buffer.append( c.write())

        for l in self.curveloops :
            F.Buffer.append( l.write())

        for s in self.surfaces :
            F.Buffer.append(s.write())

        for l in self.surfaceloops :
            F.Buffer.append( l.write() )

        for v in self.volumes :
            F.Buffer.append( v.write() )
        F.save(path)




