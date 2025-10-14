from frenet.AsciiFile import AsciiFile
from frenet.Basecurve import Basecurve
from frenet.CrossSection import CrossSection
from frenet.Tape import Tape
from frenet.TapeBlock import TapeBlock

class Geometry :

    def __init__(self, basecurve: Basecurve, cross_section: CrossSection ):
        self.basecurve = basecurve
        self.cross_section = cross_section
        self.tapes = []
        self.tape_blocks = []
        self.points = []
        self.curves = []
        self.loops = []
        self.surfaces = []

        self._create_tapes()
        self._create_tape_blocks()

    def _create_tapes(self):

        n = self.cross_section.numtapes
        for k in range(n):
            self.tapes.append(Tape(self.basecurve,self.cross_section,k))

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
        lid = 0
        sid = 0
        for t in self.tapes:
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.loops :
                lid+=1
                l.id = lid
                self.loops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)

        for t in self.tape_blocks :
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.loops :
                lid+=1
                l.id = lid
                self.loops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)
    def save(self,path: str):
        self._collect_points()
        self._collect_geometry()

        F = AsciiFile()
        for p in self.points :
            F.Buffer.append( p.write())

        for c in self.curves :
            F.Buffer.append( c.write())

        for l in self.loops :
            F.Buffer.append( l.write())

        for s in self.surfaces :
            F.Buffer.append(s.write())
        F.save(path)




