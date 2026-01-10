# Parametric Driver
import sympy, numpy, re, itertools
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import RigidMotion.PyTree as R
import Converter
import OCC

#============================================================
# name server for creating entities with unique names
#============================================================
__NameServer__ = {}
def getName(proposedName):
    global __NameServer__
    (name, __NameServer__) = getUniqueName(proposedName, __NameServer__)
    return name

# Retourne proposedName#count
def getUniqueName2(proposedName, server, sep='#'):
    namespl = proposedName.rsplit(sep, 1)
    if len(namespl) == 2:
        try: c = int(namespl[1]); name = namespl[0]
        except: name = proposedName
    else: name = proposedName

    if name not in server:
        name2 = name+sep+'1'
        server[name2] = 0
        server[name] = 1
        return (name2, server)
    else:
        c = server[name]; ret = 1
        while ret == 1:
            name2 = '%s'+sep+'%d'%(name,c)
            if name2 not in server: ret = 0
            else: ret = 1
            c += 1
        server[name2] = 0
        server[name] = c
        return (name2, server)

# Retourne proposedNamecount
def getUniqueName(proposedName, server):
    if proposedName in server:
        c = server[proposedName]+1
    else:
        c = 1
    server[proposedName] = c
    return (proposedName+str(c), server)

#============================================================
# Helpers
#============================================================
# from arrays, export zone
def exportEdges(edges):
    t = C.newPyTree(['EDGES'])
    b = Internal.getNodeFromName1(t, 'EDGES')
    for c, e in enumerate(edges):
        z = Internal.createZoneNode('edge%03d'%(c+1), e, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        b[2].append(z)
    return t

#============================================================
class Scalar( sympy.core.symbol.Symbol ):
    """Define a parametric scalar"""
    def __new__(cls, name=None, value=0., **assumptions):
        if name is None: name = getName("scalar")
        obj = sympy.core.symbol.Symbol.__new__(cls, name, **assumptions)
        return obj

    def __init__(self, name=None, value=0.):
        # scalar name is symbol name
        self.name = super().name

        # symbol sympy: self by derivation

        # instantiated value
        self.v = value
        # range
        self.range = None
        # register
        DRIVER.registerScalar(self)

    # return value
    def v(self):
        return self.v

    # print content
    def print(self, shift=0):
        #print(" "*shift, "name", self.name)
        print(" "*shift, "value", self.v)
        if self.range is not None:
            print(" "*shift, "range", self.range)

    # check if instantiated value are in range, return 1 if ok, 0 else
    def check(self):
        # check if value is in range
        if self.range is not None:
            if self.v < self.range[0]: return 0
            if self.v > self.range[1]: return 0
        return 1

    # is parameter free?
    def isFree(self):
        if self.range is not None: return True
        else: return False

Vec1 = Scalar # alias

#============================================================
class Vec2:
    """Define a parametric vector of two components"""
    def __init__(self, name=None, value=(0.0,0.0)):
        # name
        if name is not None: self.name = name
        else: self.name = getName("vec2")
        # parameters
        self.x = Scalar("x", value[0])
        self.y = Scalar("y", value[1])

    # return value
    def v(self):
        return (self.x.v, self.y.v)

    # print content
    def print(self, shift=0):
        print(" "*shift, "x")
        self.x.print(shift+4)
        print(" "*shift, "y")
        self.y.print(shift+4)

    # check instantiated values
    def check(self):
        ret = self.x.check()
        if ret == 0: return 0
        ret = self.y.check()
        if ret == 0: return 0
        return 1

#============================================================
class Point:
    """Define a parametric point"""
    def __init__(self, name=None, value=(0.0,0.0,0.0)):
        # name
        if name is not None: self.name = name
        else: self.name = getName("point")
        # parameters
        self.x = Scalar(self.name+".x", value[0])
        self.y = Scalar(self.name+".y", value[1])
        self.z = Scalar(self.name+".z", value[2])
        # register
        DRIVER.registerPoint(self)

    # return value
    def v(self):
        return (self.x.v, self.y.v, self.z.v)

    def print(self, shift=0):
        print(" "*shift, self.x.name)
        self.x.print(shift+4)
        print(" "*shift, self.y.name)
        self.y.print(shift+4)
        print(" "*shift, self.z.name)
        self.z.print(shift+4)

    def check(self):
        ret = self.x.check()
        if ret == 0: return 0
        ret = self.y.check()
        if ret == 0: return 0
        ret = self.z.check()
        if ret == 0: return 0
        return 1

Vec3 = Point # alias

#============================================================
class Grid:
    """Define a parametric grid (cartesian)"""
    def __init__(self, name=None, Xo=(0.0,0.0,0.0), Xf=(0.0,0.0,0.0), N=(2,2,2)):
        # name
        if name is not None: self.name = name
        else: self.name = getName("grid")
        # keep
        self.ni = N[0]
        self.nj = N[1]
        self.nk = N[2]
        self.xo = Xo[0]
        self.yo = Xo[1]
        self.zo = Xo[2]
        self.dx = Xf[0]-Xo[0]
        self.dy = Xf[1]-Xo[1]
        self.dz = Xf[2]-Xo[2]
        # parameters
        self.P = [[[Point("%s.P.%d.%d.%d"%(self.name, i,j,k), (self.xo+i*self.dx, self.yo+j*self.dy, self.zo+k*self.dz)) for k in range(self.nk)] for j in range(self.nj)] for i in range(self.ni)]
        # register
        DRIVER.registerGrid(self)

    # return value
    def v(self, T):
        (i,j,k) = T
        return (self.P[i][j][k].x.v, self.P[i][j],[k].y.v, self.P[i][j][k].z.v)

    def print(self, shift=0):
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.P[i][j][k].print(shift+4)

    def check(self):
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    ret = self.P[i][j][k].check()
                    if ret == 0: return 0
        return 1

#============================================================
# Entities
#============================================================
class Entity:
    """Define a parametric entity"""
    def __init__(self, name=None, listP=[], type=None, mesh=None):
        # name
        if name is not None: self.name = name
        else:
            if type is not None: self.name = getName(type)
            else: self.name = getName("entity")
        # entity type
        self.type = type
        # optional reference mesh
        self.mesh = mesh
        # parameters
        self.P = []; c = 1
        for P in listP:
            if isinstance(P, Scalar):
                self.P.append(P)
            elif isinstance(P, float):
                self.P.append(Scalar(name+'.P%d'%c, P))
            elif isinstance(P, int):
                self.P.append(Scalar(name+'.P%d'%c, P))
            elif isinstance(P, tuple) and len(P) == 3:
                self.P.append(Point(name+'.P%d'%c, P))
            elif isinstance(P, Point):
                self.P.append(P)
            elif isinstance(P, tuple) and len(P) == 2:
                self.P.append(Vec2(name+'.P%d'%c, P))
            elif isinstance(P, Vec2):
                self.P.append(P)
            elif isinstance(P, Grid):
                self.P.append(P)
            else:
                raise(ValueError, "Wrong argument.")
            c += 1

        # hook on cad
        self.update()

        # register
        DRIVER.registerEdge(self)

    def __del__(self):
        OCC.occ.freeHook(self.hook)

    def update(self):
        self.hook = OCC.createEmptyCAD()
        if self.type == "line":
            OCC._addLine(self.hook, self.P[0].v(), self.P[1].v())
        elif self.type == "polyline":
            s = len(self.P)
            n = numpy.zeros((3,s), dtype=numpy.float64)
            for c, p in enumerate(self.P): n[:,c] = p.v()
            OCC._addSpline(self.hook, n, 0, 1)
        elif self.type == "spline1": # by control points
            s = len(self.P)
            n = numpy.zeros((3,s), dtype=numpy.float64)
            for c, p in enumerate(self.P): n[:,c] = p.v()
            OCC._addSpline(self.hook, n, 0, 3)
        elif self.type == "spline2": # by approximated points
            s = len(self.P)
            n = numpy.zeros((3,s), dtype=numpy.float64)
            for c, p in enumerate(self.P): n[:,c] = p.v()
            OCC._addSpline(self.hook, n, 1, 3)
        elif self.type == "spline3": # by free form control points + mesh
            # self.P[0] is a Grid, mesh is an array

            # rebuild control points
            import Generator, Transform
            grid = self.P[0]
            ni, nj, nk = grid.ni, grid.nj, grid.nk
            cp = Generator.cart((0,0,0), (1,1,1), (ni,nj,nk))
            cp = Converter.addVars(cp, ['dx','dy','dz'])

            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        ind = i+j*ni+k*ni*nj
                        cp[1][0,ind] = grid.xo + i*grid.dx
                        cp[1][1,ind] = grid.yo + j*grid.dy
                        cp[1][2,ind] = grid.zo + k*grid.dz
                        cp[1][3,ind] = grid.P[i][j][k].x.v - grid.xo - i*grid.dx
                        cp[1][4,ind] = grid.P[i][j][k].y.v - grid.yo - j*grid.dy
                        cp[1][5,ind] = grid.P[i][j][k].z.v - grid.zo - k*grid.dz

            mesh = Transform.freeForm(self.mesh, cp)
            mesh = Transform.deform(mesh, ['dx','dy','dz'])

            # spline from approximated point
            OCC._addSpline(self.hook, mesh[1], 1, 3)

        elif self.type == "circle":
            OCC._addCircle(self.hook, self.P[0].v(), (0,0,1), self.P[1].v, 0)
        elif self.type == "arc":
            OCC._addArc(self.hook, self.P[0].v(), self.P[1].v(), self.P[2].v())
        elif self.type == "superellipse":
            OCC._addSuperEllipse(self.hook, self.P[0].v(), self.P[1].v(), self.P[2].v(), self.P[3].v, self.P[4].v)
        elif self.type == "naca4":
            import Geom
            naca = Geom.naca("%01d%01d%02d"%(int(self.P[0].v),int(self.P[1].v),int(self.P[2].v)), N=51, sharpte=True)
            OCC._addSpline(self.hook, naca[1], 1, 3)
        else:
            raise(ValueError, "Unknown entity type %s."%self.type)

    def print(self, shift=0):
        for c, P in enumerate(self.P):
            print(" "*shift, P.name)
            P.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        OCC.writeCAD(self.hook, fileName, format)

    # check parameters validity
    def check(self):
        for P in self.P:
            ret = P.check()
            if ret == 0: return 0
        return 1

#============================================================
# line from two parametric points
def Line(name=None, P1=(0.,0.,0.), P2=(0.,0.,0.)):
    return Entity(name, [P1, P2], type="line")

# polyline from list of parametric points
def Polyline(name=None, Points=[]):
    return Entity(name, Points, type="polyline")

# spline from parametric control points
def Spline1(name=None, CPs=[]):
    return Entity(name, CPs, type="spline1")

# spline from parametric approximated points
def Spline2(name=None, Ps=[]):
    return Entity(name, Ps, type="spline2")

# spline from parametric grid
def Spline3(name=None, PGrid=None, mesh=None):
    return Entity(name, [PGrid], mesh=mesh, type="spline3")

# circle from parametric center and radius
def Circle(name=None, C=(0.,0.,0.), R=1.):
    return Entity(name, [C, R], type="circle")

# superellipse from parametric center and R1, R2
def SuperEllipse(name=None, C=(0.,0.,0.), R1=1., R2=1., N=4, samples=36):
    return Entity(name, [C, R1, R2, N, samples], type="superellipse")

# arc from 3 parametric points
def Arc(name=None, P1=(0.,0.,0.), P2=(0.,0.,0.), P3=(0.,0.,0.)):
    return Entity(name, [P1, P2, P3], type="arc")

# naca from parametric 4 digits
def Naca(name=None, M=0., P=0., e=12.):
    return Entity(name, [M, P, e], type="naca4")

#============================================================
class Sketch():
    """Define a parametric sketch from a list of entities."""
    def __init__(self, name=None, listEntities=[], h=None):
        # name
        if name is not None: self.name = name
        else: self.name = getName("sketch")
        # dependant entities
        self.entities = listEntities
        # type
        self.type = "sketch"
        # hook
        self.hook = None
        # global parameters (always added)
        self.P = []
        P = Vec3('%s.position'%self.name, (0.,0.,0.))
        self.P.append(P)
        self.position = P
        P = Point('%s.rotCenter'%self.name, (0.,0.,0.))
        self.P.append(P)
        self.rotCenter = P
        P = Vec3('%s.rotAxis'%self.name, (0.,0.,1.))
        self.P.append(P)
        self.rotAxis = P
        P = Scalar('%s.rotAngle'%self.name, 0.)
        self.P.append(P)
        self.rotAngle = P
        # update hook
        self.update()
        # register
        DRIVER.registerSketch(self)
        # meshing: (hmin,hmax,hausd)
        if h is not None: self.h = h
        # meshing: list of distribs for each entity
        self.distribs = None
        # reference mesh for Dmesh
        self.refMesh = None # arrays
        self.RefMesh = None # zone list (pytree)

    def add(self, entity):
        """Add an entity to sketch"""
        self.entities.append(entity)
        self.update()

    # update the CAD hook from parameters
    def update(self):
        """Update CAD hook."""
        if self.hook is not None: OCC.occ.freeHook(self.hook)
        self.hook = OCC.createEmptyCAD('sketch.step')
        hooks = []
        for e in self.entities: hooks.append(e.hook)
        self.hook = OCC.mergeCAD(hooks)
        # global positionning
        OCC._rotate(self.hook, self.P[1].v(), self.P[2].v(), self.P[3].v)
        OCC._translate(self.hook, self.P[0].v())

    # check if parameters are valid
    def check(self):
        """Check parameters validity."""
        for e in self.entities:
            ret = e.check()
            if ret == 1: return 1
        return 0

    # print information
    def print(self, shift=0):
        """Print informations."""
        for e in self.entities:
            print(" "*shift, e.name)
            e.print(shift+4)
        for e in [self.position, self.rotCenter, self.rotAxis, self.rotAngle]:
            print(" "*shift, e.name)
            e.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        """Write CAD to file."""
        OCC.writeCAD(self.hook, fileName, format)

    # mesh sketch, return arrays
    def mesh(self):
        """Mesh edges."""
        if self.distribs is not None:
            edges = []
            for c, e in enumerate(self.distrib):
                e = OCC.occ.meshOneEdge(self.hook, c+1, -1, -1, -1, -1, e)
            edges.append(e)
        elif self.h is not None:
            edges = OCC.meshAllEdges(self.hook, self.h[0], self.h[1], self.h[2], -1)
        else:
            raise ValueError("mesh: no meshing settings in sketch.")
        return edges

    # mesh sketch, return zones
    def Mesh(self):
        """Mesh edges."""
        edges = self.mesh()
        out = []
        for c, e in enumerate(edges):
            z = Internal.createZoneNode('%s%03d'%(self.name, c+1), e, [],
                                        Internal.__GridCoordinates__,
                                        Internal.__FlowSolutionNodes__,
                                        Internal.__FlowSolutionCenters__)
            out.append(z)
        return out

    # mesh sketch (arrays)
    def meshAsReference(self):
        """Mesh and store mesh as reference."""
        self.refMesh = self.mesh()
        return self.refMesh

    # mesh sketch (zones)
    def MeshAsReference(self):
        """Mesh and store mesh as reference."""
        self.RefMesh = self.Mesh()
        return self.RefMesh

    # Compute a rmesh identically to a reference mesh that is topologically
    # equivalent (same names). copy distributions, return arrays, remesh on CAD.
    def rmesh(self, refEdges):
        import Geom
        s = Geom.getCurvilinearAbscissa(refEdges)
        for c, e in enumerate(refEdges):
            e = Converter.rmVars(e, ['u'])
            refEdges[c] = Converter.addVars([e, s[c]])
        out = []
        for c, e in enumerate(refEdges):
            e = OCC.occ.meshOneEdge(self.hook, c+1, -1, -1, -1, -1, e)
            out.append(e)
        return out

    # Compute a rmesh identically to reference mesh that is topologically
    # equivalent (same names). copy distributions, return arrays
    def Rmesh(self, RefEdges):
        arrays = C.getAllFields(RefEdges, 'nodes', api=3)
        edges = self.rmesh(arrays)
        out = []
        for c, e in enumerate(edges):
            z = Internal.createZoneNode('%s%03d'%(self.name, c+1), e, [],
                                        Internal.__GridCoordinates__,
                                        Internal.__FlowSolutionNodes__,
                                        Internal.__FlowSolutionCenters__)
            out.append(z)
        return out

    def dmesh(self):
        """Remesh using reference mesh distributions."""
        return self.rmesh(self.refMesh)

    def Dmesh(self):
        """Remesh using reference mesh distributions."""
        return self.Rmesh(self.RefMesh)

#============================================================
class Surface():
    """Define a parametric surface."""
    def __init__(self, name=None, listSketches=[], listSketches2=[], listSurfaces=[], listSurfaces2=[], data={}, h=None, type="loft"):
        # name
        if name is not None: self.name = name
        else: self.name = getName("surface")
        # type
        self.type = type
        # dependant sketches
        self.sketches = listSketches
        # optional sketches
        self.sketches2 = listSketches2
        # dependant surfaces
        self.surfaces = listSurfaces
        # optional surfaces
        self.surfaces2 = listSurfaces2
        # optional data
        self.data = data
        # hook
        self.hook = None
        # global parameters (always added)
        self.P = []
        P = Vec3('%s.position'%self.name, (0,0,0))
        self.P.append(P)
        self.position = P
        P = Point('%s.rotCenter'%self.name, (0,0,0))
        self.P.append(P)
        self.rotCenter = P
        P = Vec3('%s.rotAxis'%self.name, (0,0,1))
        self.P.append(P)
        self.rotAxis = P
        P = Scalar('%s.rotAngle'%self.name, 0.)
        self.P.append(P)
        self.rotAngle = P
        # update hook
        self.update()
        # register
        DRIVER.registerSurface(self)
        # meshing: (hmin,hmax,hausd). supersedes sketch settings.
        if h is not None: self.h = h
        # reference mesh for Dmesh
        self.RefMesh = None # mesh (pytree)
        self.RefMeshUV = None # mesh in uv space (pytree)
        self.RefMeshUV2 = None # mesh in uv for deformation
        self.refEdges = None # edges of ref mesh (arrays)
        self.inds = None # indices of borders in refMesh BC
        self.inds2 = None # indices of borders in refMesh BC with kplane
        self.DefTree = None # def tree (pytree)

    def add(self, sketch):
        """Add a sketch to the surface definition."""
        self.sketches.append(sketch)

    # update the CAD from parameters
    def update(self):
        """Update CAD hook from parameters."""
        if self.hook is not None: OCC.occ.freeHook(self.hook)
        self.hook = OCC.createEmptyCAD('surface.step')
        if self.type == "loft":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            n1 = len(self.sketches)
            edgeList = [i for i in range(1, n1+1)]
            # optional guides
            for e in self.sketches2: hooks.append(e.hook)
            n2 = n1 + len(self.sketches2)
            guideList = [i for i in range(n1+1, n2+1)]
            self.hook = OCC.mergeCAD(hooks)
            OCC._loft(self.hook, edgeList, guideList)
            if 'close' in self.data and self.data['close']:
                h0 = hooks[0]; h1 = hooks[-1]
                OCC._fillHole(h0, [1], [], 0)
                OCC._fillHole(h1, [1], [], 0)
                self.hook = OCC.occ.mergeCAD([h0,self.hook,h1])
        if self.type == "loftSet":
            hooks = []
            for e in self.sketches:
                hooks.append(e.hook)
            n = len(self.sketches)
            out = []
            for i in range(1,n):
                h0 = hooks[i-1]
                h1 = hooks[i]
                hook = OCC.mergeCAD(hooks)
                OCC._loft(hook, [1,2])
                out.append(hook)
            if len(out) > 1:
                self.hook = OCC.mergeCAD(out)
            else: self.hook = out[0]
            if 'close' in self.data and self.data['close']:
                h0 = hooks[0]; h1 = hooks[-1]
                OCC._fillHole(h0, [1], [], 0)
                OCC._fillHole(h1, [1], [], 0)
                self.hook = OCC.occ.mergeCAD([h0,self.hook,h1])
        elif self.type == "revolve":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.mergeCAD(hooks)
            nedges = OCC.getNbEdges(self.hook)
            edgeList = [i for i in range(1, nedges+1)]
            OCC._revolve(self.hook, edgeList, self.data['center'], self.data['axis'], self.data['angle'])
        elif self.type == "merge":
            hooks = []
            for e in self.surfaces: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
        elif self.type == "fill":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            nedges = OCC.getNbEdges(self.hook)
            edgeList = [i for i in range(1, nedges+1)]
            OCC._fillHole(self.hook, edgeList, [], self.data['continuity'])
        elif self.type == "mergeEdges": # for debug
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
        elif self.type == "union":
            hooks = []; n1 = 0; n2 = 0
            for e in self.surfaces:
                n1 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            for e in self.surfaces2:
                n2 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            rev1 = self.data.get('rev1',1)
            rev2 = self.data.get('rev2',1)
            OCC._boolean(self.hook, [i for i in range(1,n1+1)], [i for i in range(n1+1,n1+n2+1)], 0, rev1, rev2)
        elif self.type == "inter":
            hooks = []; n1 = 0; n2 = 0
            for e in self.surfaces:
                n1 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            for e in self.surfaces2:
                n2 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            rev1 = self.data.get('rev1',1)
            rev2 = self.data.get('rev2',1)
            OCC._boolean(self.hook, [i for i in range(1,n1+1)], [i for i in range(n1+1,n1+n2+1)], 2, rev1, rev2)
        elif self.type == "sub":
            hooks = []; n1 = 0; n2 = 0
            for e in self.surfaces:
                n1 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            for e in self.surfaces2:
                n2 += OCC.getNbFaces(e.hook)
                hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            rev1 = self.data.get('rev1',1)
            rev2 = self.data.get('rev2',1)
            OCC._boolean(self.hook, [i for i in range(1,n1+1)], [i for i in range(n1+1,n1+n2+1)], 1, rev1, rev2)

        # global positionning
        OCC._rotate(self.hook, self.P[1].v(), self.P[2].v(), self.P[3].v)
        OCC._translate(self.hook, self.P[0].v())

    # print information
    def print(self, shift=0):
        """Print surface information."""
        for e in self.sketches:
            print(" "*shift, e.name)
            e.print(shift+4)
        for e in [self.position, self.rotCenter, self.rotAxis, self.rotAngle]:
            print(" "*shift, e.name)
            e.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        """Export to CAD file."""
        OCC.writeCAD(self.hook, fileName, format)

    # mesh surface, return arrays
    def mesh(self, close=True):
        """Mesh surface."""
        if self.h is None: raise ValueError("mesh: h settings are undefined.")
        (hmin,hmax,hausd) = self.h
        edges = OCC.meshAllEdges(self.hook, hmin, hmax, hausd, -1)
        nbFaces = OCC.getNbFaces(self.hook)
        faceList = range(1, nbFaces+1)
        if hausd < 0: hList = [(hmax,hmax,hausd)]*len(faceList)
        else: hList = [(hmin,hmax,hausd)]*len(faceList)
        faces = OCC.meshAllFacesTri(self.hook, edges, True, faceList, hList, close)
        return faces

    # mesh surface, return zones
    def Mesh(self, close=True):
        """Mesh surface."""
        faces = self.mesh(close)
        out = []
        for c, e in enumerate(faces):
            z = Internal.createZoneNode('%s%03d'%(self.name, c+1), e, [],
                                        Internal.__GridCoordinates__,
                                        Internal.__FlowSolutionNodes__,
                                        Internal.__FlowSolutionCenters__)
            out.append(z)
        return out

    def MeshAsReference(self):
        """Mesh surface and store the mesh for future Dmesh."""
        import Transform.PyTree as T
        import Transform
        import Post.PyTree as P

        # keep ref edges in x,y
        if self.h is None: raise ValueError("mesh: h settings are undefined.")
        (hmin,hmax,hausd) = self.h
        dedges = OCC.meshAllEdges(self.hook, hmin, hmax, hausd, -1)
        self.refEdges = dedges

        # build refMesh in x,y (not closed)
        self.RefMesh = self.Mesh(close=False)

        # save refMesh in UV
        self.RefMeshUV = Internal.copyRef(self.RefMesh)

        self.inds = []
        self.inds2 = []
        self.inds = {}
        self.inds2 = {}
        self.wires = {}
        self.dx1 = {}
        self.dy1 = {}
        self.dz1 = {}
        zones = Internal.getZones(self.RefMeshUV)
        for i, z in enumerate(zones): # pour chaque face

            # init inds
            self.inds[i+1] = []
            self.inds2[i+1] = []
            self.wires[i+1] = []

            # mesh each edges
            wires = OCC.occ.meshEdgesOfFace(self.hook, i+1, dedges)
            # switch edge to u,v
            for loops in wires:
                for w in loops:
                    w[1][0,:] = w[1][3,:]
                    w[1][1,:] = w[1][4,:]
                    w[1][2,:] = 0.
            self.wires[i+1] = wires

            # switch z to uv
            pu = Internal.getNodeFromName2(z, 'u')
            pv = Internal.getNodeFromName2(z, 'v')

            px = Internal.copyNode(pu)
            py = Internal.copyNode(pv)
            pz = Internal.copyNode(pv)
            pz[1][:] = 0.

            Internal.getNodeFromName2(z, 'CoordinateX')[1] = px[1]
            Internal.getNodeFromName2(z, 'CoordinateY')[1] = py[1]
            Internal.getNodeFromName2(z, 'CoordinateZ')[1] = pz[1]

            C._convertArray2NGon(z, recoverBC=False)
            p = P.exteriorFaces(z)
            C._addBC2Zone(z, 'wall', 'BCWall', subzone=p)
            Internal.getNodeFromType(z, 'BC_t')[0] = 'wall'
            T._addkplane(z)
            T._contract(z, (0,0,0), (1,0,0), (0,1,0), 0.1)

            # identify edges in bc
            bc = C.extractBCOfType(z, "BCWall", reorder=False)[0]
            xc = Internal.getNodeFromName2(bc, "CoordinateX")
            self.dx1[i+1] = numpy.empty((xc[1].size), dtype=numpy.float64)
            self.dy1[i+1] = numpy.empty((xc[1].size), dtype=numpy.float64)
            self.dz1[i+1] = numpy.empty((xc[1].size), dtype=numpy.float64)
            a = C.getAllFields(bc, 'nodes', api=3)[0]
            hook = Converter.createHook(a, function='nodes')
            for loops in wires:
                out = []
                for w in loops:
                    nodes = Converter.identifyNodes(hook, w)
                    out.append(nodes)
                self.inds[i+1].append(out)
            for loops in wires:
                out = []
                for w in loops:
                    w = Transform.translate(w, (0,0,0.1))
                    nodes = Converter.identifyNodes(hook, w)
                    out.append(nodes)
                self.inds2[i+1].append(out)
            Converter.freeHook(hook)
        self.RefMeshUV2 = zones

        # complete refmesh with gridinit containing initial UV
        zones = Internal.getZones(self.RefMesh)
        for z in zones:
            gridInit = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
            if not gridInit:
                gridInit = Internal.createNode('GridCoordinates#Init', 'GridCoordinates_t', parent=z)
            cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
            pu = Internal.getNodeFromName1(cont, 'u')
            pv = Internal.getNodeFromName1(cont, 'v')
            px = Internal.copyNode(pu)
            px[0] = 'CoordinateX'
            py = Internal.copyNode(pv)
            py[0] = 'CoordinateY'
            pz = Internal.copyNode(pv)
            pz[0] = 'CoordinateZ'
            pz[1].ravel('k')[:] = 0.
            gridInit[2] = [px,py,pz]

        # build defTree
        import Ael.Quantum as KDG
        DeformationArgs={
            "Approach"          :  "Quaternions",
            "Epsilon"           :  0.15,
            "Leafsize"          :  8,
            "OmpAllInOne"       :  True,
            "Ndivision"         :  100,
            "NullDisplacements" :  "Weighted",
            "Smoothing"         :  False,
            "printLevel"        :  0 }
        self.DefTree = {}
        for i, z in enumerate(self.RefMeshUV2):
            self.DefTree[i+1] = KDG.KeDefGrid(z, **DeformationArgs)
            self.DefTree[i+1].set_Amplitude(1.)

        return self.RefMesh

    def Dmesh(self):
        """Mesh by deformation."""
        # build new edges in x,y
        nedges = []
        import Geom
        for c, e in enumerate(self.refEdges):
            s = Geom.getCurvilinearAbscissa(e)
            e = Converter.rmVars(e, ['u'])
            e = Converter.addVars([e, s])
            b = OCC.occ.meshOneEdge(self.hook, c+1, -1, -1, -1, -1, e)
            nedges.append(b)

        zones = Internal.getZones(self.RefMeshUV)
        for i, z in enumerate(zones): # pour chaque face
            self.dx1[i+1][:] = 0.
            self.dy1[i+1][:] = 0.
            self.dz1[i+1][:] = 0.

            wires = OCC.occ.meshEdgesOfFace(self.hook, i+1, nedges)
            for c, loops in enumerate(wires):
                for d, w in enumerate(loops):
                    # switch edge to u,v
                    b = w[1]
                    b[0,:] = b[3,:]
                    b[1,:] = b[4,:]
                    b[2,:] = 0.
                    bo = self.wires[i+1][c][d][1] # already in uv
                    inds = self.inds[i+1][c][d]
                    inds2 = self.inds2[i+1][c][d]

                    self.dx1[i+1][inds[:]-1] = b[0,:] - bo[0,:]
                    self.dy1[i+1][inds[:]-1] = b[1,:] - bo[1,:]
                    self.dz1[i+1][inds[:]-1] = b[2,:] - bo[2,:]
                    self.dx1[i+1][inds2[:]-1] = b[0,:] - bo[0,:]
                    self.dy1[i+1][inds2[:]-1] = b[1,:] - bo[1,:]
                    self.dz1[i+1][inds2[:]-1] = b[2,:] - bo[2,:]

            # set displacement on surfaces (zoneName#bcName)
            self.DefTree[i+1].setBndSurfTo("%s#wall"%z[0], "imposed", [self.dx1[i+1],self.dy1[i+1],self.dz1[i+1]]) # [dx,dy,dz]

        # deform
        for i, z in enumerate(self.RefMeshUV2):
            self.DefTree[i+1].makeSources()
            self.DefTree[i+1].computeMeshDisplacement()

        # copie Displacement#0/DisplacementX dans UV
        zones1 = Internal.getZones(self.RefMeshUV)
        zones2 = Internal.getZones(self.RefMeshUV2)
        zones3 = Internal.getZones(self.RefMesh)
        for i, z1 in enumerate(zones3):
            cont = Internal.getNodeFromName1(z1, "GridCoordinates#Init")
            XA1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
            XA2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
            XA3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]
            cont = Internal.getNodeFromName1(z1, "GridCoordinates")
            XB1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
            XB2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
            XB3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]

            z2 = zones2[i]
            defcont = Internal.getNodeFromName1(z2, "Displacement#0")
            DA1 = Internal.getNodeFromName1(defcont, "DisplacementX")[1]
            DA2 = Internal.getNodeFromName1(defcont, "DisplacementY")[1]
            DA3 = Internal.getNodeFromName1(defcont, "DisplacementZ")[1]

            XB1[:] = XA1[:] + DA1[:len(DA1)//2]
            XB2[:] = XA2[:] + DA2[:len(DA2)//2]
            XB3[:] = XA3[:] + DA3[:len(DA3)//2]

            a = C.getFields(Internal.__GridCoordinates__, z1, api=3)[0]
            o = OCC.occ.evalFace(self.hook, a, i+1)
            zones3[i] = C.setFields([o], zones3[i], 'nodes')

            C.convertPyTree2File(self.RefMeshUV2, 'out2.cgns')

        return self.RefMesh

def Loft(name="loft", listSketches=[], listGuides=[], close=True, h=None):
    """Create a loft surface from sketches."""
    return Surface(name=name, listSketches=listSketches, listSketches2=listGuides,
                   type="loft", data={'close':close}, h=h)

def LoftSet(name="loftset", listSketches=[], listGuides=[], close=True, h=None):
    """Create a set of loft surfaces from sketches."""
    return Surface(name=name, listSketches=listSketches, listSketches2=listGuides,
                   type="loft", data={'close':close}, h=h)

def Revolve(name='revolve', sketch=None, center=(0,0,0), axis=(0,0,1), angle=360., h=None):
    """Create a revolution surface from a sketch."""
    return Surface(name=name, listSketches=[sketch],
                   data={'center':center, 'axis':axis, 'angle':angle},
                   type="revolve", h=h)

def Merge(name="compound", listSurfaces=[], h=None):
    """Create a compound surface from a list of surfaces."""
    return Surface(name=name, listSurfaces=listSurfaces,
                   type="merge", h=h)

def MergeEdges(name="mergeEdges", listSketches=[], h=None):
    """Merge edges. Not a surface."""
    return Surface(name=name, listSketches=listSketches, type="mergeEdges", h=h)

def Fill(name="fill", sketch=None, continuity=0, h=None):
    """Create a surface that fill a sketch."""
    return Surface(name=name, listSketches=[sketch],
                   data={'continuity':continuity},
                   type="fill", h=h)

def Union(name="union", listSurfaces1=[], listSurfaces2=[], h=None):
    """Boolean union."""
    return Surface(name=name, listSurfaces=listSurfaces1,
                   listSurfaces2=listSurfaces2,
                   type="union", h=h)

def Sub(name="sub", listSurfaces1=[], listSurfaces2=[], h=None):
    """Boolean difference."""
    return Surface(name=name, listSurfaces=listSurfaces1,
                   listSurfaces2=listSurfaces2,
                   type="sub", h=h)

def Inter(name="inter", listSurfaces1=[], listSurfaces2=[], h=None):
    """Boolean intersection."""
    return Surface(name=name, listSurfaces=listSurfaces1,
                   listSurfaces2=listSurfaces2,
                   type="inter", h=h)

#============================================================
class Volume2D():
    """Define a parametric 2D volume."""
    def __init__(self, name=None, listSketches=[], orders=[]):
        # name
        if name is not None: self.name = name
        else: self.name = getName("vol")
        # sketches that define the bounded volume
        self.sketches = listSketches
        # optional ordering of sketches
        self.orders = orders
        # reference mesh for Dmesh
        self.RefMesh = None # mesh (pytree)
        self.RefMesh2 = None # mesh with kplane (pytree)
        self.refBorders = None # borders of ref mesh (arrays)
        self.inds = None # indices of borders in refMesh BC
        self.inds2 = None # indices of borders in refMesh BC with kplane
        self.DefTree = None # def tree (pytree)

    def mesh(self):
        """Call the volume mesher."""
        import Generator, Transform
        # call sketch mesher
        borders = []
        for c, s in enumerate(self.sketches):
            m = s.mesh()
            if c < len(self.orders) and self.orders[c] == -1: m = Transform.reorder(m, (-1,2,3))
            borders += m
        borders = Converter.convertArray2Tetra(borders)
        borders = Transform.join(borders)
        m = Generator.T3mesher2D(borders, triangulateOnly=0, grading=1.1, metricInterpType=0)
        m = Transform.reorder(m, (1,))
        return m

    def Mesh(self):
        m = self.mesh()
        z = Internal.createZoneNode(self.name, m, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        return z

    def MeshAsReference(self):
        """Mesh and store the mesh for future Dmesh."""
        # keep borders
        self.refBorders = []
        for s in self.sketches: self.refBorders.append(s.mesh())

        # build refMesh
        self.RefMesh = self.Mesh()
        import Transform.PyTree as T
        import Transform
        import Post.PyTree as P
        a = C.convertArray2NGon(self.RefMesh, recoverBC=False)
        p = P.exteriorFaces(a)
        C._addBC2Zone(a, 'wall', 'BCWall', subzone=p)
        a = T.addkplane(a)
        a = T.contract(a, (0,0,0), (1,0,0), (0,1,0), 0.1)
        self.RefMesh2 = a
        R._copyGrid2GridInit(self.RefMesh, mode=1)

        # identify borders in bc
        bc = C.extractBCOfName(self.RefMesh2, "wall", reorder=False)[0]
        xc = Internal.getNodeFromName2(bc, "CoordinateX")
        self.dx1 = numpy.empty((xc[1].size), dtype=numpy.float64)
        self.dy1 = numpy.empty((xc[1].size), dtype=numpy.float64)
        self.dz1 = numpy.empty((xc[1].size), dtype=numpy.float64)
        a = C.getAllFields(bc, 'nodes', api=3)[0]
        hook = Converter.createHook(a, function='nodes')
        self.inds = []
        for borders in self.refBorders:
            out = []
            for b in borders:
                nodes = Converter.identifyNodes(hook, b)
                out.append(nodes)
            self.inds.append(out)
        self.inds2 = []
        for borders in self.refBorders:
            out = []
            for b in borders:
                b = Transform.translate(b, (0,0,0.1))
                nodes = Converter.identifyNodes(hook, b)
                out.append(nodes)
            self.inds2.append(out)
        Converter.freeHook(hook)

        # build defTree
        import Ael.Quantum as KDG
        DeformationArgs={
            "Approach"          :  "Quaternions",
            "Epsilon"           :  0.15,
            "Leafsize"          :  8,
            "OmpAllInOne"       :  True,
            "Ndivision"         :  100,
            "NullDisplacements" :  "Weighted",
            "Smoothing"         :  False,
            "printLevel"        :  0 }
        self.DefTree = KDG.KeDefGrid(self.RefMesh2, **DeformationArgs)
        self.DefTree.set_Amplitude(1.)
        return self.RefMesh

    def Dmesh(self):
        """Mesh by deformation."""
        self.dx1[:] = 0.; self.dy1[:] = 0.; self.dz1[:] = 0.

        # call sketch mesher on borders
        for c, s in enumerate(self.sketches):
            out = s.rmesh(self.refBorders[c])
            for d, b in enumerate(out):
                b = out[d][1]
                bo = self.refBorders[c][d][1]
                inds = self.inds[c][d]
                inds2 = self.inds2[c][d]
                self.dx1[inds[:]-1] = b[0,:] - bo[0,:]
                self.dy1[inds[:]-1] = b[1,:] - bo[1,:]
                self.dz1[inds[:]-1] = b[2,:] - bo[2,:]
                self.dx1[inds2[:]-1] = b[0,:] - bo[0,:]
                self.dy1[inds2[:]-1] = b[1,:] - bo[1,:]
                self.dz1[inds2[:]-1] = b[2,:] - bo[2,:]
        # set displacement on surfaces (zoneName#bcName)
        self.DefTree.setBndSurfTo("%s#wall"%self.name, "imposed", [self.dx1,self.dy1,self.dz1]) # [dx,dy,dz]
        # deform
        self.DefTree.makeSources()
        self.DefTree.computeMeshDisplacement()

        # copie Displacement#0/DisplacementX dans Coordinates
        zones = Internal.getZones(self.RefMesh)
        z1 = zones[0]
        cont = Internal.getNodeFromName1(z1, "GridCoordinates#Init")
        XA1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
        XA2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
        XA3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]
        cont = Internal.getNodeFromName1(z1, "GridCoordinates")
        XB1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
        XB2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
        XB3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]

        zones = Internal.getZones(self.RefMesh2)
        z2 = zones[0]
        defcont = Internal.getNodeFromName1(z2, "Displacement#0")
        DA1 = Internal.getNodeFromName1(defcont, "DisplacementX")[1]
        DA2 = Internal.getNodeFromName1(defcont, "DisplacementY")[1]
        DA3 = Internal.getNodeFromName1(defcont, "DisplacementZ")[1]

        XB1[:] = XA1[:] + DA1[:len(DA1)//2]
        XB2[:] = XA2[:] + DA2[:len(DA2)//2]
        XB3[:] = XA3[:] + DA3[:len(DA3)//2]

        return self.RefMesh

#============================================================
class Volume3D():
    """Define a parametric 3D volume."""
    def __init__(self, name=None, listSurfaces=[], orders=[]):
        # name
        if name is not None: self.name = name
        else: self.name = getName("vol")
        # surfaces that define the bounded volume
        self.surfaces = listSurfaces
        # optional ordering of surfaces
        self.orders = orders

    def mesh(self):
        """Call the volume mesher."""
        import Generator, Transform
        # call surface mesher
        borders = []
        for c, s in enumerate(self.surfaces):
            m = s.mesh()
            if c < len(self.orders) and self.orders[c] == -1: m = Transform.reorder(m, (-1,))
            borders += m
        borders = Transform.join(borders)
        m = Generator.TetraMesher(borders)
        return m

    def Mesh(self):
        m = self.mesh()
        z = Internal.createZoneNode(self.name, m, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        return z

#============================================================
class Eq:
    """Equation"""
    def __init__(self, expr1, expr2=None):
        # references sur l'equation sympy
        self.s = sympy.Eq(expr1, expr2)
        DRIVER.registerEquation(self)

    def analyse(self):
        """Analyse equation to return vars and symbols."""
        keywords = ["=", "length", "\*", "\+", "\-", "\/", "cos", "sin", "\(", "\)"]
        pattern = f"({'|'.join(keywords)})"
        segments = re.split(pattern, self.expr)
        segments = [s.strip() for s in segments if s]  # nettoyage
        # replace by their id
        out = ''; vars = []
        for s in segments:
            if s in DRIVER.scalars:
                id = DRIVER.scalars[s].name
                out += id
                vars.append(id)
            else: out += s
        return vars, out

#============================================================
class Lt:
    """Constraint inequation"""
    def __init__(self, expr1, expr2=None):
        # references sur l'inequation sympy
        self.s = sympy.Lt(expr1, expr2)
        DRIVER.registerInequation(self)

#============================================================
class Le:
    """Constraint inequation"""
    def __init__(self, expr1, expr2=None):
        # references sur l'inequation sympy
        self.s = sympy.Le(expr1, expr2)
        DRIVER.registerInequation(self)

#============================================================
class Driver:
    """Driver is Model"""
    def __init__(self):
        # all parameters
        self.scalars = {} # id -> scalar
        self.points = {} # points
        self.grids = {} # grids
        # db
        self.db = None # optional db
        # all entities
        self.edges = {} # edges
        self.sketches = {} # wires
        self.surfaces = {} # shapes
        # all equations
        self.equationCount = 0
        self.equations = {} # equations
        # all inequations
        self.inequationCount = 0
        self.inequations = {} # inequations
        # updated by solve
        self.solution = None # solution of system in sympy symbols
        self.params = None # all model params in sympy symbols
        self.freeParams = None # all model free params in sympy symbols (set order)
        # DOE
        self.doeFileName = 'doe.hdf'
        self.doeRange = [] # list free param->discretized range
        self.doeSize = [] # list free param->size of discretization
        self.iter = None # itertools iterator on DOE
        # ROM
        self.romFileName = 'rom.hdf'
        self.K = 0 # reduced dimension
        self.mean = None
        self.Phi = None # POD vectors
        self.ak = None # POD coefficients of each mesh in doe

    def registerScalar(self, s):
        """Register parametric scalar."""
        self.scalars[s.name] = s # name -> scalar

    def registerPoint(self, p):
        """Register parametric point."""
        self.points[p.name] = p

    def registerGrid(self, p):
        """Register parametric grid."""
        self.grids[p.name] = p

    def registerEdge(self, e):
        """Register parametric entity."""
        self.edges[e.name] = e

    def registerSketch(self, e):
        """Register parametric sketch."""
        self.sketches[e.name] = e

    def registerSurface(self, e):
        """Register parametric surface."""
        self.surfaces[e.name] = e

    def registerEquation(self, eq):
        """Register equation."""
        self.equations["EQUATION%04d"%self.equationCount] = eq
        self.equationCount += 1
        # all concerned Scalar are tagged as free parameters
        symbols = eq.s.free_symbols
        for s in symbols:
            scalar = self.scalars[s.name]
            if scalar.range is None:
                scalar.range = [-999.99, 999.99] # range ajustable

    def registerInequation(self, eq):
        """Register inequation."""
        self.inequations["INEQUATION%04d"%self.inequationCount] = eq
        self.inequationCount += 1
        # all concerned Scalar are tagged as free parameters
        symbols = eq.s.free_symbols
        for s in symbols:
            scalar = self.scalars[s.name]
            if scalar.range is None:
                scalar.range = [-999.99, 999.99] # range ajustable

    def print(self):
        """Print registered entities."""
        for k in self.scalars: print(k)
        for k in self.points: print(k)
        for k in self.edges: print(k)
        for k in self.sketches: print(k)
        for k in self.surfaces: print(k)
        for k in self.equations: print(k)

    def update(self):
        """Update all entities from parameters."""
        for k in self.edges: self.edges[k].update()
        for k in self.sketches: self.sketches[k].update()
        for k in self.surfaces: self.surfaces[k].update()

    def solve(self):
        """Solve equations to get free parameters."""
        # get params
        params = []
        for s in self.scalars:
            mu = self.scalars[s]
            if mu.isFree(): params.append(mu)
        params.reverse() # reverse order to solve for explicit variables
        print('SOLVE: params=', params)

        # get equations, sub fixed params
        equations = []
        for e in self.equations:
            eq = self.equations[e]
            equations.append(eq.s)
            for s in self.scalars:
                mu = self.scalars[s]
                if not mu.isFree(): eq.s.subs(mu, mu.v)
        print('SOLVE: eqs=', equations)

        # solve([eq0,eq1], [x0,x1])
        solution = sympy.solve(equations, params, dict=True)
        print('SOLVE: sol=', solution)
        if len(solution) == 0:
            print('SOLVE: no solution')
        elif len(solution) > 1:
            print('SOLVE: many solutions, taking first')
            solution = solution[0]
        else:
            solution = solution[0]

        # number of free vars
        nparams = len(params)
        neqs = len(equations)
        nd = nparams - neqs
        print("SOLVE: nparams=", nparams)
        print("SOLVE: neqs=", neqs)
        print("SOLVE: free params=", nd)

        # who is free and valid at the end?
        freeParams = params[:]
        for s in solution:
            if solution[s].is_Float or solution[s].is_Integer or solution[s].is_Rational:
                print('SOLVE: fixed', s, 'to', solution[s])
                self.scalars[s.name].v = solution[s]
                if self.scalars[s.name].check(): print('=> valid')
                else: print('=> invalid')
            freeParams.remove(s)
        print('SOLVE: free vars=', freeParams)

        self.solution = solution
        self.params = params
        self.freeParams = freeParams

        return solution, freeParams

    # instantiation of free parameters
    # IN: paramValues: dict of free parameters given values
    # OUT: return True if valid, False if invalid
    def instantiate(self, paramValues):
        """Instantiate all from given paramValues."""
        valid = True # return
        # set free params from input
        error = False
        for f in self.freeParams:
            if f.name not in paramValues:
                print("Error: instantiate: you should specify: ", f.name)
                error = True
            else:
                self.scalars[f.name].v = paramValues[f.name]
                print('SET: fixed', f, '=', paramValues[f.name])
                if self.scalars[f.name].check(): print('SET: => valid')
                else: print('SET: => invalid'); valid = False

        if error: raise ValueError("instantiate: stopping.")

        # set other vars with equations
        soli = self.solution.copy()
        for s in soli:
            for f in self.freeParams:
                soli[s] = soli[s].subs(f, paramValues[f.name])

        # check validity for ranges
        for s in soli:
            if soli[s].is_Float or soli[s].is_Integer or soli[s].is_Rational:
                print('SET: fixed', s, '=', soli[s])
                self.scalars[s.name].v = soli[s]
                if self.scalars[s.name].check(): print('SET: => valid')
                else: print('SET: => invalid'); valid = False
            else: print('SET: some variables were not instantiated'); valid = False

        # Check validity for inequations
        params = {}
        for f in paramValues: params[self.scalars[f]] = paramValues[f]
        for s in soli:
            if soli[s].is_Float or soli[s].is_Integer or soli[s].is_Rational:
                params[s] = self.scalars[s.name].v

        for c, e in enumerate(self.inequations):
            ret = self.inequations[e].s.subs(params)
            #ret = self.inequations[e].s.evalf()
            if ret: print('SET: => ineq %d is valid'%c)
            else: print("SET: => ineq %d is invalid"%c); valid = False

        # update geometries
        self.update()

        # return True if valid in range and inequation constraints
        return valid

    # FD of free parameters on discrete mesh
    # if freevars is None, derivate for all free parameters else derivate for given parameters
    def _dXdmu(self, entity, mesh, freeParams=None, deps=1.e-10):
        """Compute derivatives dX/dmu on entity."""
        import KCore

        if len(self.freeParams) == 0:
            print("Warning: no free vars.")
            return None

        if freeParams is None: # no param given
            listVars = self.freeParams
        elif isinstance(freeParams, str): # free param given by name
            listVars = []
            for f in self.freeParams:
                if self.scalars[f.name].name == freeParams: listVars.append(f)
        elif isinstance(freeParams, list): # suppose list of names
            listVars = []
            for f in self.freeParams:
                if self.scalars[f.name].name in freeParams: listVars.append(f)
        else:
            raise TypeError("Warning: dXdmu: incorrect freevars.")

        mesho = Converter.copy(mesh)

        for c, f in enumerate(listVars):
            # free vars value dict
            d = {}
            for q in self.freeParams:
                d[q.name] = self.scalars[q.name].v
            d[f.name] += deps

            print("DIFF on: ", f.name)

            # update CAD at param+eps
            self.instantiate(d)
            mesho = entity.rmesh(mesh)
            # get derivatives
            Converter._addVars(mesh, ['dXd%d'%c, 'dYd%d'%c, 'dZd%d'%c])
            for p, m in enumerate(mesh):
                pos1 = KCore.isNamePresent(m, 'dXd%d'%c)
                pos2 = KCore.isNamePresent(m, 'dYd%d'%c)
                pos3 = KCore.isNamePresent(m, 'dZd%d'%c)
                p1x = m[1]
                p2x = mesho[p][1]
                p1x[pos1,:] = (p2x[0,:]-p1x[0,:])/deps
                p1x[pos2,:] = (p2x[1,:]-p1x[1,:])/deps
                p1x[pos3,:] = (p2x[2,:]-p1x[2,:])/deps

        # remet le hook original
        d = {}
        for q in self.freeParams:
            d[q.name] = self.scalars[q.name].v
        self.instantiate(d)

        return None

    # connect driver to db
    def connect(self, db):
        """Connect driver to db."""
        self.db = db

    # set DOE deltas for free parameters
    # It is better to set them in scalar.range
    # IN: dict of deltas for each desired free parameter
    # OUT: arange dict
    def setDOE(self, deltas={}):
        # set default
        self.doeRange = []; self.doeSize = []
        for f in self.freeParams: # give order
            p = self.scalars[f.name]
            if len(p.range) == 3: # disc given in range
                mean = 0.5*(p.range[1] - p.range[0])
                tol = 1.e-11 + 1.e-14*mean
                self.doeRange.append(numpy.arange(p.range[0], p.range[1]+tol, p.range[2]))
                self.doeSize.append(p.range[2])
            else: # set to 2 points
                self.doeRange.append(numpy.linspace(p.range[0], p.range[1], 2))
                self.doeSize.append(2)

        # set dictionary (optional)
        for k in deltas: # free param names
            for c, f in enumerate(self.freeParams):
                if self.scalars[f.name].name == k:
                    p = self.scalars[f.name]
                    mean = 0.5*(p.range[1] - p.range[0])
                    tol = 1.e-11 + 1.e-14*mean
                    self.doeRange[c] = numpy.arange(p.range[0], p.range[1]+tol, deltas[k])
                    self.doeSize[c] = deltas[k]
        return None

    # walk DOE (iterateur), instantiate
    # return the next valid free parameters point
    def walkDOE(self):
        if self.iter is None:
            # set range
            self.setDOE()
            # create iterator
            ranges = []; size = 0
            for k in self.doeRange:
                ranges.append(k)
                size += k.size
            self.iter = itertools.product(*ranges)
        # iterate
        try:
            p = next(self.iter)
        except:
            return None # end of DOE
        # compute parametric point
        pt = {}
        for c, s in enumerate(self.freeParams):
            pt[self.scalars[s.name].name] = p[c]
        # instantiate
        st = "DOE: Checking point: { "
        for k in pt: st += k+' = %g '%pt[k]
        st += ' }'
        print(st, flush=True)
        valid = self.instantiate(pt)
        if valid:
            if self.db is not None:
                exist = self.db.exist(pt)
                if not exist: return pt
                else:
                    print("DOE: => Already in db. Skipped.")
                    return self.walkDOE()
            else: return pt
        else: return self.walkDOE()

    # walk DOE1, instantiate, parallel CFD but sequential on parameters
    def walkDOE1(self):
        if Cmpi.rank == 0:
            pt = self.walkDOE()
        else:
            pt = self.walkDOE()
        return pt

    # walk DOE2, instantiate, parallel tasks (no on proc 0)
    def walkDOE2(self):
        if Cmpi.rank == 0:
            if Cmpi.size > 1:
                free = Cmpi.recv()
                pt = self.walkDOE()
                if pt is None:
                    for i in range(1,Cmpi.size):
                        Cmpi.isend(None, dest=i)
                    return None
                else:
                    Cmpi.isend(pt, dest=free)
                    return 1
            else:
                pt = self.walkDOE()
                return pt

        else:
            Cmpi.isend(Cmpi.rank, dest=0, tag=1) # i am free
            pt = Cmpi.recv(source=0) # wait for task
        return pt

    # walk DOE, instantiate, mesh, append snapshots to file, parallel
    def walkDOE3(self, entity):
        self.setDOE()
        ranges = []; size = 0
        for k in self.doeRange:
            ranges.append(range(k.size))
            size += k.size
        raf = size - size%Cmpi.size # seq reste a faire

        for indexes in itertools.product(*ranges):
            # create value dict
            values = {}
            hash = self.getHash(indexes)
            for c, f in enumerate(self.freeParams):
                val = self.doeRange[c][indexes[c]]
                f = self.freeParams[c]
                p = self.scalars[f.name]
                values[p.name] = val
            if Cmpi.rank == hash%Cmpi.size and hash < raf:
                # instantiate and mesh
                self.instantiate(values)
                mesh = entity.mesh()
                if Cmpi.rank == 0:
                    self.addSnapshot(hash, mesh)
                    if Cmpi.size > 1: Cmpi.send(1, dest=1)
                else:
                    go = Cmpi.recv(source=Cmpi.rank-1)
                    self.addSnapshot(hash, mesh)
                    if Cmpi.rank < Cmpi.size-1: Cmpi.send(Cmpi.rank+1, dest=Cmpi.rank+1)
                Cmpi.barrier()
            elif Cmpi.rank == 0 and hash >= raf:
                # instantiate and mesh
                self.instantiate(values)
                mesh = entity.mesh()
                self.addSnapshot(hash, mesh)

    # IN: list of indexes (i,j,k,...) one for each param
    # OUT: single hash integer (flatten)
    def getHash(self, indexes):
        hash = 0
        for c, i in enumerate(indexes):
            if c == 0: hash = i
            else: hash += i*self.doeSize[c-1]
        return hash

    # IN: hash
    # OUT: return (i,j,k,...)
    def getInd(self, hash):
        hashcode = hash
        np = len(self.doeSize)
        out = []
        for c in range(np):
            prod = 1
            for s in self.doeSize[:np-c]: prod *= s
            h = hashcode // prod
            out.append(h)
            hashcode = hashcode - h*prod
        return out

    # remesh input mesh to match nv points using refine
    # used nly in snapshots read (obsolete)
    def remesh(self, mesh, nv):
        import Generator
        nm = mesh[1].shape[1]
        if nm != nv:
            power = (nv-1)*1./(nm-1)
            m = Generator.refine(mesh, power, dir=1)
        else: m = mesh
        return m

    # DOE in file (to be replaced by DB)
    def createDOE(self, fileName):
        self.doeFileName = fileName
        if Cmpi.rank > 0: return None
        t = C.newPyTree(['Parameters', 'Snapshots'])
        for c, k in enumerate(self.doeRange):
            node = ["p%05d"%c, k, [], 'parameter_t']
            t[2][1][2].append(node)
        C.convertPyTree2File(t, self.doeFileName)
        return None

    # add snapshot to file (to be replaced by DB)
    def addSnapshot(self, hashcode, msh):
        import Converter.Distributed as Distributed
        import Transform
        msh = Converter.extractVars(msh, ['x','y','z'])
        msh = Transform.join(msh) # merge mesh
        node = ["%05d"%hashcode, msh[1], [], 'snapshot_t']
        Distributed.writeNodesFromPaths(self.doeFileName, 'CGNSTree/Snapshots', node)
        print("ADD: snapshot %d added."%hashcode, flush=True)
        return None

    # read a snapshot, return a mesh array (to be replaced by DB)
    def readSnaphot(self, hashcode):
        import Converter.Distributed as Distributed
        nodes = Distributed.readNodesFromPaths(self.doeFileName, ['CGNSTree/Snapshots/%05d'%hashcode])
        msh = nodes[0][1]
        msh = ['x,y,z', msh, msh.shape[1], 1, 1]
        return msh

    # read all snapshots, return a flatten matrix
    # IN: paramSlab: optional ( (5,120), (3,20), ...) for each parameters
    def readAllSnapshots(self, paramSlab=None):
        ranges = []; np = 0
        if paramSlab is None: # full range
            for k in self.doeRange:
                ranges.append(range(k.size))
                np += k.size
        else: # given range
            for k in paramSlab:
                ranges.append(range(k[0],k[1]))
                np += k[1]-k[0]
        m = self.readSnaphot(0)
        nv = m[1].shape[1]
        F = numpy.empty( (nv*3, np), dtype=numpy.float64)

        for indexes in itertools.product(*ranges):
            hash = self.getHash(indexes)
            m = self.readSnaphot(hash)
            m = self.remesh(m, nv)
            m = m[1].ravel('k')
            F[:,hash] = m[:]
        return F

    # ROM (Model)
    def writeROM(self, fileName):
        self.romFileName = fileName
        if Cmpi.rank > 0: return None
        t = C.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'phi_t']
        t[2][1][2].append(node)
        C.convertPyTree2File(t, self.romFileName)
        return None

    # build POD from full matrix F, keep K modes
    def createROM(self, F, K=-1):
        # on deformation from
        mean = numpy.mean(F, axis=1, keepdims=True)
        self.mean = mean # maillage moyen
        #F = F - mean
        Phi, S, Vt = numpy.linalg.svd(F, full_matrices=False)
        # energy of each modes
        #energy = S**2 / numpy.sum(S**2)
        if K > 0: self.K = K
        else: self.K = Phi.shape[1]
        self.Phi = Phi[:, 0:self.K]
        return Phi, S, Vt

    # get modes as meshes
    def getMode(self, i):
        m = self.Phi[:,i]
        np = m.size//3
        m = m.reshape( (3,np) )
        m = ['x,y,z', m, np, 1, 1]
        return m

    # get coords of mesh on POD
    def addCoefs(self, hashcode, msh):
        import Converter.Distributed as Distributed
        coords = numpy.empty( (self.K), dtype=numpy.float64 )
        m = msh[1].ravel('k')
        for i in range(self.K):
            c0 = numpy.dot(self.Phi[:,i], m)
            coords[i] = c0
        node = ["%05d"%hashcode, coords, [], 'coeffs_t']
        Distributed.writeNodesFromPaths(self.romFileName, 'CGNSTree/POD', node)
        print("ADD: coeffs %d added."%hashcode)
        return coords

    def addAllCoefs(self):
        ranges = []; size = 0
        for k in self.doeRange:
            ranges.append(range(k.size))
            size += k.size
        raf = size - size%Cmpi.size # seq reste a faire

        m = self.readSnaphot(0)
        nv = m[1].shape[1]

        for indexes in itertools.product(*ranges):
            hash = self.getHash(indexes)
            if Cmpi.rank == hash%Cmpi.size and hash < raf:
                m = self.readSnaphot(hash)
                m = self.remesh(m, nv)
                if Cmpi.rank == 0:
                    self.addCoefs(hash, m)
                    if Cmpi.size > 1: Cmpi.send(1, dest=1)
                else:
                    go = Cmpi.recv(source=Cmpi.rank-1)
                    self.addCoefs(hash, m)
                    if Cmpi.rank < Cmpi.size-1: Cmpi.send(Cmpi.rank+1, dest=Cmpi.rank+1)
                Cmpi.barrier()
            elif Cmpi.rank == 0 and hash >= raf:
                m = self.readSnaphot(hash)
                m = self.remesh(m, nv)
                self.addCoefs(hash, m)

    def evalROM(self, coords):
        m = self.Phi @ coords
        m = m.reshape((3, m.size//3))
        msh = ['x,y,z', m, m.shape[1], 1, 1]
        return msh

    def readCoefs(self, hashcode):
        import Converter.Distributed as Distributed
        nodes = Distributed.readNodesFromPaths(self.romFileName, ['CGNSTree/POD/%05d'%hashcode])
        coord = nodes[0][1]
        return coord

    # rebuild samples from POD
    def rebuildF(self, Phi, S, Vt):
        # Convert S to a diagonal matrix
        Sigma = numpy.diag(S)
        # Multiply to get back A
        Fr = Phi @ Sigma @ Vt
        return Fr

#============================================================
# Global
DRIVER=Driver()
