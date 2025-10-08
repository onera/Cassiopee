# Parametric CAD driver
import OCC
import sympy
import numpy
import re
import Converter.Mpi as Cmpi

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
class Scalar:
    """Define a parametric scalar"""
    def __init__(self, T, name=None):
        # name
        if name is not None: self.name = name
        else: self.name = getName("scalar")
        # Id pour unicite
        #self.id = getName("Id")
        self.id = name

        # symbol sympy
        self.s = sympy.Symbol(self.id)
        # instantiated value
        self.v = T
        # range
        self.range = None
        # register
        DRIVER.registerScalar(self)

    # return value
    def v(self):
        return self.v

    # print content
    def print(self, shift=0):
        print(" "*shift, "id", self.id)
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
    def __init__(self, T, name=None):
        # name
        if name is not None: self.name = name
        else: self.name = getName("vec2")
        # parameters
        self.x = Scalar(T[0], name="x")
        self.y = Scalar(T[1], name="y")

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
    def __init__(self, T, name=None):
        # name
        if name is not None: self.name = name
        else: self.name = getName("point")
        # parameters
        self.x = Scalar(T[0], name=self.name+".x")
        self.y = Scalar(T[1], name=self.name+".y")
        self.z = Scalar(T[2], name=self.name+".z")
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
    def __init__(self, Xo, Xf, N, name=None):
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
        self.P = [[[Point((self.xo+i*self.dx, self.yo+j*self.dy, self.zo+k*self.dz), name="P.%d.%d.%d"%(i,j,k)) for k in range(self.nk)] for j in range(self.nj)] for i in range(self.ni)]
        # register
        DRIVER.registerGrid(self)

    # return value
    def v(self, T):
        (i,j,k) = T
        return (self.P[i,j,k].x.v, self.P[i,j,k].y.v, self.P[i,j,k].z.v)

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
    def __init__(self, listP, type=None, name=None, mesh=None):
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
                self.P.append(Scalar(P, name=name+'.P%d'%c))
            elif isinstance(P, int):
                self.P.append(Scalar(P), name=name+'.P%d'%c)
            elif isinstance(P, tuple) and len(P) == 3:
                self.P.append(Point(P, name=name+'.P%d'%c))
            elif isinstance(P, Point):
                self.P.append(P)
            elif isinstance(P, tuple) and len(P) == 2:
                self.P.append(Vec2(P, name=name+'.P%d'%c))
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
        self.hook = OCC.occ.createEmptyCAD("unknown.stp", "fmt_step")
        if self.type == "line":
            OCC.occ.addLine(self.hook, self.P[0].v(), self.P[1].v())
        elif self.type == "spline1": # by control points
            s = len(self.P)
            n = numpy.zeros((3,s), dtype=numpy.float64)
            for c, p in enumerate(self.P): n[:,c] = p.v()
            OCC.occ.addSpline(self.hook, n, 0)
        elif self.type == "spline2": # by approximated points
            s = len(self.P)
            n = numpy.zeros((3,s), dtype=numpy.float64)
            for c, p in enumerate(self.P): n[:,c] = p.v()
            OCC.occ.addSpline(self.hook, n, 1)
        elif self.type == "spline3": # by free form control points + mesh
            # self.P[0] is a Grid, mesh is an array

            # rebuild control points
            import Generator, Transform, Converter
            grid = self.P[0]
            ni, nj, nk = grid.ni, grid.nj, grid.nk
            cp = Generator.cart((0,0,0), (1,1,1), (ni,nj,nk))
            cp = Converter.addVars(cp, ['dx','dy','dz'])

            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        ind = i+j*ni+k*ni*nj
                        cp[1][0,ind] = grid.P[i][j][k].x.v
                        cp[1][1,ind] = grid.P[i][j][k].y.v
                        cp[1][2,ind] = grid.P[i][j][k].z.v
                        cp[1][3,ind] = grid.P[i][j][k].x.v - grid.xo - i*grid.dx
                        cp[1][4,ind] = grid.P[i][j][k].y.v - grid.yo - j*grid.dy
                        cp[1][5,ind] = grid.P[i][j][k].z.v - grid.zo - k*grid.dz

            mesh = Transform.freeForm(self.mesh, cp)
            mesh = Transform.deform(mesh, ['dx','dy','dz'])

            # spline from approximated point
            OCC.occ.addSpline(self.hook, mesh[1], 1)

        elif self.type == "circle":
            OCC.occ.addCircle(self.hook, self.P[0].v(), (0,0,1), self.P[1].v, 0)
        elif self.type == "arc":
            OCC.occ.addArc(self.hook, self.P[0].v(), self.P[1].v(), self.P[2].v())
        else:
            raise(ValueError, "Unknown entity type %s."%self.type)

    def print(self, shift=0):
        for c, P in enumerate(self.P):
            print(" "*shift, P.name)
            P.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        OCC.occ.writeCAD(self.hook, fileName, format)

    # mesh sketch
    def mesh(self, hmin, hmax, hausd):
        edges = OCC.meshAllEdges(self.hook, hmin, hmax, hausd, -1)
        return edges

    # check parameters validity
    def check(self):
        for P in self.P:
            ret = P.check()
            if ret == 0: return 0
        return 1

#============================================================
# line from two parametric points
def Line(P1, P2, name=None):
    return Entity([P1, P2], type="line", name=name)

# spline from parmaetric control points
def Spline1(CPs, name=None):
    return Entity(CPs, type="spline1", name=name)

# spline from parametric approximated points
def Spline2(Ps, name=None):
    return Entity(Ps, type="spline2", name=name)

# spline from parametric grid
def Spline3(PGrid, mesh, name=None):
    return Entity([PGrid], mesh=mesh, type="spline3", name=name)

# circle from parametric center and radius
def Circle(C, R, name=None):
    return Entity([C, R], type="circle", name=name)

# arc from 3 parametric points
def Arc(P1, P2, P3, name=None):
    return Entity([P1, P2, P3], type="arc", name=name)

#============================================================
class Sketch():
    """Define a parametric sketch from a list of entities."""
    def __init__(self, listEntities=[], name="sketch"):
        # name
        self.name = getName(name)
        # dependant entities
        self.entities = listEntities
        # type
        self.type = "sketch"
        # hook
        self.hook = None
        # global parameters (always added)
        self.P = []
        P = Vec3((0,0,0), name='%s.position'%self.name)
        self.P.append(P)
        self.position = P
        P = Point((0,0,0), name='%s.rotCenter'%self.name)
        self.P.append(P)
        self.rotCenter = P
        P = Vec3((0,0,1), name='%s.rotAxis'%self.name)
        self.P.append(P)
        self.rotAxis = P
        P = Scalar(0., name='%s.rotAngle'%self.name)
        self.P.append(P)
        self.rotAngle = P
        # update hook
        self.update()
        # register
        DRIVER.registerSketch(self)

    def add(self, entity):
        self.entities.append(entity)

    # update the CAD from parameters
    def update(self):
        if self.hook is not None: OCC.occ.freeHook(self.hook)
        self.hook = OCC.occ.createEmptyCAD('unknown.step', 'fmt_step')
        hooks = []
        for e in self.entities: hooks.append(e.hook)
        self.hook = OCC.occ.mergeCAD(hooks)
        # global positionning
        OCC._rotate(self.hook, self.P[1].v(), self.P[2].v(), self.P[3].v)
        OCC._translate(self.hook, self.P[0].v())

    # check if parameters are valid
    def check(self):
        for e in self.entities:
            ret = e.check()
            if ret == 1: return 1
        return 0

    # print information
    def print(self, shift=0):
        for e in self.entities:
            print(" "*shift, e.name)
            e.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        OCC.occ.writeCAD(self.hook, fileName, format)

    # mesh sketch
    def mesh(self, hmin, hmax, hausd):
        edges = OCC.meshAllEdges(self.hook, hmin, hmax, hausd, -1)
        return edges

#============================================================
class Surface():
    """Define a parametric surface."""
    def __init__(self, listSketches=[], listSurfaces=[], data={}, name="surface", type="loft"):
        # name
        self.name = getName(name)
        # type
        self.type = type
        # dependant sketches
        self.sketches = listSketches
        # dependant surfaces
        self.surfaces = listSurfaces
        # optional data
        self.data = data
        # hook
        self.hook = None
        # global parameters (always added)
        self.P = []
        P = Vec3((0,0,0), name='%s.position'%self.name)
        self.P.append(P)
        self.position = P
        P = Point((0,0,0), name='%s.rotCenter'%self.name)
        self.P.append(P)
        self.rotCenter = P
        P = Vec3((0,0,1), name='%s.rotAxis'%self.name)
        self.P.append(P)
        self.rotAxis = P
        P = Scalar(0., name='%s.rotAngle'%self.name)
        self.P.append(P)
        self.rotAngle = P
        # update hook
        self.update()
        # register
        DRIVER.registerSurface(self)

    def add(self, sketch):
        """Add a sketch to the surface definition."""
        self.sketches.append(sketch)

    # update the CAD from parameters
    def update(self):
        """Update CAD hook from parameters."""
        if self.hook is not None: OCC.occ.freeHook(self.hook)
        self.hook = OCC.occ.createEmptyCAD('unknown.step', 'fmt_step')
        if self.type == "loft":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            OCC.occ.loft(self.hook, [i for i in range(1,len(hooks)+1)], [])
        elif self.type == "revolve":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            OCC.occ.revolve(self.hook, [1], self.data['center'], self.data['axis'], self.data['angle'])
        elif self.type == "compound":
            hooks = []
            for e in self.surfaces: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
        elif self.type == "fill":
            hooks = []
            for e in self.sketches: hooks.append(e.hook)
            self.hook = OCC.occ.mergeCAD(hooks)
            self.hook = OCC.occ.fillHole(self.hook, [i for i in range(1,len(hooks)+1)], [], self.data['continuity'])

        # global positionning
        OCC._rotate(self.hook, self.P[1].v(), self.P[2].v(), self.P[3].v)
        OCC._translate(self.hook, self.P[0].v())

    # print information
    def print(self, shift=0):
        """Print surface information."""
        for e in self.sketches:
            print(" "*shift, e.name)
            e.print(shift+4)

    # export CAD to file
    def writeCAD(self, fileName, format="fmt_step"):
        """Export to CAD file."""
        OCC.occ.writeCAD(self.hook, fileName, format)

    # mesh surface
    def mesh(self, hmin, hmax, hausd):
        """Mesh surface."""
        edges = OCC.meshAllEdges(self.hook, hmin, hmax, hausd, -1)
        nbFaces = OCC.getNbFaces(self.hook)
        faceList = range(1, nbFaces+1)
        if hausd < 0: hList = [(hmax,hmax,hausd)]*len(faceList)
        else: hList = [(hmin,hmax,hausd)]*len(faceList)
        faces = OCC.meshAllFacesTri(self.hook, edges, True, faceList, hList)
        return faces

def loft(listSketches=[], name="loft"):
    """Create a loft surface from sketches."""
    return Surface(listSketches=listSketches, name=name, type="loft")

def revolve(sketch, center=(0,0,0), axis=(0,0,1), angle=360., name="revolve"):
    """Create a revolution surface from a sketch."""
    return Surface(listSketches=[sketch],
                   data={'center':center, 'axis':axis, 'angle':angle},
                   name=name, type="revolve")

def compound(listSurfaces=[], name="compound"):
    """Create a compound surface from a list of surfaces."""
    return None # a faire

def fill(sketch, continuity=0, name="fill"):
    """Create a surface that fill a sketch."""
    return Surface(listSketches=[sketch],
                   data={'continuity':continuity},
                   name=name, type="fill")

def boolean(listSurfaces=[], name="boolean"):
    return None # a faire

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
                id = DRIVER.scalars[s].id
                out += id
                vars.append(id)
            else: out += s
        return vars, out

#============================================================
class Driver:
    """Driver is Model"""
    def __init__(self):
        # all parameters
        self.scalars = {} # id -> scalar
        self.scalars2 = {} # symbol -> scalar
        self.points = {} # points
        self.grids = {} # grids
        # all entities
        self.edges = {} # edges
        self.sketches = {} # wires
        self.surfaces = {} # shapes
        # all equations
        self.equationCount = 0
        self.equations = {} # equations
        # updated by solve
        self.solution = None # solution of system in sympy symbols
        self.params = None # all model params in sympy symbols
        self.freeParams = None # all model free params in sympy symbols (order)
        # DOE
        self.doeFileName = 'doe.hdf'
        self.doeRange = [] # list param->discretized range
        self.doeSize = [] # list param->size of discretization
        # ROM
        self.romFileName = 'rom.hdf'
        self.K = 0 # reduced dimension
        self.mean = None
        self.Phi = None # POD vectors
        self.ak = None # POD coefficients of each mesh in doe

    def registerScalar(self, s):
        """Register parametric scalar."""
        self.scalars[s.id] = s # id -> scalar
        self.scalars2[s.s] = s # symbol -> scalar

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

    def print(self):
        """Print registered entities."""
        for k in self.scalars: print(k)
        for k in self.points: print(k)
        for k in self.edges: print(k)
        for k in self.sketches: print(k)
        for k in self.surfaces: print(k)
        for k in self.equations: print(k)

    def update(self):
        """Update allfrom parameters."""
        for k in self.edges: self.edges[k].update()
        for k in self.sketches: self.sketches[k].update()
        for k in self.surfaces: self.surfaces[k].update()

    def solve2(self):
        """Solve equations to get free parameters."""

        # get free params
        params = []
        for s in self.scalars:
            mu = self.scalars[s]
            if mu.isFree(): params.append(mu.s)
        print('SOLVE: params=', params)

        # get equations, sub fixed params
        equations = []
        for e in self.equations:
            eq = self.equations[e]
            equations.append(eq.s)
            for s in self.scalars:
                mu = self.scalars[s]
                if not mu.isFree(): eq.s.subs(mu.s, mu.v)
        print('SOLVE: eqs=', equations)

        # solve([eq0,eq1], [x0,x1])
        solution = sympy.solve(equations, params)
        print('SOLVE: sol=', solution)

        # number of free vars
        nparams = len(params)
        neqs = len(equations)
        nd = nparams - neqs
        print("SOLVE: nparams=", nparams)
        print("SOLVE: neqs=", neqs)
        print("SOLVE: free params=", nd)

        # who is free at the end?
        freeParams = params[:]
        for s in solution:
            if solution[s].is_Float:
                print('SOLVE: fixed', s, 'to', solution[s])
                self.scalars2[s].v = solution[s]
                if self.scalars2[s].check(): print('=> valid')
                else: print('=> invalid')
            freeParams.remove(s)
        print('SOLVE: free vars=', freeParams)

        self.solution = solution
        self.params = params
        self.freeParams = freeParams
        return solution, freeParams

    # instantiation of free parameters
    # IN: paramValues: dict of free parameters given values
    def instantiate(self, paramValues):
        """Instantiate all from given paramValues."""
        # set freeParams
        for f in self.freeParams:
            self.scalars2[f].v = paramValues[f.name]
            print('SET: fixed', f, 'to', paramValues[f.name])

        # set other vars with solution
        soli = self.solution.copy()
        for k in soli:
            for f in self.freeParams:
                print("SET: set ", f, " to ", paramValues[f.name])
                soli[k] = soli[k].subs(f, paramValues[f.name])

        for s in soli:
            if soli[s].is_Float:
                print('SET: fixed', s, 'to', soli[s])
                self.scalars2[s].v = soli[s]
                if self.scalars2[s].check(): print('SET: => valid')
                else: print('SET: => invalid')

        # update geometries
        self.update()

    # diff (finite difference) of free parameters on discrete mesh
    # if freevars is None, derivate for all free parameters else derivate for given parameters
    def _diff(self, entity, mesh, freeParams=None, deps=1.e-6):
        """Compute all derivatives dX/dmu on entity."""
        import Converter, KCore

        if len(self.freeParams) == 0:
            print("Warning: no free vars.")
            return None

        if freeParams is None: # no param given
            listVars = self.freeParams
        elif isinstance(freeParams, str): # free param given by name
            listVars = []
            for f in self.freeParams:
                if self.scalars2[f].name == freeParams: listVars.append(f)
        elif isinstance(freeParams, list): # suppose list of names
            listVars = []
            for f in self.freeParams:
                if self.scalars2[f].name in freeParams: listVars.append(f)
        else:
            raise TypeError("diff: incorrect freevars.")

        mesho = Converter.copy(mesh)

        for c, f in enumerate(listVars):
            # free vars value dict
            d = {}
            for q in self.freeParams:
                d[q.name] = self.scalars2[q].v
            d[f.name] += deps

            print("DIFF on: ", f.name)

            # update CAD at param+eps
            self.instantiate(d)
            # project mesh on modified CAD
            if entity.type == "surface":
                OCC._projectOnFaces(entity.hook, mesho)
            else:
                OCC._projectOnEdges(entity.hook, mesho)
            #Converter.convertArrays2File(mesh+mesho, 'diff.plt')
            # get derivatives
            Converter._addVars(mesh, ['dx%d'%c, 'dy%d'%c, 'dz%d'%c])
            for p, m in enumerate(mesh):
                pos1 = KCore.isNamePresent(m, 'dx%d'%c)
                pos2 = KCore.isNamePresent(m, 'dy%d'%c)
                pos3 = KCore.isNamePresent(m, 'dz%d'%c)
                p1x = m[1]
                p2x = mesho[p][1]
                p1x[pos1,:] = (p2x[0,:]-p1x[0,:])/deps
                p1x[pos2,:] = (p2x[1,:]-p1x[1,:])/deps
                p1x[pos3,:] = (p2x[2,:]-p1x[2,:])/deps

        # remet le hook original
        d = {}
        for q in self.freeParams:
            d[q.name] = self.scalars2[q].v
        self.instantiate(d)

        return None

    # set DOE deltas for free parameters
    # it is better to set them in scalar.range
    # IN: dict of deltas for each desired free parameter
    # OUT: arange dict
    def setDOE(self, deltas={}):
        # set default
        self.doeRange = []; self.doeSize = []
        for f in self.freeParams: # give order
            p = self.scalars2[f]
            if len(p.range) == 3: # disc given
                self.doeRange.append(numpy.arange(p.range[0], p.range[1], p.range[2]))
                self.doeSize.append(p.range[2])
            else: # set to 2 points
                self.doeRange.append(numpy.linspace(p.range[0], p.range[1], 2))
                self.doeSize.append(2)

        # set dictionary (optional)
        for k in deltas: # free param names
            for c, f in enumerate(self.freeParams):
                if self.scalars2[f].name == k:
                    p = self.scalars2[f]
                    self.doeRange[c] = numpy.arange(p.range[0], p.range[1], deltas[k])
                    self.doeSize[c] = deltas[k]
        return None

    # walk DOE, append snapshots to file
    def walkDOE(self, entity, hmin, hmax, hausd):
        import itertools
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
                p = self.scalars2[f]
                values[p.name] = val
            if Cmpi.rank == hash%Cmpi.size and hash < raf:
                # instantiate and mesh
                self.instantiate(values)
                mesh = entity.mesh(hmin, hmax, hausd)
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
                mesh = entity.mesh(hmin, hmax, hausd)
                self.addSnapshot(hash, mesh)

    # IN: list of indexes for each param
    # OUT: single hash integer (flatten)
    def getHash(self, indexes):
        hash = 0
        for c, i in enumerate(indexes):
            if c == 0: hash = i
            else: hash += i * self.doeSize[c-1]
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

    # Compute a dmesh from a mesh and a deps on a parameter freevar
    def dmesh(self, entity, mesh, freevar, deps=0.1):
        import Converter, Transform
        self._diff(entity, mesh, freevar, deps)
        Converter._initVars(mesh, '{dx0} = {dx0}*%g'%deps)
        Converter._initVars(mesh, '{dy0} = {dy0}*%g'%deps)
        Converter._initVars(mesh, '{dz0} = {dy0}*%g'%deps)
        mesh2 = Transform.deform(mesh, ['dx0','dy0','dz0'])
        return mesh2

    # remesh input mesh to match nv points using refine
    def remesh(self, mesh, nv):
        import Generator
        nm = mesh[1].shape[1]
        if nm != nv:
            power = (nv-1)*1./(nm-1)
            m = Generator.refine(mesh, power, dir=1)
        else: m = mesh
        return m

    # DOE in file
    def createDOE(self, fileName):
        self.doeFileName = fileName
        if Cmpi.rank > 0: return None
        import Converter.PyTree
        t = Converter.PyTree.newPyTree(['Parameters', 'Snapshots'])
        for c, k in enumerate(self.doeRange):
            node = ["p%05d"%c, k, [], 'parameter_t']
            t[2][1][2].append(node)
        Converter.PyTree.convertPyTree2File(t, self.doeFileName)
        return None

    # add snapshot to file
    def addSnapshot(self, hashcode, msh):
        import Converter.Distributed as Distributed
        import Transform, Converter
        msh = Converter.extractVars(msh, ['x','y','z'])
        msh = Transform.join(msh) # merge mesh
        node = ["%05d"%hashcode, msh[1], [], 'snapshot_t']
        Distributed.writeNodesFromPaths(self.doeFileName, 'CGNSTree/Snapshots', node)
        print("ADD: snapshot %d added."%hashcode, flush=True)
        return None

    # read a snapshot, return an mesh array
    def readSnaphot(self, hashcode):
        import Converter.Distributed as Distributed
        nodes = Distributed.readNodesFromPaths(self.doeFileName, ['CGNSTree/Snapshots/%05d'%hashcode])
        msh = nodes[0][1]
        msh = ['x,y,z', msh, msh.shape[1], 1, 1]
        return msh

    # read all snapshots, return a flatten matrix
    # IN: paramSlab: optional ( (5,120), (3,20), ...) for each parameters
    def readAllSnapshots(self, paramSlab=None):
        import itertools
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

    # ROM
    def writeROM(self, fileName):
        self.romFileName = fileName
        if Cmpi.rank > 0: return None
        import Converter.PyTree
        t = Converter.PyTree.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'phi_t']
        t[2][1][2].append(node)
        Converter.PyTree.convertPyTree2File(t, self.romFileName)
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
        import itertools
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
