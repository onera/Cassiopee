# Parametric CAD driver
import OCC
import sympy
import numpy
import re

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

        # global parameters (always added)
        P = Vec3((0,0,0), name='%s.position'%self.name)
        self.P.append(P)
        P = Point((0,0,0), name='%s.rotCenter'%self.name)
        self.P.append(P)
        P = Vec3((0,0,1), name='%s.rotAxis'%self.name)
        self.P.append(P)
        P = Scalar(0., name='%s.rotAngle'%self.name)
        self.P.append(P)

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

        # global positionning
        OCC._translate(self.hook, self.P[-4].v())
        OCC._rotate(self.hook, self.P[-3].v(), self.P[-2].v(), self.P[-1].v)

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
    """Define a parametric sketch."""
    def __init__(self, listEdges=[], name="sketch"):
        # name
        self.name = getName("sketch")
        # entities
        self.edges = listEdges
        self.hook = None
        self.update()
        # register
        DRIVER.registerSketch(self)

    def add(self, edge):
        self.edges.append(edge)

    # update the CAD from parameters
    def update(self):
        if self.hook is not None: OCC.occ.freeHook(self.hook)
        self.hook = OCC.occ.createEmptyCAD('unknown.step', 'fmt_step')
        hooks = []
        for e in self.edges: hooks.append(e.hook)
        self.hook = OCC.occ.mergeCAD(hooks)

    # check if parameters are valid
    def check(self):
        for e in self.edges:
            ret = e.check()
            if ret == 1: return 1
        return 0

    # print information
    def print(self, shift=0):
        for e in self.edges:
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
class Eq:
    """Equation"""
    def __init__(self, expr1, expr2=None):
        # references sur l'equation sympy
        self.s = sympy.Eq(expr1, expr2)
        DRIVER.registerEquation(self)

    def analyse(self):
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
        # updated when creating scalar, points, entities
        self.scalars = {} # id -> scalar
        self.scalars2 = {} # symbol -> scalar
        self.points = {} # points
        self.grids = {} # grids

        self.edges = {} # edges
        self.sketches = {} # wires
        self.surfaces = {} # shapes
        self.equationCount = 0
        self.equations = {} # equations

        # updated by solve
        self.solution = None # solution of system in sympy symbols
        self.vars = None # all model vars in sympy symbols
        self.freevars = None # all model free vars in sympy symbols

    def registerScalar(self, s):
        self.scalars[s.id] = s # id -> scalar
        self.scalars2[s.s] = s # symbol -> scalar

    def registerPoint(self, p):
        self.points[p.name] = p

    def registerGrid(self, p):
        self.grids[p.name] = p

    def registerEdge(self, e):
        self.edges[e.name] = e

    def registerSketch(self, e):
        self.sketches[e.name] = e

    def registerEquation(self, eq):
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
        # update from parameters
        for k in self.edges: self.edges[k].update()
        for k in self.sketches: self.sketches[k].update()

    def solve2(self):
        # solve all equations and all parameters

        # get free vars
        vars = []
        for s in self.scalars:
            mu = self.scalars[s]
            if mu.isFree(): vars.append(mu.s)
        print('SOLVE: vars=', vars)

        # get equations, sub fixed vars
        equations = []
        for e in self.equations:
            eq = self.equations[e]
            equations.append(eq.s)
            for s in self.scalars:
                mu = self.scalars[s]
                if not mu.isFree(): eq.s.subs(mu.s, mu.v)
        print('SOLVE: eqs=', equations)

        # solve([eq0,eq1], [x0,x1])
        solution = sympy.solve(equations, vars)
        print('SOLVE: sol=', solution)

        # number of free vars
        nvars = len(vars)
        neqs = len(equations)
        nd = nvars - neqs
        print("SOLVE: vars=", nvars)
        print("SOLVE: neqs=", neqs)
        print("SOLVE: free vars=", nd)

        # who is free at the end?
        freevars = vars[:]
        for s in solution:
            if solution[s].is_Float:
                print('SOLVE: fixed', s, 'to', solution[s])
                self.scalars2[s].v = solution[s]
                if self.scalars2[s].check(): print('=> valid')
                else: print('=> invalid')
            freevars.remove(s)
        print('SOLVE: free vars=', freevars)

        self.solution = solution
        self.vars = vars
        self.freevars = freevars
        return solution, freevars

    # instantiation of free vars
    # IN: freevalues: dict
    def instantiate(self, freevalues):

        # set freevars
        for f in self.freevars:
            self.scalars2[f].v = freevalues[f.name]
            print('SET: fixed', f, 'to', freevalues[f.name])

        # set other vars with solution
        soli = self.solution.copy()
        for k in soli:
            for f in self.freevars:
                print("SET: set ", f, " to ", freevalues[f.name])
                soli[k] = soli[k].subs(f, freevalues[f.name])

        #print(soli)
        for s in soli:
            if soli[s].is_Float:
                print('SET: fixed', s, 'to', soli[s])
                self.scalars2[s].v = soli[s]
                if self.scalars2[s].check(): print('SET: => valid')
                else: print('SET: => invalid')

        # update geometry
        self.update()

    # diff (finite difference) of free parameters on discrete mesh
    def _diff(self, entity, mesh, deps=1.e-6):
        import Converter, KCore

        freevars = self.freevars
        if len(freevars) == 0:
            print("Warning: no free vars.")
            return None

        mesho = Converter.copy(mesh)

        for c, f in enumerate(freevars):
            # free vars value dict
            d = {}
            for q in freevars:
                d[q.name] = self.scalars2[q].v
            d[f.name] += deps
            print("DIFF on: ", f.name)

            # update CAD at param+eps
            self.instantiate(d)
            # project mesh on modified CAD
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
                print(pos1, pos2, pos3)
                p1x[pos1,:] = (p2x[0,:]-p1x[0,:])/deps
                p1x[pos2,:] = (p2x[1,:]-p1x[1,:])/deps
                p1x[pos3,:] = (p2x[2,:]-p1x[2,:])/deps

        # remet le hook original
        d = {}
        for q in freevars:
            d[q.name] = self.scalars2[q].v
        self.instantiate(d)

        return None

#============================================================
# Global
DRIVER=Driver()
