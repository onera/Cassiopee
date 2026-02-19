# - tkDraw -
"""Applet to draw curves."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Geom.PyTree as D
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.Vector as Vector
import time, math
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

CURRENTPOLYLINE = []
CURRENTZONE = None
ALLZONES = []

#==============================================================================
# Recupere les vecteurs unitaires du canvas (si il existe)
#==============================================================================
def getVectorsFromCanvas():
    e1 = (1,0,0); e2 = (0,1,0)
    node = Internal.getNodeFromName1(CTK.t, 'CANVAS')
    if node is None: return e1,e2
    zones = Internal.getZones(node)
    if zones == []: return e1,e2
    zone = zones[0]
    if (Internal.getZoneType(zone) != 1): return e1,e2
    [x1,y1,z1] = C.getValue(zone, Internal.__GridCoordinates__, (1,1,1))
    [x2,y2,z2] = C.getValue(zone, Internal.__GridCoordinates__, (2,1,1))
    [x3,y3,z3] = C.getValue(zone, Internal.__GridCoordinates__, (1,2,1))
    e1 = (x2-x1, y2-y1, z2-z1)
    e2 = (x3-x1, y3-y1, z3-z1)
    e1 = Vector.normalize(e1)
    e2 = Vector.normalize(e2)
    return e1,e2

#==============================================================================
def setSurface():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    VARS[2].set(selected)

#=============================================================================
def getSurfaces():
    name = VARS[2].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in zones:
                if z[0] == sname[1]: surfaces.append(z)
    return surfaces

#==============================================================================
def draw():
    if CTK.t == []: return
    type = VARS[0].get()
    npts = CTK.varsFromWidget(VARS[1].get(), 2)
    if len(npts) != 1:
        CTK.TXT.insert('START', 'Invalid number of points.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    npts = npts[0]
    if type == 'Polyline': drawPolyline()
    elif type == 'Line': drawLine(npts)
    elif type == 'Circle': drawCircle(npts)
    elif type == 'Circular arc': drawArc(npts)
    elif type == 'Rectangle': drawRectangle(npts)
    elif type == 'Cubic': drawCubic(npts)
    elif type == 'Free hand': drawFreeHand()

#==============================================================================
def drawLine(npts):
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    nodes = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(nodes[0], CTK.t)
    CTK.TXT.insert('START', 'Click first point...\n')
    prev = []
    if CTK.__BUSY__ == False:
        CTK.__BUSY__ = True
        TTK.sunkButton(WIDGETS['draw'])
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            CTK.saveTree()
            surfaces = getSurfaces()
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                WIDGETS['draw'].update()
                if CTK.__BUSY__ == False: break
            if CTK.__BUSY__:
                if prev == []:
                    prev = l
                    CTK.TXT.insert('START', 'Click second point...\n')
                elif prev != l:
                    line = D.line(prev, l, npts)
                    if surfaces != []: line = T.projectOrthoSmooth(line, surfaces)
                    CTK.add(CTK.t, nob, -1, line)
                    CTK.TXT.insert('START', 'Line created.\n')
                    CTK.__BUSY__ = False
                    TTK.raiseButton(WIDGETS['draw'])
                    #C._fillMissingVariables(CTK.t)
                    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                    CTK.TKTREE.updateApp()
                    CPlot.render()
                    CPlot.setState(cursor=0)
                    prev = []
                    return
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['draw'])
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['draw'])
        CPlot.setState(cursor=0)
    return

#==============================================================================
def drawCircle(npts):
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    nodes = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(nodes[0], CTK.t)
    CTK.TXT.insert('START', 'Click first point...\n')
    w = WIDGETS['draw']
    prev = []; second = []
    if CTK.__BUSY__ == False:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            CTK.saveTree()
            surfaces = getSurfaces()
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                if prev == []:
                    prev = l
                    CTK.TXT.insert('START', 'Click second point...\n')
                elif second == [] and prev != l:
                    second = l
                    CTK.TXT.insert('START', 'Click third point...\n')
                elif prev != l and second != l:
                    x1 = l[0]; y1 = l[1]; z1 = l[2]
                    x2 = prev[0]; y2 = prev[1]; z2 = prev[2]
                    x3 = second[0]; y3 = second[1]; z3 = second[2]
                    xa = x2 - x1; ya = y2 - y1; za = z2 - z1
                    xb = x3 - x1; yb = y3 - y1; zb = z3 - z1
                    xc = x3 - x2; yc = y3 - y2; zc = z3 - z2
                    a2 = xa*xa + ya*ya + za*za
                    b2 = xb*xb + yb*yb + zb*zb
                    c2 = xc*xc + yc*yc + zc*zc
                    A = 2*b2*c2 + 2*c2*a2 + 2*a2*b2 - a2*a2 - b2*b2 - c2*c2
                    R = math.sqrt( a2*b2*c2 / A )

                    nx = ya*zb - za*yb
                    ny = za*xb - xa*zb
                    nz = xa*yb - ya*xb
                    tx = ya*nz - za*ny
                    ty = za*nx - xa*nz
                    tz = xa*ny - ya*nx
                    norm = tx*tx + ty*ty + tz*tz
                    normi = 1./math.sqrt(norm)
                    tx = tx*normi; ty = ty*normi; tz = tz*normi;
                    alpha = R*R - (xa*xa+ya*ya+za*za)*0.25
                    if alpha >= 0: alpha = math.sqrt(alpha)
                    else: alpha = 0.
                    center = [0,0,0]
                    center[0] = 0.5*(x1+x2) + alpha*tx
                    center[1] = 0.5*(y1+y2) + alpha*ty
                    center[2] = 0.5*(z1+z2) + alpha*tz
                    l = (center[0]-x3)*(center[0]-x3) + \
                        (center[1]-y3)*(center[1]-y3) + \
                        (center[2]-z3)*(center[2]-z3)
                    if (abs(l - R*R) > 1.e-10):
                        center[0] = 0.5*(x1+x2) - alpha*tx
                        center[1] = 0.5*(y1+y2) - alpha*ty
                        center[2] = 0.5*(z1+z2) - alpha*tz
                        l = (center[0]-x3)*(center[0]-x3) + \
                            (center[1]-y3)*(center[1]-y3) + \
                            (center[2]-z3)*(center[2]-z3)
                    circle = D.circle( (center[0],center[1],center[2]), R, N=npts)
                    e1 = [x1-center[0], y1-center[1], z1-center[2]]
                    e2 = [x2-center[0], y2-center[1], z2-center[2]]
                    e3 = Vector.cross(e1, e2)
                    if (e3[0]*e3[0]+e3[1]*e3[1]+e3[2]*e3[2]) < 1e-24:
                        e2 = [x3-center[0], y3-center[1], z3-center[2]]
                        e3 = Vector.cross(e1, e2)
                    e4 = Vector.cross(e1, e3)
                    circle = T.rotate(circle,
                                      (center[0], center[1], center[2]),
                                      ((1,0,0), (0,1,0), (0,0,1)),
                                      (e1, e4, e3))
                    if surfaces != []:
                        circle = T.projectOrthoSmooth(circle, surfaces)
                    CTK.add(CTK.t, nob, -1, circle)
                    CTK.TXT.insert('START', 'Circle created.\n')
                    CTK.__BUSY__ = False
                    TTK.raiseButton(w)
                    CPlot.setState(cursor=0)
                    #C._fillMissingVariables(CTK.t)
                    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                    CTK.TKTREE.updateApp()
                    CPlot.render()
                    CPlot.setState(cursor=0)
                    prev = []
                    return
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
def drawArc(npts):
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    nodes = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(nodes[0], CTK.t)
    CTK.TXT.insert('START', 'Click first point...\n')
    w = WIDGETS['draw']
    prev = []; second = []
    if not CTK.__BUSY__:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            CTK.saveTree()
            surfaces = getSurfaces()
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                if prev == []:
                    prev = l
                    CTK.TXT.insert('START', 'Click second point...\n')
                elif second == [] and prev != l:
                    second = l
                    CTK.TXT.insert('START', 'Click third point...\n')
                elif prev != l and second != l:
                    x1 = l[0]; y1 = l[1]; z1 = l[2]
                    x2 = prev[0]; y2 = prev[1]; z2 = prev[2]
                    x3 = second[0]; y3 = second[1]; z3 = second[2]
                    xa = x2 - x1; ya = y2 - y1; za = z2 - z1
                    xb = x3 - x1; yb = y3 - y1; zb = z3 - z1
                    xc = x3 - x2; yc = y3 - y2; zc = z3 - z2
                    a2 = xa*xa + ya*ya + za*za
                    b2 = xb*xb + yb*yb + zb*zb
                    c2 = xc*xc + yc*yc + zc*zc
                    A = 2*b2*c2 + 2*c2*a2 + 2*a2*b2 - a2*a2 - b2*b2 - c2*c2
                    if A > 1.e-48: R = math.sqrt(a2*b2*c2 / A)
                    else: R = 0.

                    nx = ya*zb - za*yb
                    ny = za*xb - xa*zb
                    nz = xa*yb - ya*xb
                    tx = ya*nz - za*ny
                    ty = za*nx - xa*nz
                    tz = xa*ny - ya*nx
                    norm = tx*tx + ty*ty + tz*tz
                    if norm > 1.e-24: normi = 1./math.sqrt(norm)
                    else: normi = 1.e24

                    tx = tx*normi; ty = ty*normi; tz = tz*normi;
                    alpha = R*R - (xa*xa+ya*ya+za*za)*0.25
                    if alpha >= 0: alpha = math.sqrt(alpha)
                    else: alpha = 0.
                    center = [0,0,0]
                    center[0] = 0.5*(x1+x2) + alpha*tx
                    center[1] = 0.5*(y1+y2) + alpha*ty
                    center[2] = 0.5*(z1+z2) + alpha*tz
                    dx3 = center[0]-x3; dy3 = center[1]-y3; dz3 = center[2]-z3
                    l = dx3*dx3 + dy3*dy3 + dz3*dz3
                    if abs(l - R*R) > 1.e-10:
                        center[0] = 0.5*(x1+x2) - alpha*tx
                        center[1] = 0.5*(y1+y2) - alpha*ty
                        center[2] = 0.5*(z1+z2) - alpha*tz
                        dx3 = center[0]-x3; dy3 = center[1]-y3; dz3 = center[2]-z3
                        l = dx3*dx3 + dy3*dy3 + dz3*dz3

                    e1 = [x1-center[0], y1-center[1], z1-center[2]]
                    e2 = [x2-center[0], y2-center[1], z2-center[2]]
                    e3 = Vector.cross(e1, e2)
                    if (e3[0]*e3[0]+e3[1]*e3[1]+e3[2]*e3[2]) < 1e-24:
                        e2 = [x3-center[0], y3-center[1], z3-center[2]]
                        e3 = Vector.cross(e1, e2)
                    e4 = Vector.cross(e1, e3)

                    # Images des pts dans le plan xyz
                    pt1 = D.point((x1,y1,z1))
                    pt2 = D.point((x2,y2,z2))
                    pt3 = D.point((x3,y3,z3))
                    pt1 = T.rotate(pt1,
                                   (center[0], center[1], center[2]),
                                   (e1, e4, e3),
                                   ((1,0,0), (0,1,0), (0,0,1)) )
                    pt2 = T.rotate(pt2,
                                   (center[0], center[1], center[2]),
                                   (e1, e4, e3),
                                   ((1,0,0), (0,1,0), (0,0,1)))
                    pt3 = T.rotate(pt3,
                                   (center[0], center[1], center[2]),
                                   (e1, e4, e3),
                                   ((1,0,0), (0,1,0), (0,0,1)))
                    xp1 = C.getValue(pt1, 'CoordinateX', 0)
                    yp1 = C.getValue(pt1, 'CoordinateY', 0)
                    zp1 = C.getValue(pt1, 'CoordinateZ', 0)
                    xp2 = C.getValue(pt2, 'CoordinateX', 0)
                    yp2 = C.getValue(pt2, 'CoordinateY', 0)
                    zp2 = C.getValue(pt2, 'CoordinateZ', 0)
                    xp3 = C.getValue(pt3, 'CoordinateX', 0)
                    yp3 = C.getValue(pt3, 'CoordinateY', 0)
                    zp3 = C.getValue(pt3, 'CoordinateZ', 0)

                    dx1 = (xp1-center[0])/R; dy1 = (yp1-center[1])/R
                    if dx1 > 1.: dx1 = 1.
                    if dx1 < -1.: dx1 = -1.
                    if dy1 > 0: teta1 = math.acos(dx1)
                    else: teta1 = 2*math.pi - math.acos(dx1)
                    teta1 = teta1*180./math.pi; teta1 = 360.

                    dx2 = (xp2-center[0])/R; dy2 = (yp2-center[1])/R
                    if dx2 > 1.: dx2 = 1.
                    if dx2 < -1.: dx2 = -1.
                    if dy2 > 0: teta2 = math.acos(dx2)
                    else: teta2 = 2*math.pi - math.acos(dx2)
                    teta2 = teta2*180./math.pi

                    dx3 = (xp3-center[0])/R; dy3 = (yp3-center[1])/R
                    if dx3 > 1.: dx3 = 1.
                    if dx3 < -1.: dx3 = -1.
                    if dy3 > 0: teta3 = math.acos(dx3)
                    else: teta3 = 2*math.pi - math.acos(dx3)
                    teta3 = teta3*180./math.pi

                    if teta3 > teta2: teta1 = 360.
                    else: teta1 = 0.

                    circle = D.circle((center[0],center[1],center[2]), R,
                                      tetas=teta2, tetae=teta1, N=npts)
                    circle = T.rotate(circle,
                                      (center[0], center[1], center[2]),
                                      ((1,0,0), (0,1,0), (0,0,1)),
                                      (e1, e4, e3))
                    if surfaces != []:
                        circle = T.projectOrthoSmooth(circle, surfaces)
                    CTK.add(CTK.t, nob, -1, circle)
                    CTK.TXT.insert('START', 'Circle created.\n')
                    CTK.__BUSY__ = False
                    TTK.raiseButton(w)
                    CPlot.setState(cursor=0)
                    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                    CTK.TKTREE.updateApp()
                    CPlot.render()
                    CPlot.setState(cursor=0)
                    prev = []
                    return
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
def drawRectangle(npts):
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    nodes = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(nodes[0], CTK.t)
    CTK.TXT.insert('START', 'Click left/lower corner...\n')
    w = WIDGETS['draw']
    prev = []; second = []
    if not CTK.__BUSY__:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            CTK.saveTree()
            surfaces = getSurfaces()
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                if prev == []:
                    prev = l
                    CTK.TXT.insert('START', 'Click right/up corner...\n')
                elif prev != l:
                    e1,e2 = getVectorsFromCanvas()
                    e1n = Vector.norm(e1)
                    e2n = Vector.norm(e2)
                    if e2n > e1n: e1 = e2
                    P1 = l; P2 = prev
                    P1P2 = Vector.sub(P2, P1)
                    P1P2n = Vector.norm(P1P2)
                    Q = Vector.norm(Vector.cross(e1, P1P2))
                    L = math.sqrt( P1P2n*P1P2n - Q*Q )
                    sign = Vector.dot(e1, P1P2)
                    if sign > 0: e1 = Vector.mul(L, e1)
                    else: e1 = Vector.mul(-L, e1)
                    P3 = Vector.add(P1, e1)
                    P4 = Vector.sub(P2, e1)
                    l1 = D.line(P1, P3, npts)
                    l2 = D.line(P3, P2, npts)
                    l3 = D.line(P2, P4, npts)
                    l4 = D.line(P4, P1, npts)
                    rect = T.join([l1,l2,l3,l4])
                    if surfaces != []:
                        rect = T.projectOrthoSmooth(rect, surfaces)
                    CTK.add(CTK.t, nob, -1, rect)
                    CTK.TXT.insert('START', 'Rectangle created.\n')
                    CTK.__BUSY__ = False
                    TTK.raiseButton(w)
                    CPlot.setState(cursor=0)
                    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                    CTK.TKTREE.updateApp()
                    CPlot.render()
                    CPlot.setState(cursor=0)
                    prev = []
                    return
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
def drawPolyline():
    global CURRENTZONE; global CURRENTPOLYLINE
    if CTK.t == []: return
    w = WIDGETS['draw']

    if not CTK.__BUSY__:
        CPlot.unselectAllZones()
        CTK.saveTree()
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                CPlot.unselectAllZones()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                CURRENTPOLYLINE.append((l[0],l[1],l[2]))
                if CURRENTZONE is None:
                    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
                    base = Internal.getNodeFromName1(CTK.t, 'CONTOURS')
                    nob = C.getNobOfBase(base, CTK.t)
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.add(CTK.t, nob, -1, a)
                    ret = Internal.getParentOfNode(CTK.t, CURRENTZONE)
                    noz = ret[1]
                else:
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.replace(CTK.t, nob, noz, a)
                (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                CTK.TKTREE.updateApp()
                CPlot.render()
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        surfaces = getSurfaces()
        if surfaces != []:
            ret = Internal.getParentOfNode(CTK.t, CURRENTZONE)
            nob = C.getNobOfBase(ret[0], CTK.t)
            a = T.projectOrthoSmooth(CURRENTZONE, surfaces)
            noz = ret[1]
            CTK.replace(CTK.t, nob, noz, a)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        CURRENTZONE = None
        CURRENTPOLYLINE = []
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
def drawCubic(npts):
    global CURRENTZONE; global CURRENTPOLYLINE
    if CTK.t == []: return
    w = WIDGETS['draw']
    if not CTK.__BUSY__:
        CPlot.unselectAllZones()
        CTK.saveTree()
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                CPlot.unselectAllZones()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                CURRENTPOLYLINE.append((l[0],l[1],l[2]))
                if CURRENTZONE is None:
                    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
                    base = Internal.getNodeFromName1(CTK.t, 'CONTOURS')
                    nob = C.getNobOfBase(base, CTK.t)
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.add(CTK.t, nob, -1, a)
                    ret = Internal.getParentOfNode(CTK.t, CURRENTZONE)
                    noz = ret[1]
                else:
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.replace(CTK.t, nob, noz, a)
                (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                CTK.TKTREE.updateApp()
                CPlot.render()
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        ret = Internal.getParentOfNode(CTK.t, CURRENTZONE)
        a = D.polyline(CURRENTPOLYLINE)
        d = G.cart( (0,0,0), (1./(npts-1),1,1), (npts,1,1) )
        a = G.map(a, d)
        surfaces = getSurfaces()
        if surfaces != []: a = T.projectOrthoSmooth(a, surfaces)
        nob = C.getNobOfBase(ret[0], CTK.t)
        CTK.replace(CTK.t, nob, ret[1], a)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        CURRENTZONE = None
        CURRENTPOLYLINE = []
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
def drawFreeHand():
    global CURRENTZONE; global CURRENTPOLYLINE; global ALLZONES
    w = WIDGETS['draw']
    prev = []; first = []
    if CTK.__BUSY__ == False:
        CPlot.unselectAllZones()
        CTK.saveTree()
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        buttonState = 0
        while CTK.__BUSY__:
            if prev == []: # first point
                l = []
                while l == []:
                    l = CPlot.getActivePoint()
                    if l != []: prev = l; first = l
                    time.sleep(CPlot.__timeStep__)
                    w.update()
                    if not CTK.__BUSY__: break
            else: # next points
                diff = -1.
                while (diff < 1.e-10):
                    (buttonState,x,y,z) = CPlot.getMouseState()
                    l = (x,y,z)
                    diff = Vector.norm2(Vector.sub(l,prev))
                    diff1 = Vector.norm2(Vector.sub(l,first))
                    if diff1 < 1.e-10: l = first
                    if buttonState == 5: break
                    time.sleep(CPlot.__timeStep__)
                    w.update()
                    if not CTK.__BUSY__: break

            prev = l
            CPlot.unselectAllZones()
            if buttonState == 5: # button released
                ALLZONES.append(CURRENTZONE)
                CURRENTZONE = None; prev = []; first = []
                CURRENTPOLYLINE = []
                CTK.TKTREE.updateApp()

            if (CTK.__BUSY__ == True and buttonState != 5):
                CURRENTPOLYLINE.append((l[0],l[1],l[2]))
                if CURRENTZONE is None:
                    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
                    base = Internal.getNodeFromName1(CTK.t, 'CONTOURS')
                    nob = C.getNobOfBase(base, CTK.t)
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.add(CTK.t, nob, -1, a)
                    ret = Internal.getParentOfNode(CTK.t, CURRENTZONE)
                    noz = ret[1]
                else:
                    a = D.polyline(CURRENTPOLYLINE)
                    CURRENTZONE = a
                    CTK.replace(CTK.t, nob, noz, a)
                (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
                CPlot.render()
            buttonState = 0
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        surfaces = getSurfaces()
        if surfaces != []:
            if CURRENTZONE is not None: ALLZONES += [CURRENTZONE]
            for s in ALLZONES:
                ret = Internal.getParentOfNode(CTK.t, s)
                nob = C.getNobOfBase(ret[0], CTK.t)
                a = T.projectOrthoSmooth(s, surfaces)
                noz = ret[1]
                CTK.replace(CTK.t, nob, noz, a)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        CURRENTZONE = None; ALLZONES = []
        CURRENTPOLYLINE = []
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDraw  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Drawing tools.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=1)

    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkDraw')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Figure type -
    V = TK.StringVar(win); V.set('Line'); VARS.append(V)
    if 'tkDrawType' in CTK.PREFS: V.set(CTK.PREFS['tkDrawType'])
    # -1- Npts -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkDrawNpts' in CTK.PREFS: V.set(CTK.PREFS['tkDrawNpts'])
    # -2- underlaying surface
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Surface -
    B = TTK.Button(Frame, text="Surf", command=setSurface,
                   image=iconics.PHOTO[8], compound=TK.RIGHT, padx=0, pady=0)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set underlaying model surfaces (you draw onto).')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Underlaying model surfaces.')

    # - Figure type -
    B = TTK.OptionMenu(Frame, VARS[0], 'Line', 'Circle', 'Circular arc',
                       'Rectangle', 'Polyline', 'Cubic', 'Free hand')
    BB = CTK.infoBulle(parent=B, text='Type of drawing.')
    B.grid(row=1, column=0, sticky=TK.EW)

    # - Npts -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of points.')

    # - Draw -
    B = TTK.Button(Frame, text="Draw", command=draw)
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Click to start to draw. Click again to end.')
    WIDGETS['draw'] = B

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['EdgeNoteBook'].add(WIDGETS['frame'], text='tkDraw')
    except: pass
    CTK.WIDGETS['EdgeNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['EdgeNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkDrawType'] = VARS[0].get()
    CTK.PREFS['tkDrawNpts'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('Line')
    VARS[1].set('10')
    CTK.PREFS['tkDrawType'] = VARS[0].get()
    CTK.PREFS['tkDrawNpts'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkDraw '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
