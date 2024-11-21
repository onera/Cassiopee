# Building mesher

import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Transform.PyTree as T

#=====================================================================
# Maillage en O autour d'un batiment
# IN: base: liste de 4 BARs decrivant la base du building
# IN: height: la hauteur du batiment
# IN: h: la taille des mailles sur le batiment
# IN: hp: la taille de maille paroi (y+=1)
#=====================================================================
def building(base, height, h, hp):
    if len(base) != 4: raise ValueError('building: the base must have 4 sides.')

    # line 0 with h
    a01 = D.uniformize(base[0], h=h)
    # line 1 with h
    a12 = D.uniformize(base[1], h=h)
    # line 2 with a0 npts
    a23 = D.uniformize(base[2], N=C.getNPts(a01))
    # line 3 with a1 npts
    a30 = D.uniformize(base[3], N=C.getNPts(a12))
    b0 = G.TFI([a01,a12,a23,a30])

    # vertical line
    P0 = C.getValue(a01, 'GridCoordinates', ind=0)
    Pz0 = (P0[0],P0[1],P0[2]+height)
    a0z0 = D.line(P0,Pz0,N=2)
    a0z0 = D.uniformize(a0z0, h=h)
    npts = C.getNPts(a0z0)
    P1 = C.getValue(a12, 'GridCoordinates', ind=0)
    Pz1 = (P1[0],P1[1],P1[2]+height)
    a1z1 = D.line(P1,Pz1,N=npts)
    P2 = C.getValue(a23, 'GridCoordinates', ind=0)
    Pz2 = (P2[0],P2[1],P2[2]+height)
    a2z2 = D.line(P2,Pz2,N=npts)
    P3 = C.getValue(a30, 'GridCoordinates', ind=0)
    Pz3 = (P3[0],P3[1],P3[2]+height)
    a3z3 = D.line(P3,Pz3,N=npts)

    h01 = D.line(Pz0,Pz1,N=C.getNPts(a01))
    h12 = D.line(Pz1,Pz2,N=C.getNPts(a12))
    h23 = D.line(Pz2,Pz3,N=C.getNPts(a23))
    h30 = D.line(Pz3,Pz0,N=C.getNPts(a30))

    b1 = G.TFI([h01,h12,h23,h30])
    b1[0] = C.getZoneName('top')
    T._reorder(b1, (2,1,3))
    b2 = G.TFI([a01,a1z1,h01,a0z0])
    b2[0] = C.getZoneName('front')
    #T._reorder(b2, (2,1,3))
    b3 = G.TFI([a12,a2z2,h12,a1z1])
    b3[0] = C.getZoneName('right')
    #T._reorder(b3, (2,1,3))
    b4 = G.TFI([a23,a3z3,h23,a2z2])
    b4[0] = C.getZoneName('back')
    T._reorder(b4, (2,1,3))
    b5 = G.TFI([a30,a0z0,h30,a3z3])
    b5[0] = C.getZoneName('left')
    T._reorder(b5, (2,1,3))

    s = [b1,b2,b3,b4,b5]
    d = G.cart((0,0,0), (h,1,1), (9,1,1)) # hauteur = 9*h
    d = G.enforcePlusX(d, hp, 9, 40)
    [b1,b2,b3,b4,b5] = G.addNormalLayers(s, d, niter=50)

    # remap near wall
    li = T.subzone(b2, (1,1,1), (-1,1,1))
    L = D.getLength(li)
    D._getCurvilinearAbscissa(li)
    C._initVars(li, '{CoordinateX}={s}')
    li = G.enforcePlusX(li, hp/L, 9, 20)
    b2 = G.map(b2, li, 1)
    b3 = G.map(b3, li, 1)
    b4 = G.map(b4, li, 2)
    b5 = G.map(b5, li, 2)
    s = [b1,b2,b3,b4,b5]

    # Boundary conditions
    C._addBC2Zone(s, 'wall', 'BCWall', 'kmin')
    C._addBC2Zone(s, 'overlap', 'BCOverlap', 'kmax')
    tp = C.newPyTree(['TOWER', s])
    tp = X.connectMatch(tp)
    C._fillEmptyBCWith(tp, 'walls', 'BCWall')

    return Internal.getZones(tp)

#=====================================================================
# Create computational domain
# OUT: single block
#=====================================================================
def domain(base, height, h, hp):
    if len(base) != 4: raise ValueError('building: the base must have 4 sides.')
    # line 0 with h
    a01 = D.uniformize(base[0], h=h)
    # line 1 with h
    a12 = D.uniformize(base[1], h=h)
    # line 2 with a0 npts
    a23 = D.uniformize(base[2], N=C.getNPts(a01))
    # line 3 with a1 npts
    a30 = D.uniformize(base[3], N=C.getNPts(a12))
    b0 = G.TFI([a01,a12,a23,a30])

    # vertical line
    P0 = C.getValue(a01, 'GridCoordinates', ind=0)
    Pz0 = (P0[0],P0[1],P0[2]+height)
    a0z0 = D.line(P0,Pz0,N=2)
    a0z0 = D.uniformize(a0z0, h=h)
    npts = C.getNPts(a0z0)
    P1 = C.getValue(a12, 'GridCoordinates', ind=0)
    Pz1 = (P1[0],P1[1],P1[2]+height)
    a1z1 = D.line(P1,Pz1,N=npts)
    P2 = C.getValue(a23, 'GridCoordinates', ind=0)
    Pz2 = (P2[0],P2[1],P2[2]+height)
    a2z2 = D.line(P2,Pz2,N=npts)
    P3 = C.getValue(a30, 'GridCoordinates', ind=0)
    Pz3 = (P3[0],P3[1],P3[2]+height)
    a3z3 = D.line(P3,Pz3,N=npts)

    h01 = D.line(Pz0,Pz1,N=C.getNPts(a01))
    h12 = D.line(Pz1,Pz2,N=C.getNPts(a12))
    h23 = D.line(Pz2,Pz3,N=C.getNPts(a23))
    h30 = D.line(Pz3,Pz0,N=C.getNPts(a30))

    b1 = G.TFI([h01,h12,h23,h30])

    b2 = G.TFI([a01,a1z1,h01,a0z0])
    b3 = G.TFI([a12,a2z2,h12,a1z1])
    b4 = G.TFI([a23,a3z3,h23,a2z2])
    b5 = G.TFI([a30,a0z0,h30,a3z3])

    # vol zone
    vol = G.TFI([b0,b1,b2,b3,b4,b5])
    vol = T.reorder(vol, (3,2,1))

    # remap en z
    npts = C.getNPts(a0z0)
    d = G.cart((0,0,0), (h,1,1), (npts,1,1))
    D._getCurvilinearAbscissa(d)
    C._initVars(d, '{CoordinateX} = {s}')
    d = G.enforcePlusX(d, hp, min(npts//2,40), 40)

    vol = G.map(vol, d, dir=3)

    # BCs
    C._addBC2Zone(vol, 'wall', 'BCWall', 'kmin')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'kmax')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'jmin')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'jmax')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'imax')
    C._addBC2Zone(vol, 'in', 'BCInflow', 'imin')    

    return [vol]

#=====================================================================
# Return a square from two points
# OUT: list of lines
#=====================================================================
def square(P0, P1):
    (x0,y0,z0) = P0
    (x1,y1,z1) = P1
    P2 = (x1,y0,z0)
    P3 = (x0,y1,z1)
    l0 = D.line(P0,P2,N=2)
    l1 = D.line(P2,P1,N=2)
    l2 = D.line(P1,P3,N=2)
    l3 = D.line(P3,P0,N=2)
    return [l0,l1,l2,l3]

#=====================================================================
def map(P0, P1, fileName):
    s = square(P0,P1)
    #CPlot._addRender2Zone(a, material='texmat')
