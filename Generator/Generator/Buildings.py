# Building mesher

import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Transform.PyTree as T

# Maillage en O autour de la tour
# IN: base: liste de 4 BARs decrivant la base du building
# IN: height: la hauteur
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
    T._reorder(b1, (2,1,3))
    b2 = G.TFI([a01,a1z1,h01,a0z0])
    #T._reorder(b2, (2,1,3))
    b3 = G.TFI([a12,a2z2,h12,a1z1])
    #T._reorder(b3, (2,1,3))
    b4 = G.TFI([a23,a3z3,h23,a2z2])
    T._reorder(b4, (2,1,3))
    b5 = G.TFI([a30,a0z0,h30,a3z3])
    T._reorder(b5, (2,1,3))

    s = [b1,b2,b3,b4,b5]

    d = G.cart((0,0,0), (h,1,1), (9,1,1)) # hauteur = 9*h
    d = G.enforcePlusX(d, hp, 9, 40)
    s = G.addNormalLayers(s, d, niter=50)

    C._addBC2Zone(s, 'wall', 'BCWall', 'kmin')
    C._addBC2Zone(s, 'ov', 'BCOverlap', 'kmax')
    tp = C.newPyTree(['TOWER',s])
    tp = X.connectMatch(tp)
    C._fillEmptyBCWith(tp, 'wall', 'BCWall')
    return Internal.getZones(tp)

# Create domain - a mettre en cartRx
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
    d = G.enforcePlusX(d, hp, npts//2, 40)

    vol = G.map(vol, d, dir=3)

    # BCs
    C._addBC2Zone(vol, 'wall', 'BCWall', 'kmin')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'kmax')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'jmin')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'jmax')
    C._addBC2Zone(vol, 'far', 'BCFarfield', 'imax')
    C._addBC2Zone(vol, 'in', 'BCInflow', 'imin')    

    return [vol]

# Return a square from two points
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
    
def map(P0, P1, fileName):
    s = square(P0,P1)
    #CPlot._addRender2Zone(a, material='texmat')
