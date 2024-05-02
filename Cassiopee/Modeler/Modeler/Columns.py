# - All columns -
import Geom as D
import Transform as T
import Generator as G
import Converter as C
from . import Boxes
import math

#==============================================================================
# Column avec tete et pied carre
# IN: R: rayon
# IN: N: nbre de pts
# IN: h: hauteur de la colonne
#==============================================================================
def column(R=0.2, N=10, h=1.):
    c = D.circle((0,0,0), R, N=N)
    p = G.fittingPlaster(c)
    b = G.gapfixer(c, p)
    b2 = T.translate(b, (0,0,h))
    l = D.line((0,0,0), (0,0,h), N=5)
    c2 = D.lineDrive(c, l)
    c2 = C.convertArray2Tetra(c2)
    o = T.join([b,b2,c2])
    o = G.close(o)
    o = T.reorder(o, (-1,))
    dh = R*math.sqrt(2.)
    box1 = Boxes.box((-dh,-dh,-R),(dh,dh,0), chamfer=0.05*R)
    box2 = Boxes.box((-dh,-dh,h), (dh,dh,h+R), chamfer=0.05*R)
    box1 = C.convertArray2Tetra(box1)
    box2 = C.convertArray2Tetra(box2)
    o = T.join([o,box1,box2])
    return o

#==============================================================================
# Colonne avec tete et pied ronds
# IN: R1: rayon colonne au pied
# IN: R2: rayon colonne a la tete
# IN: N: nbre de pts
# IN: h: hauteur de la colonne
#==============================================================================
def column2(R1=0.2, R2=0.2, N=10, h=1.):
    l = D.line((R1,0,0),(R2,0,h), N=N//2+2)
    Rc = R2; Rc2 = Rc*0.8
    p = D.polyline([(R1,0,0), (R1+Rc,0,0), (R1+Rc,0,-Rc2), (R1,0,-Rc2)])
    s1 = D.spline(p, N=N, M=N)
    p = D.polyline([(R2,0,h), (R2+Rc,0,h), (R2+Rc,0,h+Rc2), (R2,0,Rc2+h)])
    s2 = D.spline(p, N=N, M=N)
    c = C.convertArray2Tetra([s1,l,s2])
    c = T.join(c)
    c = G.close(c)
    o = D.axisym(c, (0,0,0), (0,0,1), 360, N)
    return o

#==============================================================================
def column3(R1=0.2, R2=0.2, N=10, h=1.):
    from . import Circles
    a = Circles.circle1(1., 0.9, Nd=N//2, fracD=0.2, N=N)
    l = D.line((R1,0,0),(R2,0,h), N=N//2+2)
    o = D.axisym(l, (0,0,0), (0,0,1), rmod=a)
    Rc = R2; Rc2 = Rc*0.8
    p = D.polyline([(R1,0,0), (R1+Rc,0,0), (R1+Rc,0,-Rc2), (R1,0,-Rc2)])
    s1 = D.spline(p, N=N, M=N)
    p = D.polyline([(R2,0,h), (R2+Rc,0,h), (R2+Rc,0,h+Rc2), (R2,0,Rc2+h)])
    s2 = D.spline(p, N=N, M=N)
    o1 = D.axisym(s1, (0,0,0), (0,0,1), 360, N)
    o2 = D.axisym(s2, (0,0,0), (0,0,1), 360, N)
    g = [o,o1,o2]
    g = C.convertArray2Hexa(g)
    g = T.join(g)
    g = G.close(g)
    return g
