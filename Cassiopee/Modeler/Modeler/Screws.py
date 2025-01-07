"""All models of Screws."""
import Geom as D
import Transform as T
import Generator as G
import Converter as C
import Post as P
import Intersector
from .Boxes import box

#==============================================================================
# Hexa head of screw
# IN: X: center position
# IN: r: radius of screw
# IN: h: height of screw
#==============================================================================
def hexaScrew(X, r, h, chamfer=-1):
    if chamfer < 0:
        a = D.circle(X, r, N=7)
        b = T.translate(a, (0,0,h))
        c = G.stack([a,b])
        a = G.tetraMesher(a)
        a = T.reorder(a, (-1,))
        b = G.tetraMesher(b)
        c = C.convertArray2Tetra(c)
        a = T.join([a,c,b])
    else:
        a = D.circle(X, r, N=7)
        b = T.translate(a, (0,0,h-chamfer))
        d = T.translate(a, (0,0,h))
        d = T.scale(d, 1.-chamfer)
        c = G.stack([a,b,d])
        a = G.tetraMesher(a)
        a = T.reorder(a, (-1,))
        d = G.tetraMesher(d)
        c = C.convertArray2Tetra(c)
        c = T.reorder(c, (-1,))
        a = T.join([a,c,d])
    return a

#==============================================================================
# Rounded head of screw
# IN: X: position
# IN: r: radius
# IN: h: height
# IN: drive: type of drive. 0: None, 1: straight
#==============================================================================
def roundScrew(X, r, h, drive=0):
    p = D.polyline([(X[0],X[1],X[2]+h), (X[0]+0.8*r,X[1],X[2]+h*0.8), (X[0]+r,X[1],X[2])])
    p = D.spline(p, N=10)
    p = D.axisym(p, X, axis=(0,0,1), Ntheta=20)
    p = C.convertArray2Tetra(p)
    p = G.close(p)
    p = T.reorder(p, (-1,))
    e = P.exteriorFaces(p)
    e = G.tetraMesher(e)
    e = T.reorder(e, (-1,))
    a = T.join(p, e)
    if drive:
        bb = box((X[0]-2*r, X[1]-0.3*h, X[2]+0.5*h), (X[0]+2*r, X[1]+0.3*h, X[2]+2.5*h))
        a = Intersector.booleanMinus(a, bb)
        return a
    return a