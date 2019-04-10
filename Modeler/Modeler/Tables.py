# - all tables -

from . import Boxes
import Geom as D
import Transform as T
import Post as P
import Generator as G
import Converter as C

#==============================================================================
# Table rectangulaire simple
# L: length: 
# W: width
# H: height
#==============================================================================
def table1(L=2., W=1., H=1.):
    """Simple rectangular table.
    Usage: table1(L, W, H)"""
    b = Boxes.box((0,0,0),(L,W,H*0.1))
    eps = min(L,W)*0.1
    eps2 = eps*0.3
    p1 = Boxes.box((eps2,eps2,-H),(eps2+eps,eps2+eps,0))
    p2 = Boxes.box((eps2,W-eps2,-H),(eps2+eps,W-eps2-eps,0))
    p3 = Boxes.box((L-eps2,eps2,-H),(L-eps2-eps,eps2+eps,0))
    p4 = Boxes.box((L-eps2,W-eps2,-H),(L-eps2-eps,W-eps2-eps,0))
    return [b,p1,p2,p3,p4]

#==============================================================================
# Table ronde de cafe
# L: rayon de la table
# H: height
#==============================================================================
def table2(L=0.6, H=1., N=20):
    """Simple round table.
    Usage: table2(L, H, N)"""
    e1 = 0.1*H # epaisseur plateau
    e2 = 0.05*H # rayon pied
    p1 = D.polyline([(L,0,0), (L,e1,0), (e2,e1,0), (e2,0.8*H,0)])
    p2 = D.polyline([(0.6*L,H,0), (e2,0.8*H,0)])
    l = [p1,p2]
    #l = D.connect1D([p1,p2], sharpness=1, N=5, lengthFactor=0.2)
    a = D.axisym(l, (0,0,0), (0,1,0), 360., N)
    a = T.rotate(a, (0,0,0), (1,0,0), -90.)
    a = C.convertArray2Hexa(a)
    a = T.join(a)
    a = G.close(a)
    ex = P.exteriorFaces(a)
    ex = T.splitConnexity(ex)
    p = G.fittingPlaster(ex[0], bumpFactor=0.)
    b = G.gapfixer(ex[0], p)
    a = C.convertArray2Tetra(a)
    a = T.join([a,b])
    return [a]
