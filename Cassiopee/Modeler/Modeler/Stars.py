# - all stars -
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Etoile plane
# IN: R1: rayon interieur
# IN: R2: rayon exterieur
# IN: shift: decalage des cercles
# IN: h: epaisseur
#==============================================================================
def star(R1=1., R2=2., shift=0.5, N=10, h=1.):
    c1 = D.circle((0,0,0), R2, N=N)
    c2 = D.circle((0,0,0), R1, N=N)
    c2 = T.rotate(c2, (0,0,0), (0,0,1), shift*360./(N-1))
    points = []
    xc1 = c1[1][0]; yc1 = c1[1][1]
    xc2 = c2[1][0]; yc2 = c2[1][1]
    for i in range(N):
        points.append((xc1[i],yc1[i],0))
        points.append((xc2[i],yc2[i],0))
    c = D.polyline(points)
    p = G.fittingPlaster(c)
    b = G.gapfixer(c, p)
    b2 = T.translate(b, (0,0,h))
    l = D.line((0,0,0), (0,0,h), N=3)
    c = D.lineDrive(c, l)
    c = C.convertArray2Tetra(c)
    o = T.join([b,c,b2])
    o = G.close(o)
    o = T.reorder(o, (-1,))
    return o
