# - Cylinders/Cones -
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Cylindre
# R1: rayon en bas
# R2: rayon en haut
# h: hauteur
# Rc: rayon de fermeture, si < 0, fermeture plane
#==============================================================================
def cylinder(R1=1., R2=1., N=10, h=1., Rc=-1):
    line = D.line((R1,0,0), (R2,0,h))
    o = D.axisym(line, (0,0,0), (0,0,1), 360, N)
    o = C.convertArray2Tetra(o)

    if Rc < 0: # fermeture plane
        c = D.circle((0,0,0), R1, N=N)
        p = G.fittingPlaster(c)
        b1 = G.gapfixer(c, p)
        c = D.circle((0,0,h), R2, N=N)
        p = G.fittingPlaster(c)
        b2 = G.gapfixer(c, p)
    else: # fermeture ronde
        p = D.polyline([(R1,0,0),(R1,0,-Rc),(0,0,-Rc)])
        s = D.spline(p, order=3, N=10, M=10)
        b1 = D.axisym(s, (0,0,0), (0,0,1), 360, N)
        b1 = C.convertArray2Tetra(b1)
        p = D.polyline([(R2,0,h),(R1,0,h+Rc),(0,0,h+Rc)])
        s = D.spline(p, order=3, N=10, M=10)
        b2 = D.axisym(s, (0,0,0), (0,0,1), 360, N)
        b2 = C.convertArray2Tetra(b2)
    o = T.join([o,b1,b2])
    o = G.close(o)
    o = T.reorder(o, (+1,))
    return o
