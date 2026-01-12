# - All leaves -

import Geom as D
import Generator as G
import Transform as T
import Converter as C

#==============================================================================
# leave
# IN: w: width of leave
# IN: h: height
# IN: N: number of points of discretisation
#==============================================================================
def leave(w=1., h=2., N=10):
    l = D.line((0,0,0), (0,h,0),N=N)
    p = D.polyline([(0,0,0), (0.9*w,h*0.1,0), (w,h*0.5,0), (0.,h*0.9,0), (0,h,0)])
    s = D.spline(p, order=3, N=N, M=N)
    c = C.convertArray2Tetra([l,s])
    c = T.join(c)
    pl = G.fittingPlaster(c)
    b = G.gapfixer(c, pl)
    b2 = T.symmetrize(b, (0,0,0), (0,1,0), (0,0,1))
    o = T.join([b,b2])
    o = G.close(o)
    return o
