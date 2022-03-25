# - all Hubs -
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Hub sans mat
# IN: d: diametre
# IN: h: hauteur
# IN: d2: diametre a mi-hauteur
#==============================================================================
def hub1(d, h, d2):
    P0 = (0,0,0) # top
    P1 = (d,-h,0) # bottom
    P2 = (d2,-h/2.,0) # mid
    print([P0,P2,P1], flush=True)
    c = D.polyline([P0,P2,P1])
    s = D.spline(c, order=3, N=30)
    #s = D.axisym(s, P0, Ntheta=30)
    return s

#==============================================================================
# Hub avec mat
# IN: d: diametre
#==============================================================================
