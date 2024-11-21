# - all glasses -
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Rp, hp: rayon et hauteur pied
# Rc, hc: rayon et hauteur corps
# Rb, hb: rayon et hauteur ballon
#==============================================================================
def glass1(Rp=0.9,hp=0.2,Rc=0.12,hc=1.,Rb=1.1,hb=1.4,N=20):

    # pied
    p1 = D.polyline([(Rp,0,0),(0.8*Rp,hp*0.5,0),(Rc,hp,0)])

    # corps
    p2 = D.polyline([(Rc,hp,0),(Rc,hp+hc,0)])

    # ballon
    p3 = D.polyline([(Rc,hp+hc,0),(Rb,hp+hc,0),(Rb,hp+hc+hb*0.9,0),(0.9*Rb,hp+hc+hb,0)])
    s3 = D.spline(p3)

    a = T.join([p1,p2,s3])
    a = D.axisym(a, (0,0,0), (0,1,0), 360., N)
    a = C.convertArray2Hexa(a)
    a = G.close(a)
    a = T.rotate(a, (0,0,0), (1,0,0), 90.)
    return a
