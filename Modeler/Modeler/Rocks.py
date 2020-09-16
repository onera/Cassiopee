# - Rocks -
import Geom as D
import Transform as T
import Generator as G
import Converter as C
#import Post as P

#==============================================================================
def cobbleStone(hx=1., hy=0.3, hz=0.2, N=8):
    a = D.sphere((0,0,0), R=0.5*hx, N=N)
    #a = C.convertArray2Tetra(a)
    #a = P.refine(a, w=1./8.)
    a = T.scale(a, factor=(1.,hy,hz))
    a = C.convertArray2Tetra(a)
    a = G.close(a)
    return a
