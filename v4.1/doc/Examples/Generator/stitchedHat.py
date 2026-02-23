# - stitchedHat (array) -
import Geom as D
import Generator as G
import Transform as T
import Converter as C

c = D.circle( (0,0,0), 1., 360., 0., 100)
c = T.contract(c, (0,0,0), (0,1,0), (0,0,1), 0.1)
c = G.stitchedHat(c, (0,0,0), 1.e-3)
C.convertArrays2File([c], 'out.plt')
