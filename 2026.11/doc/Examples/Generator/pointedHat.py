# - pointedHat (array) -
import Geom as D
import Generator as G
import Converter as C
c = D.circle( (0,0,0), 1., 360., 0., 100)
surf = G.pointedHat(c,(0.,0.,1.))
C.convertArrays2File([surf], 'out.plt')
