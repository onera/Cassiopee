# - cart2Cyl (array) -
import Transform as T
import Generator as G
import Converter as C
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360, 1., (360,20,10))
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1))
C.convertArrays2File(a, 'out.plt')
