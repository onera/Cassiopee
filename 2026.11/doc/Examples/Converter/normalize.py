# - normalize (array) -
import Converter as C
import Generator as G
import Geom as D

a = D.sphere((0,0,0), 1., 50)
n = G.getNormalMap(a)
n = C.center2Node(n)
n[1] = n[1]*10
n = C.normalize(n, ['sx','sy','sz'])
a = C.addVars([a, n])
C.convertArrays2File(a, 'out.plt')
