# - projectAllDirs (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
a = D.sphere6((0,0,0), 1., 20)
b = G.cart((1.1,-0.1,-0.1),(0.03,0.03,0.03), (1,50,50))
n = G.getNormalMap(b)
n = C.center2Node(n)
b = C.addVars([b,n])
c = T.projectAllDirs([b], a, ['sx','sy','sz'])
C.convertArrays2File([b]+c, 'out.plt')
