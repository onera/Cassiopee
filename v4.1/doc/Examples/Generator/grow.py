# - grow (array) -
import Converter as C
import Generator as G
import Geom as D

a = D.sphere( (0,0,0), 1., 50 )
n = G.getNormalMap(a)
n = C.center2Node(n); n[1] = n[1]*100.
b = G.grow(a, n)
C.convertArrays2File([b], 'out.plt')
