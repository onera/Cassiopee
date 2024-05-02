# - intersection (array) -
import Intersector as XOR
import Generator as G
import Converter as C
import Geom as D

s1 = D.sphere((0,0,0), 1, N=20)
s2 = D.sphere((0.,1.,0.), 1, N=30)
s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)
x = XOR.intersection(s1, s2, tol=0.)
C.convertArrays2File([x], 'out.plt')
