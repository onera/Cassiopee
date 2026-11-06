# - conformUnstr (pyTree) -
# Conforming 1 or 2 TRI/BAR together (same type for both operands)
import Generator.PyTree as G
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Geom.PyTree as D
from Geom.Parametrics import base
import Transform.PyTree as T

s1 = D.sphere((0,0,0), 1, N=20)

s2 = D.surface(base['plane'], N=30)
s2 = T.translate(s2, (0.2,0.2,0.2))

s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)

x = XOR.conformUnstr(s1, s2, tol=0.)
C.convertPyTree2File(x, 'out.plt')

c1 = D.circle((0,0,0), 1, N=100)
c2 = D.circle((0.2,0,0), 1, N=50)

c1 = C.convertArray2Tetra(c1); c1 = G.close(c1)
c2 = C.convertArray2Tetra(c2); c2 = G.close(c2)

x = XOR.conformUnstr(c1, c2, tol=0.)
C.convertPyTree2File(x, 'out1.plt')
