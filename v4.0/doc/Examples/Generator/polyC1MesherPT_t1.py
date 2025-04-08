# - polyC1Mesher (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Transform.PyTree as T

# Definition de la geometrie
a = D.naca(12., 101)
a = T.reorder(a, (-1,2,3))

h = 0.1; hp = 0.001; density = 100.
res = G.polyC1Mesher(a, h, hp, density)
t = C.newPyTree(['Base']); t[2][1][2] += res[0]
test.testT(t, 1)
#
# test sans split au bord d'attaque
#
res = G.polyC1Mesher(a, h, hp, density, splitCrit=0.01)
t = C.newPyTree(['Base']); t[2][1][2] += res[0]
test.testT(t,2)
