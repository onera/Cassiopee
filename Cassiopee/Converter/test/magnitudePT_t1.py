# - magnitude (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

a = D.sphere((0,0,0), 1., 50 )
a = G.getNormalMap(a)
a = C.magnitude(a, ['centers:sx','centers:sy','centers:sz'])
t = C.newPyTree(['Base',2,a])
test.testT(t, 1)
