# - getLength (pyTree)-
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

# test getLength structure
a = D.line((0,0,0), (1,0,0))
l = D.getLength(a)
test.testO(l, l)

# test getLength BAR
a2 = C.convertArray2Tetra(a)
l = D.getLength(a2)
test.testO(l, 2)

# test getLength par lots structure
b = D.line((10,0,0), (1,0,0))
l = D.getLength([a,b])
test.testO(l, 3)

# test getLength par lots non structure
b2 = C.convertArray2Tetra(b)
l = D.getLength([a2,b2])
test.testO(l, 4)

# test getLength par lots mixte
l = D.getLength([a,b2])
test.testO(l, 5)
