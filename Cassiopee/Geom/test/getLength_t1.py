# - getLength -
import Geom as D
import Converter as C
import KCore.test as test

res = C.array('r', 1,1,1)

# test getLength structure
a = D.line((0,0,0), (1,0,0))
l = D.getLength(a) - 1.0
res[1][0,0] = l
test.testA([res], 1)

# test getLength BAR
a2 = C.convertArray2Tetra(a)
l = D.getLength(a2) - 1.0
res[1][0,0] = l
test.testA([res], 2)

# test getLength par lots structure
b = D.line((10,0,0), (1,0,0))
l = D.getLength([a,b])
res[1][0,0] = l
test.testA([res], 3)

# test getLength par lots non structure
b2 = C.convertArray2Tetra(b)
l = D.getLength([a2,b2])
res[1][0,0] = l
test.testA([res], 4)

# test getLength par lots mixte
l = D.getLength([a,b2])
res[1][0,0] = l
test.testA([res], 5)
