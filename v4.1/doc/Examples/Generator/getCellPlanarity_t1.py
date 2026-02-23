# - getCellPlanarity (array) -
import KCore.test as test
import Generator as G
import Geom as D
import Converter as C

a = D.sphere( (0,0,0), 1., 10)
p = G.getCellPlanarity(a)
p = C.center2Node(p)
a = C.addVars([a, p])
test.testA([a],1)
