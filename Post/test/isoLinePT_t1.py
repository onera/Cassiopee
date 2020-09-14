# - isoLine (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x, y): return x*x+y*y

# Test sur un champ en noeuds
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, 'field', F, ['CoordinateX','CoordinateY'])
iso = P.isoLine(a, 'field', 15.)
test.testT(iso, 1)

# Test sur un champ en centres
b = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
b = C.node2Center(b, ['CoordinateX', 'CoordinateY'])
b = C.initVars(b, 'centers:field', F, ['centers:CoordinateX','centers:CoordinateY'])
b = C.rmVars(b, ['centers:CoordinateX', 'centers:CoordinateY'])
iso = P.isoLine(b, 'centers:field', 15.)
t = C.newPyTree(['Base'])
t[2][1][2] += [b, iso]
#test.testA(iso, 2)
