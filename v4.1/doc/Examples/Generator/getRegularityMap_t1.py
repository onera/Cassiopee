# - getRegularityMap (array) -
import Generator as G
import Converter as C
import Geom as D
import Transform as T
import KCore.test as test

# Test 3D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
ac = C.node2Center(a)
reg = G.getRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 1)

# Test 2D structure
msh = D.naca(12., 5001)
# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = G.hyper2D(msh, distrib, "C")
a = T.reorder(a, (-3,2,1))
ac = C.node2Center(a)
reg = G.getRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 2)

# Test 1D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,10,1))
ac = C.node2Center(a)
reg = G.getRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 3)
