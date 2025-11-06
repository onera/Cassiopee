# - exteriorElts (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test
import numpy

def densityField(x, y, z):
    return numpy.sqrt((x-2)**2 + (y-1)**2 + z**2)

def pressureField(x, y, z):
    return numpy.logical_and(x < 12,  y + z >= 1.5).astype(float)

def _addVars(a):
    C._initVars(
        a, 'Density', densityField,
        ["CoordinateX", "CoordinateY", "CoordinateZ"],
        isVectorized=True
    )
    C._initVars(
        a, 'centers:pressure', pressureField,
        ["centers:CoordinateX", "centers:CoordinateY", "centers:CoordinateZ"],
        isVectorized=True
    )

# BE - 3D
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorElts(a)
test.testT(b, 1)

# ME - 2D (disconnected)
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
b = G.cartHexa((10,-2,0), (1,1,1), (8,6,1))
c = G.cartTetra((20,-3,0), (1,1,1), (5,9,1))
a = C.mergeConnectivity([a, b, c], None, boundary=0)
_addVars(a)
b = P.exteriorElts(a)
test.testT(b, 2)

# ME - 2D (all adjacent)
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
b = G.cartHexa((9,-2,0), (1,1,1), (8,6,1))
c = G.cartTetra((16,-3,0), (1,1,1), (5,9,1))
a = C.mergeConnectivity([a, b, c], None, boundary=0)
_addVars(a)
b = P.exteriorElts(a)
test.testT(b, 3)

# ME - 3D
a = G.cartTetra((0,0,-4), (1,1,1), (10,10,8))
b = G.cartPyra((9,-2,-4), (1,1,1), (8,6,5))
c = G.cartPenta((16,-3,-4), (1,1,1), (5,9,7))
d = G.cartHexa((20,-2,-4), (1,1,1), (6,8,6))
a = C.mergeConnectivity([a, b, c, d], None, boundary=0)
_addVars(a)
b = P.exteriorElts(a)
test.testT(b, 4)


"""
# NGon - 2D
a = G.cartNGon((0,0,0), (1,1,1), (4,4,1), api=1)
_addVars(a)
b = P.exteriorElts(a)
C.convertPyTree2File(b, "out.cgns")
test.testT(b, 10)

# NGon - 3D
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
_addVars(a)
b = P.exteriorElts(a)
test.testT(b, 11)

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
_addVars(a)
b = P.exteriorElts(a)
test.testT(b, 12)
"""
