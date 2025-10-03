# - computeGrad (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

import numpy
def F(x,y,z): return (numpy.sin(y), numpy.sin(z), numpy.sin(x))
# F = (sin(y), sin(z), sin(x))
# rot(F) = (-cos(z), -cos(x), -cos(y))

N = 10

#-------------------------------
# 2D
#-------------------------------

# BE (TRI)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,2)

# BE (QUAD)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,3)

# ME (QUAD+TRI)
a = G.cartHexa((0.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
b = G.cartTetra((1.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,4)

#-------------------------------
# 3D
#-------------------------------

# NGONv3
z = G.cartNGon((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N), api=1)
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,5)

# NGONv4
z = G.cartNGon((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N), api=3)
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,6)

# BE (TETRA)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,7)

# BE (PENTA)
z = G.cartPenta((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,8)

# BE (PYRA)
z = G.cartPyra((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,9)

# BE (HEXA)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,10)

# ME (HEXA+TETRA)
a = G.cartHexa((0.,0.,0.), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
b = G.cartTetra((1.,0.,0.), (1./(N-1),1./(N-1),1./(N-1)), (N,N,N))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, ['F1', 'F2', 'F3'], F, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = P.computeNormCurl(z, ['F1', 'F2', 'F3'])
test.testT(z,11)
