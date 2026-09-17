# - subzone elements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# 3D NGON
N = 51
d = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=1)
eltsL = [i for i in range(N*N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 1)

# 3D Tetra
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,N,N))
eltsL = [i for i in range(N*N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 2)

# 3D Hexa
N = 51
d = G.cartHexa((0,0,0), (1,1,1), (N,N,N))
eltsL = [i for i in range(N*N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 3)

# 2D quad
N = 51
d = G.cartHexa((0,0,0), (1,1,1), (N,N,1))
eltsL = [i for i in range(N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 4)

# 2D TRI
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,N,1))
eltsL = [i for i in range(N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 5)

# 1D BAR
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,1,1))
eltsL = [i for i in range(N//2)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 6)

# 3D NGON
N = 11
d = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=3)
eltsL = [i for i in range(N*N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
test.testT(d, 11)

# 1D ME
N = 21
a = G.cartTetra((0,0,0), (1,0,0), (N,1,1))
b = G.cartTetra((N-1,0,0), (0,1,0), (1,N,1))
c = G.cartTetra((0,N-1,0), (1,0,0), (N,1,1))
a = C.mergeConnectivity([a, b, c])
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
eltsL = [i+j*20 for j in range(3) for i in range(N//2)]  # 1st half of each BAR
a = T.subzone(a, eltsL, type='elements')
test.testT(a, 61)

# -- 2D ME: tri-quad-tri
#                |
#               tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
eltsL = (
    [i for i in range(36)] + [72+i for i in range(18)]
    + [108+i for i in range(36)] + [180+i for i in range(36)]
)
a = T.subzone(a, eltsL, type='elements')
test.testT(a, 41)

# 3D ME: hexa - pyra
#          |      |
#        pyra   hexa
a = G.cartHexa((0.,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c, d], None)
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
eltsL = (
    [i for i in range(32)] + [64+i for i in range(192)]  # 1st half of these conns
    + [640+i for i in range(192)] + [864+i for i in range(32)]  # 2nd half of these conns
)
a = T.subzone(a, eltsL, type='elements')
test.testT(a, 32)
