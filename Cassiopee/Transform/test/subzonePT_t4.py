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

"""
# 3D NGON
N = 11
d = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=3)
eltsL = [i for i in range(N*N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, eltsL, type='elements')
#C.convertPyTree2File(d, "out.cgns")
#test.testT(d, 11)
"""
