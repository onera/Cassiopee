# - subzone elements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

N = 51
d = G.cartNGon((0,0,0), (1,1,1),(N,N,N))
eltsL=[]
for i in range(N*N): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,1)
# 3D Tetra
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,N,N))
eltsL=[]
for i in range(N*N): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,2)
#  3D Hexa
N = 51
d = G.cartHexa((0,0,0), (1,1,1),(N,N,N))
eltsL=[]
for i in range(N*N): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,3)

# 2D quad
N = 51
d = G.cartHexa((0,0,0), (1,1,1),(N,N,1))
eltsL=[]
for i in range(N): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,4)

# 2D TRI
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,N,1))
eltsL=[]
for i in range(N): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,5)

# 1D BAR
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,1,1))
eltsL=[]
for i in range(N//2): eltsL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d,eltsL, type='elements')
test.testT(d,6)
