# subzone faces (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

N = 51
d = G.cartNGon((0,0,0), (1,1,1),(N,N,N))
facesL=[]
for i in range(1,N*N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,1)
# 3D Tetra
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,N,N))
facesL=[]
for i in range(1,N*N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,2)
#  3D Hexa
N = 51
d = G.cartHexa((0,0,0), (1,1,1),(N,N,N))
facesL=[]
for i in range(1,N*N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,3)

# 2D quad
N = 51
d = G.cartHexa((0,0,0), (1,1,1),(N,N,1))
facesL=[]
for i in range(N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,4)

# 2D TRI
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,N,1))
facesL=[]
for i in range(N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,5)

# 1D BAR
N = 51
d = G.cartTetra((0,0,0), (1,1,1),(N,1,1))
facesL=[]
for i in range(N): facesL.append(i+1)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,6)

# 1D struct
N = 51
d = G.cart((0,0,0), (1,1,1),(N,1,1))
facesL=[]
for i in range(N): facesL.append(i)
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,7)

# 2D struct
N = 10; ni1 = N-1
d = G.cart((0,0,0), (1,1,1),(N,N,1))
facesL=[0,ni1,ni1+1,2*ni1+1,112,116,117]
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d,8)
# 3D struct
N = 10; ni1 = N-1
d = G.cart((0,0,0), (1,1,1),(N,N,N))
C._initVars(d,'{F}={CoordinateX}')
C._initVars(d,'{centers:G}={centers:CoordinateY}')
facesL=[]
for k in range(ni1):
    for j in range(ni1):
        facesL.append(j*N+k*N*ni1)
d2 = T.subzone(d, facesL, type='faces')
test.testT(d2,9)
