# - subzone faces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test

# -- STRUCT
# 1D
N = 51
d = G.cart((0,0,0), (1,1,1), (N,1,1))
facesL = [i for i in range(N)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 7)

# 2D
N = 10; ni1 = N-1
d = G.cart((0,0,0), (1,1,1), (N,N,1))
facesL = [0, ni1, ni1+1, 2*ni1+1, 112, 116, 117]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 8)

# 3D
N = 10; ni1 = N-1
d = G.cart((0,0,0), (1,1,1), (N,N,N))
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
facesL = []
for k in range(ni1):
    for j in range(ni1):
        facesL.append(j*N+k*N*ni1)
d2 = T.subzone(d, facesL, type='faces')
test.testT(d2, 9)


# -- NGON, api 1
# 3D
N = 51
d = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=1)
facesL = [i for i in range(2, N*N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 1)

# 2D
N = 11
d = G.cartNGon((0,0,0), (1,0,1), (N,1,N), api=1)
facesL = [i for i in range(1, N+1)] + [N*(N-1) + i for i in range(1, 2*N-1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateX}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 10)

# 1D
N = 11
d = G.cartNGon((0,0,0), (1,0,0), (N,1,1), api=1)
facesL = [i for i in range(1, 5)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateX}')
#d = T.subzone(d, facesL, type='faces')  # not implemented yet, output must be NODE
#C.convertPyTree2File(d, 'out.cgns'); exit()
#test.testT(d, 11)

# -- NGON, api 3
# 3D
N = 11
d = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=3)
facesL = [i for i in range(1, N*N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 12)

# 2D
N = 11
d = G.cartNGon((0,0,0), (1,0,1), (N,1,N), api=3)
facesL = [i for i in range(1, N+1)] + [N*(N-1) + i for i in range(1, 2*N-1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateX}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 13)

# 1D
N = 11
d = G.cartNGon((0,0,0), (1,0,0), (N,1,1), api=3)
facesL = [i for i in range(1, 5)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateX}')
#d = T.subzone(d, facesL, type='faces')  # not implemented yet, output must be NODE
#C.convertPyTree2File(d, 'out.cgns'); exit()
#test.testT(d, 14)


# -- BE
# TETRA
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,N,N))
facesL = [i for i in range(2, N*N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 2)

# PYRA
N = 11
indices = []
d = G.cartPyra((0,0,0), (1,1,1), (N,N,N))
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
P.exteriorFaces(d, indices=indices)
nf = len(indices[0])
indices = indices[0][:nf//2]  # first half of all exterior faces
d = T.subzone(d, indices, type='faces')
#C.convertPyTree2File(d, 'out.cgns'); exit()
test.testT(d, 15) # TODO known bug centers:G - same as in T.join

# PENTA
N = 11
indices = []
d = G.cartPenta((0,0,0), (1,1,1), (N,N,N))
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
P.exteriorFaces(d, indices=indices)
nf = len(indices[0])
indices = indices[0][:nf//2]  # first half of all exterior faces
d = T.subzone(d, indices, type='faces')
#C.convertPyTree2File(d, 'out.cgns'); exit()
test.testT(d, 16)

# HEXA
N = 51
d = G.cartHexa((0,0,0), (1,1,1), (N,N,N))
facesL = [i for i in range(2, N*N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
#C.convertPyTree2File(d, 'out.cgns'); exit()
test.testT(d, 3)  # TODO bug centers:G

# QUAD
N = 51
d = G.cartHexa((0,0,0), (1,1,1), (N,N,1))
facesL = [i for i in range(1, N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 4)

# TRI
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,N,1))
facesL = [i for i in range(1, N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 5)

# BAR
N = 51
d = G.cartTetra((0,0,0), (1,1,1), (N,1,1))
facesL = [i for i in range(1, N+1)]
C._initVars(d, '{F}={CoordinateX}')
C._initVars(d, '{centers:G}={centers:CoordinateY}')
d = T.subzone(d, facesL, type='faces')
test.testT(d, 6)

## -- ME
# 2D: tri-quad-tri
indices = []
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
a = T.join([a, b, c])
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
indices = [i for i in range(1, 500)]
a = T.subzone(a, indices, type='faces')
#C.convertPyTree2File(a, 'out.cgns'); exit()
test.testT(a, 20)

# 3D: pyra - penta - hexa
indices = []
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = T.join([a, b, c])
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
P.exteriorFaces(a, indices=indices)
nf = len(indices[0])
indices = indices[0][:nf//2]  # first half of all exterior faces
a = T.subzone(a, indices, type='faces')
#C.convertPyTree2File(a, 'out.cgns'); exit()
test.testT(a, 21)
