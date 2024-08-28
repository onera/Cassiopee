# - extractMesh
# 1-  maillage d'interpolation cartesien tetra
# 2-   "             "         tetra+struct
# 3-  meme chose mais a l'ordre 3
# 4-  meme chose mais maillage d'extraction tri

import Converter as C
import Post as P
import Generator as G
import KCore.test as test
import math

tol = 1.e-6

# Create a function
def F(x,y,z): return math.cos(x)

# Maillage en noeuds
ni = 11; nj = 11; nk = 11;
m = G.cartTetra((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))

# init by function
m = C.initVars(m, 'F', F, ['x','y','z'])

# Cree un maillage d'extraction
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (20, 20, 1))

# Extrait la solution sur le maillage d'extraction
a2 = P.extractMesh([m], a, order=2, tol=tol)
test.testA([m,a2], 1)

# test structure
# Maillage en noeuds
ni2 = 20; nj2 = 15; nk2 = 12
m2 = G.cart((0.5,0.1,0), (1./(ni2-1),1./(nj2-1),1./(nk2-1)), (ni2,nj2,nk2))

# init by function
m2 = C.initVars(m2, 'F', F, ['x','y','z'])

#Extrait la solution sur le maillage d'extraction a partir de
# grilles structurees et non structurees
# ordre 2

a3 = P.extractMesh([m,m2], a, order=2, tol=tol)
test.testA([m,m2,a3],2)

# test ordre 3
a4 = P.extractMesh([m,m2], a, order=3, tol=tol)
test.testA([m,m2,a4], 3)

# test maillage d'extraction non structure
b = C.convertArray2Tetra(a)
a5 = P.extractMesh([m,m2], b, order=2, tol=tol)
test.testA([m,m2,a5], 4)

a6 = P.extractMesh([m,m2], [b], order=2, tol=tol)
test.testA(a6,5)
