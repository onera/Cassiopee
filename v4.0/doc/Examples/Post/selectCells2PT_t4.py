# - selectCells2 (pyTree) -
# Cas test avec un tag en noeuds
import Converter.PyTree   as C
import Converter.Internal as Internal
import Generator.PyTree   as G
import Post.PyTree        as P
import KCore.test         as test
import Intersector.PyTree as XOR

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

# Test tag aux noeuds - sans champ en centre
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = C.initVars(a, 'tag', F, ['CoordinateX','CoordinateY','CoordinateZ'])
b = P.selectCells2(a, 'tag')
test.testT(a,1)


# Test tag aux noeuds - avec champ en centre
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = C.initVars(a, 'tag', F, ['CoordinateX','CoordinateY','CoordinateZ'])
C._initVars(a,'{centers:var1}={centers:CoordinateZ}+{centers:CoordinateX}')
b = P.selectCells2(a, 'tag')
test.testT(a,2)


# Test tag aux centres - sans champ en centre
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = C.initVars(a, 'centers:tag', F, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
b = P.selectCells2(a, 'centers:tag')
test.testT(a,3)


# Test tag aux centres - avec champ en centre
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

C._initVars(a,'{centers:var1}={centers:CoordinateZ}+{centers:CoordinateX}')
C._initVars(a,'{var2}=1.')
a = C.initVars(a, 'centers:tag', F, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
b = P.selectCells2(a, 'centers:tag')
test.testT(a,4)
