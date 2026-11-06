# - selectCells (pyTree) -
import Converter.PyTree   as C
import Converter.Internal as Internal
import Generator.PyTree   as G
import Post.PyTree        as P
import KCore.test         as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False


# CAS 1 : tag au centre - existance champs aux centres et aux noeuds
# ==================================================================
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}={centers:CoordinateZ}+{centers:CoordinateX}')
C._initVars(a,'{var2}={CoordinateX}')
C._initVars(a,'{var}=({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '{centers:var1}<3.',strict=1)
test.testT(a,1)


# CAS 2 : tag au centre - existance champ aux noeuds uniquement
# =============================================================
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}={centers:CoordinateZ}+{centers:CoordinateX}')
C._initVars(a,'{var2}={CoordinateX}')
C._initVars(a,'{var}=({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '{centers:var1}<3.',strict=1)
test.testT(a,2)
