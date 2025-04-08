# - selectCells (pyTree) -
import Converter.PyTree   as C
import Converter.Internal as Internal
import Generator.PyTree   as G
import Post.PyTree        as P
import KCore.test         as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False


# CAS 1 : tag au noeud - existance champs aux centres et aux noeuds
# =================================================================
a = G.cartNGon( (0,0,0), (1,1,1), (21,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}=1.')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
C._initVars(a,'{var}=({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '{CoordinateX}+{CoordinateZ}<3.',strict=1)
test.testT(a,1)


# CAS 2 : Idem que 1 + 'reverse' pour certaines faces
# ===================================================
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}=1.')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
C._initVars(a,'{centers:var3}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{var}=({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2>25.',strict=0)
a = P.selectCells(a, '{CoordinateY}>1.5',strict=0)
test.testT(a,2)


# CAS 3 : tag au noeud + uniquement des champs aux centres
# =========================================================
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{centers:var1}=1.')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
C._initVars(a,'{centers:var3}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2>25.',strict=0)
test.testT(a,3)


# CAS 4 : tag au noeud + uniquement des champs aux noeuds
# =======================================================
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{var1}=1.')
C._initVars(a,'{var2}={centers:CoordinateX}')

a = Internal.createElsaHybrid(a,method=1,methodPE=1)
a = Internal.rmNodesByName(a, ':elsA#Hybrid')

a = P.selectCells(a, '({CoordinateX}-5.)**2+({CoordinateY}-5.)**2+({CoordinateZ}-5.)**2>25.',strict=0)
test.testT(a,4)
