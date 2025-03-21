import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree      as P
import KCore.test       as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

# Struct
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, '{centers:var1}>15.')
test.testT(a,1)

# NGon
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}=1.')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
C._initVars(a,'{centers:var3}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
a = P.selectCells(a, '{centers:var3}>10.')
test.testT(a,2)

# Non struct
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
a = C.convertArray2Hexa(a)
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, '{centers:var1}>15.')
test.testT(a,3)
