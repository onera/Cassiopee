import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree      as P
import KCore.test       as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

# Struct (strict = 0)
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, '{var1}>15.')
test.testT(a,1)

# Struct (strict = 1)
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, '{var1}>15.', strict=1)
test.testT(a,2)

# Struct (formula)
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(a,3)


# NGon
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{centers:var1}=1.')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
C._initVars(a,'{var3}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
a = P.selectCells(a, '{var3}>10.', strict=1)
test.testT(a,4)

# Non struct
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
a = C.convertArray2Hexa(a)
C._initVars(a,'{Density}={CoordinateX}')
C._initVars(a,'{var1}={centers:CoordinateX}+2.*{centers:CoordinateY}+{centers:CoordinateZ}')
C._initVars(a,'{centers:var2}={centers:CoordinateX}')
a = P.selectCells(a, '{var1}>15.')
test.testT(a,5)
