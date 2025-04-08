# - deformNormals (array) -
import Converter as C
import Geom as D
import Transform as T
import KCore.test as test

# Structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.initVars(a,'{F}={x}+{y}**2')
alpha = C.initVars(a,'alpha', 0.5)
alpha = C.extractVars(alpha,['alpha'])
b = T.deformNormals(a,alpha)
test.testA([a,b], 1)

# Non structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.convertArray2Hexa(a)
a = C.initVars(a,'{F}={x}+{y}**2')
alpha = C.initVars(a,'alpha', 0.5)
alpha = C.extractVars(alpha,['alpha'])
b = T.deformNormals(a,alpha)
test.testA([a,b], 2)

# liste de zones
a = D.sphere6( (0,0,0), 1., 10)
a = C.initVars(a,'{F}={x}+{y}**2')
alpha = C.initVars(a,'alpha', 0.5)
alpha = C.extractVars(alpha,['alpha'])
b = T.deformNormals(a,alpha)
test.testA(a+b, 3)
