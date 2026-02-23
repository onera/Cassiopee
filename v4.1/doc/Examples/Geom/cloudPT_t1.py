# - cloud (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test
import numpy as np

na, nb = 100, 50
xa = np.arange(na)/na
ya = xa
za = np.zeros_like(xa)
xb = -np.arange(1,nb+1)/na
yb = xb
zb = xb

a = D.cloud((xa,ya,za))
b = D.cloud((xb,yb,zb))
t = C.newPyTree(['Base',1,a,b])
test.testT(t, 1)

test.writeCoverage(100)
