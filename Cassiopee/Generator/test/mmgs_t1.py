# - mmgs (array) -
import Geom as D
import Generator as G
import Converter as C
import KCore.test as test

a = D.sphere6((0,0,0), 1., N=20, ntype='TRI')

# Optimisation
b = G.mmgs(a, optim=1)
test.testA(b, 1)

# Remaillage isotrope hmin/hmax/hausd
b = G.mmgs(a, hausd=0.01, hmax=0.1)
test.testA(b, 2)

# Raffinement isotrope sizemap
a = C.initVars(a, 'sizemap=0.2*abs({x})+0.05')
b = G.mmgs(a)
test.testA(b, 3)
