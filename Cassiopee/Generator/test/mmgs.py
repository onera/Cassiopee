# - mmgs (array) -
import Geom as D
import Generator as G
import Converter as C

a = D.sphere6( (0,0,0), 1., N=20, ntype='TRI' )

# Optimisation
b = G.mmgs(a, optim=1)
C.convertArrays2File(b, 'out1.plt')

# Remaillage isotrope avec hmin/hmax/hausd
b = G.mmgs(a, hausd=0.001, hmin=0.001, hmax=0.1, grow=1.3)
C.convertArrays2File(b, 'out2.plt')

# Raffinement aniso avec hmin/hmax/hausd
b = G.mmgs(a, hausd=0.001, hmin=0.001, hmax=0.1, grow=1.3, anisotropy=1)
C.convertArrays2File(b, 'out3.plt')

# Raffinement isotrope avec sizemap scalar
a = C.initVars(a, 'sizemap=0.2*abs({x})+0.05')
b = G.mmgs(a)
C.convertArrays2File(b, 'out4.plt')
