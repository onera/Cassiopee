# - mmgs (array) -
import Geom as D
import Generator as G
import Converter as C

a = D.sphere6( (0,0,0), 1., N=20, ntype='TRI' )

# Optimisation
b = G.mmgs(a, optim=1)
C.convertArrays2File(b, 'out.plt')

# Remaillage avec parametres
b = G.mmgs(a, hausd=0.01, hmax=0.1)
C.convertArrays2File(b, 'out.plt')

# Raffinement avec sizemap
a = C.initVars(a, 'sizemap=0.2*abs({x})+0.05')
b = G.mmgs(a, hausd=10.)
C.convertArrays2File(b, 'out.plt')
