# - mmgs (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = D.sphere6( (0,0,0), 1., N=20, ntype='TRI' )

# Optimisation
b = G.mmgs(a, optim=1)
C.convertPyTree2File(b, 'out1.cgns')

# Remaillage avec parametres
b = G.mmgs(a, hausd=0.01, hmax=0.1)
C.convertPyTree2File(b, 'out2.cgns')

# Raffinement avec sizemap
a = C.initVars(a, '{sizemap}=0.2*abs({CoordinateX})+0.05')
b = G.mmgs(a, hausd=10.)
C.convertPyTree2File(b, 'out3.cgns')
