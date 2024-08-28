# - getDistribution (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

Foil = D.naca(12., N=49)
a = D.getDistribution(Foil)
C.convertPyTree2File(a, 'out.cgns')
