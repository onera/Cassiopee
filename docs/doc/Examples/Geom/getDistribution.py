# - getDistribution (array) -
import Geom as D
import Converter as C

Foil = D.naca(12., N=49)
a = D.getDistribution(Foil)
C.convertArrays2File(a, 'out.plt')
