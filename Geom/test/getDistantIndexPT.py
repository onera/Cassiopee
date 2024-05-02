# - getDistantIndex (pyTree)-
import Geom.PyTree as D

a = D.line((0.,0.,0.), (1.,0.,0), 100)
print('distant Index: %d.'%D.getDistantIndex(a, 25, 0.2))
