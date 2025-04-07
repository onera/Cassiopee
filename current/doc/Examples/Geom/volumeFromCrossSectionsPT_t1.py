# - volumeFromCrossSection (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test
contours = []
for z in [0.,1.]:
    contours.append(D.circle((0,0,z),1.,N=15))
vol = D.volumeFromCrossSections(contours)
t = C.newPyTree(['Base',2,vol])
test.testT(t, 1)
