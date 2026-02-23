# - volumeFromCrossSections -
import Geom as D
import KCore.test as test
contours = []
for z in [0.,1.]:
    contours.append(D.circle((0,0,z),1.,N=15))
vol = D.volumeFromCrossSections(contours)
test.testA([vol],1)
