# - getNPts (array2) -
import Converter as C
import KCore.test as test

a = C.array('x,y,z', 4,4,4, api=2)
npts1 = C.getNPts(a)

a = C.array('x,y,z', 4,4,'HEXA', api=2)
npts2 = C.getNPts(a)

test.testO([npts1,npts2], 1)
