# - naca (array) -
import Geom as D
import KCore.test as test

a = D.naca(12.)
test.testA([a],1)

a = D.naca(12., N=501)
test.testA([a],2)

a = D.naca('0012', N=301, sharpte=1)
test.testA([a],3)

a = D.naca('23012', N=301, sharpte=1)
test.testA([a],4)

a = D.naca('0008-45', N=301, sharpte=1)
test.testA([a],5)

test.writeCoverage(100)
