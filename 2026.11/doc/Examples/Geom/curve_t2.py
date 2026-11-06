# - curve (array) -
import Geom as D
import KCore.test as test
from Geom.Parametrics import base

a = D.curve(base['circle'], N=100)
test.testA([a],1)
a = D.curve(base['line'], N=100)
test.testA([a],2)
test.writeCoverage(100)
