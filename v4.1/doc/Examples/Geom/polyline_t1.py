# - polyline -
import Geom as D
import KCore.test as test

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)])
test.testA([a],1)
