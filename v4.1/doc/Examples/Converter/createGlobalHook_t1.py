# - createGlobalHook (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = G.cartNGon((9,0,0), (1,1,1), (10,10,10))
hook = C.createGlobalHook([a,b], function='faceCenters')
test.testO(hook, 1)
