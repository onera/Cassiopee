# - createHook (array) -
import Converter as C
import Generator as G
import KCore.test as test

#test.stdTestA(C.createHook, 'extractMesh')
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
hook = C.createHook([a], function='extractMesh')
#test.testO(hook,1)
#
a = G.cart( (0,0,0), (1,1,1), (10,10,1) )
hook = C.createHook([a], function='extractMesh')
test.testO(hook,2)

a = G.cart( (0,0,0), (1,1,1), (10,10,1) )
hook = C.createHook([a], function='adt')
test.testO(hook,3)
