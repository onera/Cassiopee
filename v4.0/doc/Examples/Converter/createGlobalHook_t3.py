# - createGlobalHook (array) -
import Converter as C
import KCore.test as test

test.stdTestA(C.createGlobalHook, 'elementCenters',1)
test.stdTestA(C.createGlobalHook, 'elementCenters',0)
