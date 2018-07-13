# - createGlobalHook (array) -
import Converter as C
import Generator as G
import KCore.test as test

test.stdTestA(C.createGlobalHook, 'nodes',1)
test.stdTestA(C.createGlobalHook, 'nodes',0)
