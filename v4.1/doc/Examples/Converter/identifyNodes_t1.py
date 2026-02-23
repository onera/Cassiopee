# - identifyNodes (array) -
import Converter as C
import Generator as G
import KCore.test as test

def identify(a):
    f = G.cartHexa((2,0,0),(1,1,1),(3,1,1))
    if isinstance(a[0], list):
        a = a[0]; f = [f,f]
    hook = C.createHook(a, function='nodes')
    nodes = C.identifyNodes(hook, f)
    C.freeHook(hook)
    return nodes

test.stdTestA(identify)
