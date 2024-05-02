# - nearestNodes (array) -
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

def nearest(a):
    f = G.cartHexa((2,0,0),(1,1,1),(3,1,1))
    f = T.translate(f,(0.15,0.,0.))
    if isinstance(a[0], list):
        a=a[0]; f = [f,f]
    hook = C.createHook(a, function='nodes')
    nodes = C.nearestNodes(hook, f)
    C.freeHook(hook)
    return nodes

test.stdTestA(nearest)
