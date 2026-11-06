# - nearestFaces (array) -
import Converter as C
import Post as P
import Transform as T
import KCore.test as test

def nearest(a):
    f = P.exteriorFaces(a)
    f = T.translate(f,(0.15,0.,0.))
    if isinstance(a[0], list): a=a[0]
    hook = C.createHook(a, function='faceCenters')
    faces = C.nearestFaces(hook, f)
    C.freeHook(hook)

test.stdTestA(nearest)
