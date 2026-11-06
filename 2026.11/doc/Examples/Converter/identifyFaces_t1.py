# - identifyFaces (array) -
import Converter as C
import Post as P
import KCore.test as test

def identify(a):
    a = C.convertArray2NGon(a)
    f = P.exteriorElts(a)
    if isinstance(a[0], list): a=a[0]
    hook = C.createHook(a, function='faceCenters')
    faces = C.identifyFaces(hook, f)
    C.freeHook(hook)
    return faces

test.stdTestA(identify)
