# - nearestElements (array) -
# doit etre exact
import Converter as C
import Generator as G
import Transform as T

eps = 1.e-45
N = 50
# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = T.subzone(a, (1,1,1), (N,N//2,N//2))
hook = C.createHook(a, function='elementCenters')
faces = C.identifyElements(hook,b,tol=eps)
ret = faces[faces<0.].shape[0]
#if ret != 0.: print 'nearestElements (structured) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.'
# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = T.subzone(a, (1,1,1), (N,N//2,1))
hook = C.createHook(a, function='faceCenters')
faces = C.identifyElements(hook,b,tol=eps)
ret = faces[faces<0.].shape[0]
#if ret != 0.: print 'identifyElements/faces (structured) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.'


# EB
a = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N//2,2))
hook = C.createHook(a, function='elementCenters')
faces = C.identifyElements(hook,b,tol=eps)
ret = faces[faces<0.].shape[0]
#if ret != 0.: print 'identifyElements (EB) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.'
# NGON
a = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N//2,2))
hook = C.createHook(a, function='elementCenters')
faces = C.identifyElements(hook,b,tol=eps)
ret = faces[faces<0.].shape[0]
#if ret != 0.: print 'identifyElements (NGon) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.'
print('done.')
