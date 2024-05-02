# - nearestFaces (array) -
# doit etre exact
import Converter as C
import Generator as G
N = 50
# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cart((100000,0,0), (1,1.e-8,1), (N//2,N//2,N//2))
hook = C.createHook(a, function='faceCenters')
faces,d = C.nearestFaces(hook,b)
ret = d.sum()
#if ret != 0.: print('nearestFaces (Struct) FAILED: FPU is not correct [Struct/FACES]. Check compilation options.')
# NGON
a = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N//2,N//2,N//2))
hook = C.createHook(a, function='faceCenters')
faces,d = C.nearestFaces(hook,b)
ret = d.sum()
#if ret != 0.: print('nearestFaces (BE) FAILED: FPU is not correct [BE/FACES]. Check compilation options.')
# NGON
a = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N//2,N//2,N//2))
hook = C.createHook(a, function='faceCenters')
faces,d = C.nearestFaces(hook,b)
ret = d.sum()
#if ret != 0.: print('nearestFaces (NGON) FAILED: FPU is not correct [NGON/FACES]. Check compilation options.')
