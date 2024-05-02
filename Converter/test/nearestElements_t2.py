# - nearestElements (array) -
# doit etre exact
import Converter as C
import Generator as G
import Transform as T

N = 50
# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = T.subzone(a, (1,1,1), (N,N//2,N//2))
hook = C.createHook(a, function='elementCenters')
faces,d = C.nearestElements(hook,b)
ret = d.sum()
if ret != 0.: print('1. nearestElements (structured) FAILED: FPU is not correct [STRUCT/ELTS]. Check compilation options.')

# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = T.subzone(a, (1,1,1), (N,N//2,1))
hook = C.createHook(a, function='faceCenters')
faces,d = C.nearestElements(hook,b)
ret = d.sum()
if ret != 0.: print('2. nearestElements/faces (structured) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.')

# EB
a = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N//2,2))
hook = C.createHook(a, function='elementCenters')
faces,d = C.nearestElements(hook,b)
ret = d.sum()
if ret != 0.: print('3. nearestElements (EB) FAILED: FPU is not correct [STRUCT/ELTS]. Check compilation options.')
# NGON
a = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N//2,2))
hook = C.createHook(a, function='elementCenters')
faces,d = C.nearestElements(hook,b)
ret = d.sum()
if ret != 0.: print('4. nearestElements (NGon) FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.')
print('done.')
