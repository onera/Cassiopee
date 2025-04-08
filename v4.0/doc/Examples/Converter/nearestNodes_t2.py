# - nearestNodes (array) -
# Permet de tester que Converter a ete bien installe
import Converter as C
import Generator as G
import Transform as T

#==============================================================================
# Identify nodes structure - doit etre exact!!
a = G.cart((100000,0,0), (1,1.e-8,1), (30,30,30))
b = T.subzone(a, (1,1,1), (30,30,1))

hook = C.createHook(a, function='nodes')
ids, d = C.nearestNodes(hook, b)
ret = d.sum()
# ret doit etre 0 exact (binaire)
#if ret != 0.: print 'FAILED: FPU is not correct [STRUCT/NODES]. Check compilation options.'

#==============================================================================
# Identify faces structure - doit etre exact!!
a = G.cart((100000,0,0), (1,1.e-8,1), (30,30,30))
b = T.subzone(a, (1,1,1), (30,30,1))
hook = C.createHook(a, function='faceCenters')
ids, d = C.nearestElements(hook, b)
ret = d.sum()
#if ret != 0.: print 'FAILED: FPU is not correct [STRUCT/FACES]. Check compilation options.'

#==============================================================================
# Identify nodes on NGON - doit etre exact!!
a = G.cart((100000,0,0), (1,1.e-8,1), (30,30,30))
b = T.subzone(a, (1,1,1), (30,30,1))

a = C.convertArray2NGon(a)
b = C.convertArray2NGon(b)
hook = C.createHook(a, function='nodes')

ids, d = C.nearestNodes(hook, b)
ret = d.sum()
#if ret != 0.: print 'FAILED: FPU is not correct [NGON/NODES]. Check compilation options.'

#==============================================================================
# Identify faces sur NGON - doit etre exact!!
a = G.cart((100000,0,0), (1,1.e-8,1), (10,10,10))
b = T.subzone(a, (1,1,1), (10,10,1))

a = C.convertArray2NGon(a)
b = C.convertArray2NGon(b)
hook = C.createHook(a, function='faceCenters')

ids, d = C.nearestElements(hook, b)
ret = d.sum()
#if ret != 0.: print 'FAILED: FPU is not correct [NGON/FACES]. Check compilation options.'
