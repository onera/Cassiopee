# - identifyNodes (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
hook, indir = C.createGlobalHook([a,b], function='nodes', indir=1)
offset = [0, C.getNPts(a), C.getNPts(b)]

f = D.point((13,3,3))
nodes = C.identifyNodes(hook, f)
ind = nodes[0]
print('Le premier point de f a pour indice', ind-offset[indir[ind]], 'sur la zone', indir[ind])
#>> Le premier point de f a pour indice 332 sur la zone 1
