# - computeDiv2 (array) -
import Converter as C
import Post as P
import Generator as G
import math

m = G.cartNGon((0,0,0), (1,1,1), (4,4,4))

def Fu(a): return math.cos(a)
def Fv(a): return 4.*a
def Fw(a,b,c): return b*(c**2)

mc = C.node2Center(m)
mc = C.initVars(mc, 'fldX', Fu, ['x'])
mc = C.initVars(mc, 'fldY', Fv, ['y'])
mc = C.initVars(mc, 'fldZ', Fw, ['x', 'y', 'z'])
mv = C.extractVars(mc, ['fldX', 'fldY', 'fldZ'])

p = P.computeDiv2(m, mv) # p is defined on centers
p = C.center2Node(p) # back on initial mesh
p = C.addVars([m, p])
C.convertArrays2File([p], 'out.plt')
