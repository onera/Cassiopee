# - computeCurl (array) -
import Converter as C
import Post as P
import Generator as G

ni = 30; nj = 40; nk = 2
def F(x,y,z):
    return 12*y*y + 4

def H(x) :
    return 0.

# Maillage en noeuds
m0 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.addVars(m0,'F1')
m = C.initVars(m,'F1',F,['x','y','z'])
m = C.addVars(m,'F2')
m = C.initVars(m,'F2',0.)
m = C.addVars(m,'F3')
m = C.initVars(m,'F3',0.)
#
mc = P.node2Center(m)
varname = ['F1','F2','F3']
p = P.computeCurl(m,varname)
#
p = C.addVars([mc,p])
C.convertArrays2File([p],"out1.plt","bin_tp")

# TEST 2D
ni = 30; nj = 40; nk = 1

# Maillage en noeuds
m0 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.addVars(m0,'F1')
m = C.initVars(m,'F1',F,['x','y','z'])
m = C.addVars(m,'F2')
m = C.initVars(m,'F2',0.)
m = C.addVars(m,'F3')
m = C.initVars(m,'F3',0.)
#
mc = P.node2Center(m)
p = P.computeCurl(m,varname)
#
p = C.addVars([mc,p])
C.convertArrays2File([p],"out2.plt","bin_tp")
