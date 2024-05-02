# - computeGrad (array) -
import Converter as C
import Post as P
import Generator as G

ni = 30; nj = 40; nk = 2
def F(x):
    return 12*x*x + 4

def DF(x):
    return 12.*x*2

# Maillage en noeuds 3d
m0 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.addVars(m0,'F')
m = C.initVars(m,'F',F,['x'])
#
p = P.computeGrad(m, 'F')
c = C.node2Center(m0)
c = C.addVars(c,'DF')
c = C.initVars(c,'DF',DF,['x'])

c = C.addVars([c,p])
C.convertArrays2File([c],"out1.plt","bin_tp")
#
# CAS 2D
import Transform as T
# Maillage en noeuds
m0 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1.), (ni,nj,nk))
m = C.addVars(m0,'F')
m = C.initVars(m,'F',F,['x'])

p = P.computeGrad(m, 'F')
c = C.node2Center(m)
c = C.addVars(c,'DF')
c = C.initVars(c,'DF',DF,['x'])

c = C.addVars([c,p])
C.convertArrays2File([c],"out2.plt","bin_tp")
