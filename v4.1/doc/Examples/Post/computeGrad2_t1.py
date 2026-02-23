# - computeGrad2 (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

#-----
# 3D
#-----
ni = 30; nj = 40; nk = 20
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],1)
#
ni = 30; nj = 40; nk = 20
m = G.cartNGon((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],12)
#
ni = 30; nj = 40; nk = 2
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],2)
#
ni = 30; nj = 40; nk = 2
m = G.cartNGon((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],22)

#-----
# 2D
#-----
ni = 10; nj = 30; nk = 1
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],3)

ni = 10; nj = 30; nk = 1
m = G.cartNGon((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc,'{Density}=2*{x}+{x}*{y}')
mc = C.extractVars(mc,['Density'])
mc = P.computeGrad2(m,mc)
test.testA([mc],32)
