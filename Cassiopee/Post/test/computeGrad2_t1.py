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
mc = C.initVars(mc, '{Density}=3*{x}+2*{y}+{z}')
mc = C.extractVars(mc, ['Density'])
mc = P.computeGrad2(m, mc)
import Converter.PyTree as Co
z = Co.convertArrays2ZoneNode('Zone', [m, mc])
Co.convertPyTree2File(z, "outx.cgns")
test.testA([mc], 32)

ni = 10; nj = 1; nk = 30
m = G.cartNGon((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc, '{Density}=3*{x}+2*{y}+{z}')
mc = C.extractVars(mc, ['Density'])
mc = P.computeGrad2(m, mc)
test.testA([mc], 33)

ni = 1; nj = 10; nk = 30
m = G.cartNGon((0,0,0), (1,10./(nj-1),1), (ni,nj,nk))
mc = C.node2Center(m)
mc = C.initVars(mc, '{Density}=3*{x}+2*{y}+{z}')
mc = C.extractVars(mc, ['Density'])
mc = P.computeGrad2(m, mc)
test.testA([mc], 34)
