# - collapse (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

hi = 0.1; hj = 0.01; hk = 1.
ni = 20; nj = 2; nk = 1
a = G.cartTetra((0.,0.,0.),(hi,hj,hk),(ni,nj,nk))
#a = C.initVars(a, 'F', 2.); #a = C.initVars(a, 'centers:G', 1.)
b = T.collapse(a)
t = C.newPyTree(['BAR',1]); t[2][1][2].append(b)
test.testT(t, 1)

a = G.cartTetra((0.,0.,0.),(hi,hj,hk),(ni,nj,nk))
T._collapse(a)
test.testT(t, 2)
