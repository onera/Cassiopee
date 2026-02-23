# - node2ExtCenter (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

ni = 30; nj = 40; nk = 1
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
ac = C.node2ExtCenter(a)
t = C.newPyTree(['Base',2]); t[2][1][2].append(ac)
C.convertPyTree2File(t, 'out.cgns')
