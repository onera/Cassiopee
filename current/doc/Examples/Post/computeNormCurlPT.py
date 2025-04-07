# - computeNormCurl (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

def F(x,y,z): return 12*y*y + 4

ni = 30; nj = 40; nk = 3
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'F1',F,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.addVars(m,'F2'); m = C.addVars(m,'F3')

varname = ['F1','F2','F3']
m = P.computeNormCurl(m, varname)
C.convertPyTree2File(m, 'out.cgns')
