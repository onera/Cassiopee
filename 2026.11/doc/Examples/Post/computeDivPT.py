# - computeDiv (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, '{fldX}=2*{CoordinateX}+{CoordinateX}*{CoordinateY}')
m = C.initVars(m, '{fldY}=4*{CoordinateY}')
m = C.initVars(m, '{fldZ}={CoordinateX}+{CoordinateZ}*{CoordinateZ}')
m = P.computeDiv(m, 'fld')
C.convertPyTree2File(m, 'out.cgns')
