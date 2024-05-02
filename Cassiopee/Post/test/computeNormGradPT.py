# - computeNormGrad (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, '{Density}=2*{CoordinateX}+{CoordinateX}*{CoordinateY}')
m = P.computeNormGrad(m, 'Density')
C.convertPyTree2File(m, 'out.cgns')
