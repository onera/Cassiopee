# - computeGrad2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

m = G.cartNGon((0,0,0), (1,1,1), (4,4,4))
C._initVars(m, '{centers:Density}= 2*{centers:CoordinateX}+{centers:CoordinateY}*{centers:CoordinateZ}')
m = P.computeGrad2(m, 'centers:Density')
C.convertPyTree2File(m, 'out.cgns')
