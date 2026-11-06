# - computeDiv2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

m = G.cartNGon((0,0,0), (1,1,1), (4,4,4))
C._initVars(m, '{centers:fldX}= cos({centers:CoordinateX})')
C._initVars(m, '{centers:fldY}= 4.*{centers:CoordinateY}')
C._initVars(m, '{centers:fldZ}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')
m = P.computeDiv2(m, 'centers:fld')
C.convertPyTree2File(m, 'out.cgns')
