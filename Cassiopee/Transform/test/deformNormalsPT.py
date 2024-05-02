# - deformNormals (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T

a = D.sphere6((0,0,0), 1., 10)
a = C.convertArray2Hexa(a)
a = T.join(a); a = G.close(a)

a = C.initVars(a, '{alpha}=0.5*{CoordinateX}')
a = T.deformNormals(a, 'alpha', niter=2)
C.convertPyTree2File(a, 'out.cgns')
