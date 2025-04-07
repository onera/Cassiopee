# - extractPlane (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G

m = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (50,50,50))
m = T.rotate(m, (0,0,0), (1,0,0), 35.)
m = C.initVars(m,'Density',1); m = C.initVars(m,'centers:cellN',1)
z = P.extractPlane(m, (0.5, 1., 0., 1), 2)
C.convertPyTree2File(z, 'out.cgns')
