# - stitchedHat (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
c = D.circle((0,0,0), 1., 360., 0., 100)
c = T.contract(c, (0,0,0), (0,1,0), (0,0,1), 0.1)
c = G.stitchedHat(c, (0,0,0), 1.e-4)
C.convertPyTree2File(c, 'out.cgns')
