# - enforceCurvature2 (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G

a = D.naca(12.)
ni = 101; distrib = G.cart((0,0,0),(1./(ni-1),1,1),(ni,1,1))
distrib2 = G.enforceCurvature2(distrib,a)
a = G.map(a, distrib2)
C.convertPyTree2File(a, "out.cgns")
