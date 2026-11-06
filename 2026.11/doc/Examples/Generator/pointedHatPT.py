# - pointedHat (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
c = D.circle( (0,0,0), 1., 360., 0., 100)
surf = G.pointedHat(c,(0.,0.,1.))
t = C.newPyTree(['Base']); t[2][1][2].append(surf)
C.convertPyTree2File(t, 'out.cgns')
