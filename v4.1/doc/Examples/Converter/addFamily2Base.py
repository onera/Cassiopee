# - addFamily2Base (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0, 360, 1, (50,20,20))
t = C.newPyTree(['Base', a])
# Add family name referencing a BCWall BC type
C._addFamily2Base(t[2][1], 'flap', 'BCWall')
# Add just a family name
C._addFamily2Base(t[2][1], 'component1')
C.convertPyTree2File(t, 'out.cgns')
