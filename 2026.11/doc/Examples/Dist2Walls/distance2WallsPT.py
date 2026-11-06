# - distance2Walls (pyTree) -
import Dist2Walls.PyTree as Dist2Walls
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
sphere = D.sphere((1.2,0.,0.),0.2,100)
t = C.newPyTree(['Base',a])
t = Dist2Walls.distance2Walls(t, sphere)
C.convertPyTree2File(t, 'out.cgns')
