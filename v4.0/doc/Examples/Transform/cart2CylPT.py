# - cart2Cyl (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360., 1., (360,20,10))
T._cart2Cyl(a, (0.,0.,0.),(0,0,1))
C.convertPyTree2File(a, 'out.cgns')
