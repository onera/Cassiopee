# - magnitude (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G

a = D.sphere((0,0,0), 1., 50)
a = G.getNormalMap(a)
a = C.magnitude(a, ['centers:sx','centers:sy','centers:sz'])
C.convertPyTree2File(a, 'out.cgns')
