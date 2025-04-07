# - gapfixer (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

# Fix the gap inside a circle drawn on a plane
a = D.circle((0,0,0), 1, N=100)
a = C.convertArray2Tetra(a); a = G.close(a)
b = G.cart((-2.,-2.,0.), (0.1,0.1,1.), (50,50,1))
a1 = G.gapfixer(a, b)
C.convertPyTree2File(a1, 'out.cgns')

# Fill the gap in the circle, using one defined point
hp = D.point((0.5, 0.5, 0.))
a2 = G.gapfixer(a, b, hp, refine=0)
C.convertPyTree2File(a2, 'outHP.cgns')
