# - axisym (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

# Axisym a curve
a0 = D.line((0.5,0,0), (0.6,0,1))
a = D.axisym(a0, (0.,0.,0.), (0.,0.,1.), 360., 360)
C.convertPyTree2File(a, "out.cgns")

# Axisym a curve with varying r
a0 = D.line((1.0,0,0), (0.,0,1))
a1 = D.circle((0,0,0), 2.)
a = D.axisym(a0, (0.,0.,0.), (0.,0.,1.), rmod=a1)
C.convertPyTree2File([a,a0,a1], "out1.cgns")

# Axisym a 2D cart grid
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = D.axisym(a, (1.,0.,0.), (0.,1.,0.), 30., 4)
C.convertPyTree2File(a, 'out2.cgns')
