# - axisym (array) -
import Generator as G
import Converter as C
import Geom as D

# Axisym a curve
a0 = D.line((0.5,0,0), (0.6,0,1))
a = D.axisym(a0, (0.,0.,0.), (0.,0.,1.), 360., 360)
C.convertArrays2File(a, "out.plt")

# Axisym a curve with varying r
a0 = D.line((1.0,0,0), (0.,0,1))
a1 = D.circle((0,0,0), 2.)
import Modeler.Models as Models
a1 = Models.circle2(1, 0.8)
a = D.axisym(a0, (0.,0.,0.), (0.,0.,1.), rmod=a1)
C.convertArrays2File([a,a0,a1], "out1.plt")

# Axisym a 2D cart grid
a0 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = D.axisym(a0, (1.,0.,0.), (0.,1.,0.), 30., 4)
C.convertArrays2File(a, "out2.plt")
