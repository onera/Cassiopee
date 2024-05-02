# - addSeparationLine (array) -
import Geom as D
import Converter as C

# Add a line to a circle
a1 = D.circle((0,0,0), 1, 0., 360, 1000)
a2 = D.line((0.,1.,0.), (0.,2.,0), 100)
arrays = D.addSeparationLine(a1, a2)
C.convertArrays2File(arrays, "out.plt")
