# - bezier (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

# Bezier 1D
pts = D.polyline([(0.,0.,0.), (0.,1.,0.), (2.,1.,0.), (2.,0.,0.),
                  (4.,-1.,0.), (5.,6.,0.),])
a = D.bezier(pts, 100); a[0] = 'bezier'
C.convertPyTree2File(a, 'out.cgns')
