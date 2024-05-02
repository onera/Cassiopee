# - display (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import CPlot.PyTree as CPlot

a = D.sphere((0,0,0), 1.)
a = C.convertArray2Node(a)
C._initVars(a, '{radius}=0.01*{CoordinateX}')
CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[1,0])

CPlot.display(a, mode='render')
