# - display1D (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import numpy

# 1D data defined in zones
b = G.cart((0,0,0), (0.1,1,1), (50,1,1))
c = G.cart((5,0,0), (0.1,1,1), (50,1,1))
B = [b, c]
CPlot.setState(gridSize=(1,2))
for i in range(100):
    C._initVars(B, '{f}=sin({CoordinateX}+0.01*%d)'%i)
    C._initVars(B, '{g}={CoordinateX}')
    CPlot.display1D(B, slot=0, bgBlend=1., gridPos=(0,0), var1='CoordinateX', var2='f')

# 1D data defined in numpys
import numpy
x = numpy.linspace(0, 2*numpy.pi)
y = numpy.sin(x)
CPlot.display1D([x,y], slot=1, var1='x', var2='y', gridPos=(0,1), bgBlend=0.8)
