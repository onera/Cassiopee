# - display1D (array) -
import Generator as G
import CPlot
import Converter as C
import numpy

# 1D data defined in arrays
b = G.cart((0,0,0), (0.1,1,1), (50,1,1))
c = G.cart((5,0,0), (0.1,1,1), (50,1,1))
B = [b, c]
CPlot.setState(gridSize=(1,2))
for i in range(100):
    B = C.initVars(B, '{f}=sin({x}+0.02*%d)'%i)
    B = C.initVars(B, '{g}={x}')
    CPlot.display1D(B, slot=0, bgBlend=1., gridPos=(0,0),
                    var1='x', var2='f')

# 1D data defined in numpys
x = numpy.linspace(0, 2*numpy.pi)
y = numpy.sin(x)
CPlot.display1D([x,y], slot=1, var1='x', var2='y', gridPos=(0,1), bgBlend=0.)
