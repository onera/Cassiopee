# - display (array) -
# Pour les NGons
import CPlot
import Generator as G
import Converter as C

# 3D
a = G.cartNGon((0,0,0), (1,1,1), (20,20,20))
a = C.initVars(a, '{f}={x}')

# 2D
b = G.cartNGon((30,0,0), (1,1,1), (20,20,1))
b = C.initVars(b, '{f}={x}')

# 1D
c = G.cartNGon((0,30,0), (1,1,1), (20,1,1))
c = C.initVars(c, '{f}={x}')

CPlot.display([a,b,c], mode=0, meshStyle=2)
