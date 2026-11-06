# - scanner (pyTree) -
import OCC.Scanner as Scanner
import OCC.PyTree as OCC
import Geom.PyTree as D

# model
a = D.sphere( (0,0,0), 1., N=100 )
surf = Scanner.scan(a, (0,0,1), 10)
OCC.writeCAD(surf, 'out.step')