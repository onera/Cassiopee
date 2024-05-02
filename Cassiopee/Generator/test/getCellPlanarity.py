# - getCellPlanarity (array) -
import Converter as C
import Generator as G
import Geom as D

a = D.sphere( (0,0,0), 1., 10)
p = G.getCellPlanarity(a)
p = C.center2Node(p); a = C.addVars([a, p])
C.convertArrays2File([a], 'out.plt')
