# - getCurvilinearAbscissa (array) -
import Converter as C
import Geom as D
import Transform as T

a = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a = T.join (a, a2)
a3 = D.getCurvilinearAbscissa(a)
a = C.addVars([a, a3])
C.convertArrays2File(a, "out.plt")
