# - tunnel (array) -
import Modeler.Models as Models
import Geom as D
import Converter as C

l = D.line((0,0,0), (1,0,1), N=10)
a = Models.tunnel1(l)
C.convertArrays2File(a, 'out.plt')
