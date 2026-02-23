# - table (array) -
import Modeler.Models as Models
import Transform as T
import Converter as C

a = Models.table1(L=2., W=1., H=1.)
b = Models.table2(L=0.6, H=1., N=20); b = T.translate(b, (0,3,0))
C.convertArrays2File(a+b, 'out.plt')
