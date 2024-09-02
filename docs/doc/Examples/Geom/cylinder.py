# - cylinder (array) -
import Geom as D
import Converter as C

a = D.cylinder((0,0,0), 1., 10.)
C.convertArrays2File(a, 'out.plt')
