# - skySphere (array) -
import Converter as C
import Modeler.Models as Models

# Sky sphere
a = Models.skySphere((0,0,0), 1.)
C.convertArrays2File(a, 'out.plt')
