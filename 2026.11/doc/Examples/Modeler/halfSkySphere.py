# - halfSkySphere (array) -
import Converter as C
import Modeler.Models as Models

# Half sky sphere
a = Models.halfSkySphere((0,0,0), 1.)
C.convertArrays2File(a, 'out.plt')
