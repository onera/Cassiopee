# - hub (array) -
import Converter as C
import Transform as T
import Modeler.Models as Models

# Hub1
a = Models.hub1(1.,1.,0.5)

C.convertArrays2File(a, 'out.plt')
