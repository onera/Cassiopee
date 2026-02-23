# - wing (array) -
import Converter as C
import Modeler.Models as Models

a = Models.wing()
C.convertArrays2File(a, 'out.plt')
