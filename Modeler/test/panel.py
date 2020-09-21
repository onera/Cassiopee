# - panel (array) -
import Converter as C
import Modeler.Models as Models

a = Models.panel("Bonjour", h=5)
C.convertArrays2File(a, 'out.plt')
