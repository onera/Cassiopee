# - panel (array) -
import Converter as C
import Modeler.Models as Models

# Panneau avec un pied
a = Models.panel("Bonjour", h=5)

# Frame
b = Models.frame(5, 4, 0.1)

C.convertArrays2File(b, 'out.plt')
