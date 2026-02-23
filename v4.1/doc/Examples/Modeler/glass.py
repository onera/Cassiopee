# - glass (array) -

import Converter as C
import Modeler.Models as Models

a = Models.glass1()
C.convertArrays2File(a, 'out.plt')
