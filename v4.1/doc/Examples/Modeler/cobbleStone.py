# - cobbleStone (array) -
import Converter as C
import Modeler.Models as Models

a = Models.cobbleStone(hx=0.5, hy=0.3, hz=0.2)
C.convertArrays2File(a, 'out.plt')
