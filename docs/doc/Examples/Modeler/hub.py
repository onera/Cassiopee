# - hub (array) -
import Converter as C
import Modeler.Models as Models

h1 = Models.hub1(0.1, -2.e-3, 2.e-3, 5.e-3, 51)

h2 = Models.hub1(0.1, -2.e-3, 2.e-3, 7.e-3, 71)

h3 = Models.hub1(0.1, -2.e-3, 12.e-3, 30.e-3, 151)

C.convertArrays2File(h3, 'out.plt')
