# - blur (array) -
import Converter as C
import CPlot

a = C.convertFile2Arrays('one.png')
CPlot.blur(a, 0.8)
C.convertArrays2File(a, 'one2.png')
