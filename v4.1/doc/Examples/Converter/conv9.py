# - convertArrays2File -
import Converter as C
from numpy import *

a = zeros((2,2), float64)
# x et y du point 0
a[0,0] = 1
a[1,0] = 0.1
# x et y du point 1
a[0,1] = 2.
a[1,1] = -0.1
arrays = [['x,y', a, 10, 1, 1]]
C.convertArrays2File(arrays, "out.plt", "bin_tp")
