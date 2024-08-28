# - square -
import Converter as C
from Generator.Shapes import square

m = square((0,0), (1,1), 0.5, (10,10,5))
C.convertArrays2File(m, "out.plt", "bin_tp")
