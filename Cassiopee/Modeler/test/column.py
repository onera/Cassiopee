# - column (array) -
import Converter as C
import Modeler.Models as Models
import Transform as T

# Column with square foot
a1 = Models.column(R=1, N=10, h=10.)

# Column with round foot
a2 = Models.column2(R1=1, R2=0.8, N=20, h=10.)
a2 = T.translate(a2, (0,5,0))

# Column with
a3 = Models.column3(R1=1, R2=0.8, N=31, h=10.)
a3 = T.translate(a3, (0,10,0))

C.convertArrays2File([a1,a2,a3], 'out.plt')
