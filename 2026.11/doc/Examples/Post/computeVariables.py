# - computeVariables (array) -
import Converter as C
import Post as P
import Generator as G

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
c = C.array('ro,rou, rov,row,roE', ni, nj, 2)
c = C.initVars(c, 'ro', 1.)
c = C.initVars(c, 'rou', 1.)
c = C.initVars(c, 'rov', 0.)
c = C.initVars(c, 'row', 0.)
c = C.initVars(c, 'roE', 1.)
m = C.addVars([m, c])

#--------------------------------------------
# Pressure and Mach number extraction
# default values of rgp and gamma are used
#--------------------------------------------
p = P.computeVariables(m, ['Mach', 'Pressure'])
m = C.addVars([m, p])
C.convertArrays2File(m, "out.plt")
