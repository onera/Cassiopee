# - computeVariables2 (array) -
import Converter as C
import Post as P
import Generator as G

# test sur un array
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))

Mach = 0.4
ROINF = 1.29
UINF = 263.55
AINF = UINF/Mach
ROUINF = ROINF*UINF
PINF = ROINF*AINF*AINF/1.4
PEINF = PINF/(0.4*ROINF)
ROEINF = ROINF*PEINF

c = C.array('ro,rou, rov,row,roE', ni, nj, 2)
c = C.initVars(c, 'ro', ROINF)
c = C.initVars(c, 'rou', ROUINF)
c = C.initVars(c, 'rov', 0.)
c = C.initVars(c, 'row', 0.)
c = C.initVars(c, 'roE', ROEINF)
m = C.addVars([m,c])
A = [m]

variables = ['ViscosityMolecular']
mu = P.computeVariables(m, variables, gamma=1.4, rgp=287.053, Cs=110.4, betas=0.000001458)
print(mu)
