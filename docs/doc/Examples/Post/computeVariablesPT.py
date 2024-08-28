# - computeVariables (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
vars = ['Density','MomentumX', 'MomentumY', 'MomentumZ', \
        'EnergyStagnationDensity']
for v in vars: m = C.addVars(m, v)

# Pressure and Mach number extraction
m = P.computeVariables(m, ['Mach', 'Pressure'])
C.convertPyTree2File(m, 'out.cgns')
