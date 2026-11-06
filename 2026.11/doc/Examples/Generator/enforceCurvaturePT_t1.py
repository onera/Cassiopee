# - enforceCurvature (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Naca profile with lines
a = D.naca(12., 500)

# Distribution on the profile
Ni = 20; Nj = 20; Nk = 1; h = 1./(Ni-1)
b = G.cart((0,0,0), (h, 0.25/Nj,1), (Ni,Nj,Nk))
b = C.addBC2Zone(b, 'wall1','BCWall','jmin')
b = C.addBC2Zone(b, 'match1','BCMatch','imin',b,'imax',[1,2])
b = C.addBC2Zone(b, 'match2','BCMatch','imax',b,'imin',[1,2])
b = C.addBC2Zone(b, 'wall2','BCWall','jmax')
b = C.addVars(b,'Density'); b = C.addVars(b,'centers:cellN')
b = G.enforceCurvature(b, a, 0.6)
test.testT(b,1)
