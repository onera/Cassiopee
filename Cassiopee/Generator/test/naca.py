# - naca (array) -
import Converter as C
import Generator as G
import Geom as D
import Transform as T

# Put a naca profile in msh
msh = D.naca(12., 5001); l1 = D.getLength(msh)

# Put a line in msh2
msh2 = D.line((1.,0.,0.),(20.,0.,0.), 5001); l2 = D.getLength(msh2)

# Join the two meshes
msh = T.join(msh, msh2)

# Add another line
msh2 = D.line((20.,0.,0.),(1.,0.,0.), 5001); l3 = D.getLength(msh2)
msh = T.join(msh2, msh)

C.convertArrays2File([msh], 'naca.plt')

# distribution
Ni = 200; Nj = 100
distrib = G.cart((0,0,0), (1./(Ni-1), 20./(Nj-1),1), (Ni,Nj,1))

distrib = G.enforcePlusY(distrib, 2.e-3, 50, 30)
distrib = G.enforceMoinsY(distrib, 2., 80, -60)
distrib = G.enforcePlusX(distrib, 1.e-4, 30, 50)
distrib = G.enforceMoinsX(distrib, 0.1, 80, -30)
distrib = G.enforceX(distrib, (l1*0.5)/(0.5*l1+l2), 1.25e-3, 30, 20)

distrib = G.enforcePoint(distrib, (l1*0.5)/(0.5*l1+l2))

distrib2 = T.symetrize(distrib, (0.,0.,0.), (0,1,0), (0,0,1))
distrib2 = T.reorder(distrib2, (-1,2,3))
distrib = T.join(distrib2, distrib)
distrib = T.contract(distrib, (0,0,0), (0,1,0), (0,0,1), 0.5)
distrib = T.translate(distrib, (0.5,0,0))
C.convertArrays2File([distrib], 'distrib.plt')
msh = T.reorder(msh, (1,2,3))

msh = G.hyper2D(msh, distrib, "C")

C.convertArrays2File([msh], 'out.plt')
