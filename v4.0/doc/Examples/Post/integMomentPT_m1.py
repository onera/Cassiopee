# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Converter.Mpi as Cmpi
import Post.Mpi as Pmpi
import KCore.test as test

# Test avec une grille cartesienne sur chaque proc
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'vx', Cmpi.rank)
m = C.initVars(m, 'vy', 2.)
m = C.initVars(m, 'vz', 1.)
m = C.initVars(m, 'centers:Fx', 2.)
m = C.initVars(m, 'centers:Fy', 2.)
m = C.initVars(m, 'centers:Fz', 3.)

#res1 = Pmpi.integMoment(m, center=(0,0,0), vector=['vx','vy','vz'])
#res2 = Pmpi.integMoment(m, center=(0,0,0), vector=['cneters:Fx','centers:Fy','centers:Fz'])
#if Cmpi.rank == 0: test.testO(res1+res2, 1)

# Test avec des []
t = C.newPyTree(['Base'])
if Cmpi.rank == 0:
    m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
    m = C.initVars(m, 'vx', Cmpi.rank)
    m = C.initVars(m, 'vy', 2.)
    m = C.initVars(m, 'vz', 1.)
    t[2][1][2] += [m]

res1 = Pmpi.integMoment(t, center=(0,0,0), vector=['vx','vy','vz'])
print(res1)
#if Cmpi.rank == 0: test.testO(res1, 2)
