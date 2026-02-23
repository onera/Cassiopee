# - gencartmb (array) -
import Generator as G
import KCore.test as test

# 4 niveaux
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 1., (20,20,2))
h = 1.e-1; Dfar = 10.; nlvl = [5,5,5]
A = G.gencartmb([a], h, Dfar, nlvl)
test.testA(A, 1)
