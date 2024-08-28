# - gencartmb (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

#-------------------------------------------------------
# donnees utilisateurs
#-------------------------------------------------------
# -- maillages de corps --
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (20,20,10))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')

# In : pas de la grille cartesienne la plus fine
h = 1.e-1
# In : distance d'eloignement au corps
Dfar = 10.

# In : nb de points par niveau :
# ici 4 niveaux, mais le dernier est calcule automatiquement
nlvl = [5,5,5] # nlvl[0] : grille grossiere

# retourne la liste des zones correspondant aux grilles cartesiennes
t = C.newPyTree(['Bodies','Cart']); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.)
t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
zones = G.gencartmb(t[2][1], h, Dfar, nlvl)
t[2][2][2] += zones
test.testT(t,1)
