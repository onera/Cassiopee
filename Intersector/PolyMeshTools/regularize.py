import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I
import Generator.PyTree as G
import time
import os, sys

def nb_cells(a):
  ncellsTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = I.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      ncellsTot += ncells
  return ncellsTot

def nb_faces(a):
  ncfacesTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
  	GEl = I.getElementNodes(z)
  	NGON = 0; found = False
  	for c in GEl:
  		if c[1][0] == 22: found = True; break
  		NGON += 1
  	if found:
  		node = GEl[NGON]
  		er = I.getNodeFromName1(node, 'ElementRange')
        ncfacesTot += er[1][1]-er[1][0]+1
  return ncfacesTot

def aglomerate(t, vr, vm):
	nb_cells0 = nb_cells(t)
	carry_on=1
	i=0
	while (carry_on == 1):  
		print "iter %s"%i
		t=XOR.agglomerateSmallCells(t, vmin=vm, vratio=vr)
		nb_cells1 = nb_cells(t)
		if (nb_cells1 == nb_cells0): carry_on=0
		if (carry_on == 0) : print "no cell found."
		if (nb_cells1 != nb_cells0) : print "%d cells have been aglomerated"%(nb_cells0-nb_cells1)
		nb_cells0 = nb_cells1
		i=i+1
		#if (i == 3) : carry_on=0
	return t

def aglomerateNonStar(t):
  nb_cells0 = nb_cells(t)
  carry_on=1
  i=0
  while (carry_on == 1):  
    print " "
    print "iter %s"%i
    t=XOR.agglomerateNonStarCells(t)
    print "check closure"
    XOR.checkCellsClosure(t)
    nb_cells1 = nb_cells(t)
    if (nb_cells1 == nb_cells0): carry_on=0
    if (carry_on == 0) : print "no cell found."
    if (nb_cells1 != nb_cells0) : 
      print "%d cells have been aglomerated"%(nb_cells0-nb_cells1)
      #C.convertPyTree2File(t, "nonstar_iter_%s.cgns"%i)
    nb_cells0 = nb_cells1
    i=i+1
    #if (i == 3) : carry_on=0
  return t


if len(sys.argv) is not 5 :
    print "ARG ERROR : 2 arguments to provide : IFILE OFILE VMIN VRATIO"
    sys.exit()

ifile=sys.argv[1]
ofile=sys.argv[2]
VMIN = float(sys.argv[3])
VRATIO = float(sys.argv[4])

t = C.convertFile2PyTree(ifile)

print " "
print " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
print " 1. agglomerate"
vratio = VRATIO*VRATIO*VRATIO
print " 1.1 : aglomerate cells vr %s"%(vratio)
t=aglomerate(t, vr=vratio, vm=VMIN)
vratio = VRATIO*VRATIO
print " 1.2 : aglomerate cells vr %s"%(vratio)
t=aglomerate(t, vr=vratio, vm=VMIN)
print " 1.3 : aglomerate cells vr %s"%(VRATIO)
t=aglomerate(t, vr=VRATIO, vm=VMIN)

C.convertPyTree2File(t, "agglo1.cgns")

for i in range(5):
  print " "
  print " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
  print " 2. split. iter %s"%(i)
  print " 2.1. prepareCellsSplit by convexifying any concave PH..."
  set = 1 # 0 for concave cells or 1 for non-centroid-star_shaped cells
  policy = 0 #0 : convexify concave pgs on PH set. 1 : starify concave pgs on PH set. 2 : starify any pgs at concave-chains ends
  t = XOR.prepareCellsSplit(t, PH_set = set, split_policy = policy, PH_conc_threshold = 0.01, PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8)
  print " 2.2 splitNonStarCells..."
  t = XOR.splitNonStarCells(t)
  print " 2.3 : simplify cells"
  t = XOR.simplifyCells(t, 0)# do not treat externals
  C.convertPyTree2File(t, "split%s.cgns"%(i))

print " "
print " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
print " 3. agglomerate"
vratio = VRATIO*VRATIO*VRATIO
print " 3.1 : aglomerate cells vr %s"%(vratio)
t=aglomerate(t, vr=vratio, vm=VMIN)
vratio = VRATIO*VRATIO
print " 3.2 : aglomerate cells vr %s"%(vratio)
t=aglomerate(t, vr=vratio, vm=VMIN)
print " 3.3 : aglomerate cells vr %s"%(VRATIO)
t=aglomerate(t, vr=VRATIO, vm=VMIN)

C.convertPyTree2File(t, "agglo2.cgns")

print " "
print " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
print " 4. agglomerate non-stars"
t=aglomerateNonStar(t)

print " 5. : simplify cells"
t = XOR.simplifyCells(t, 0)# do not treat externals

C.convertPyTree2File(t, "agglo3.cgns")

print " "
print " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
print " 6. collapse"
t = XOR.collapseUncomputableFaces(t)
t = G.close(t, -1.)

C.convertPyTree2File(t, "collapse.cgns")

t=aglomerate(t, vr=VRATIO, vm=VMIN)

C.convertPyTree2File(t, ofile)
