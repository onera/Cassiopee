# Remailleur structure ASM
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.Vector as Vector

#=====================================================================
# Helpers
#=====================================================================

#=====================================================================
# IN: rac = 1,2,3,4,5,6,-1,-2,-3,-4,-5,-6
# 1,2 : suivant imin, imax
# 3,4 : suivant jmin, jmax
# 5,6 : suivant kmin, kmax
# negatif, suivant -imin, ...
# OUT: rac = 1,2,3 (direction i,j ou k)
#=====================================================================
def getRac(rac):
    racl = abs(rac)
    if racl == 1: return 1
    elif racl == 2: return 1
    elif racl == 3: return 2
    elif racl == 4: return 2
    elif racl == 5: return 3
    else: return 3

#====================================================================
# Retourne la direction de la fenetre :
# "imin", "imax", "jmin", "jmax", "kmin", "kmax"
# La fenetre est supposee full
#=====================================================================
def getWinDir(win):
    if win[0] == win[1] and win[0] == 1: return 'imin'
    elif win[0] == win[1]: return 'imax'
    if win[2] == win[3] and win[2] == 1: return 'jmin'
    elif win[2] == win[3]: return 'jmax'
    if win[4] == win[5] and win[4] == 1: return 'kmin'
    elif win[4] == win[5]: return 'kmax'


#=====================================================================
# IN: oppRac de transform (1,2,3 ou -1,-2,-3)
# IN: window donor
# OUT: oppRac suivant minmax donor (1,2,3,4,5 ou -1,...)
#=====================================================================
def getOppRac(oppRac, donor):
    if oppRac == 1:
        if donor[0] == 1: oppRac = 1
        else: oppRac = 2
    elif oppRac == -1:
        if donor[0] == 1: oppRac = -1
        else: oppRac = -2
    elif oppRac == 2:
        if donor[2] == 1: oppRac = 3
        else: oppRac = 4
    elif oppRac == -2:
        if donor[2] == 1: oppRac = -3
        else: oppRac = -4
    elif oppRac == 3:
        if donor[4] == 1: oppRac = 5
        else: oppRac = 6
    elif oppRac == -3:
        if donor[4] == 1: oppRac = -5
        else: oppRac = -6
    return oppRac

# Donne les transformations par des raccords
# rnge : win raccord sur z1
# donor : win raccord sur z2
# trf : transform de z1 a z2
# input (i,j,k) dans z1
# Renvoie le (i,j,k) dans z2
def trfIndex(ind, rnge, donor, trf):
    (imin1, imax1, jmin1, jmax1, kmin1, kmax1) = rnge
    (imin2, imax2, jmin2, jmax2, kmin2, kmax2) = donor
    (trfi, trfj, trfk) = trf
    (i,j,k) = ind

    # indice i
    if trfi == 1: ip = i-imin1+imin2
    elif trfi == -1: ip = imax1-i+imin2
    elif trfj == 1: ip = j-jmin1+jmin2
    elif trfj == -1: ip = jmax1-j+jmin2
    elif trfk == 1: ip = k-kmin1+kmin2
    elif trfk == -1: ip = kmax1-k+kmin2

    # indice j
    if trfi == 2: jp = i-imin1+imin2
    elif trfi == -2: jp = imax1-i+imin2
    elif trfj == 2: jp = j-jmin1+jmin2
    elif trfj == -2: jp = jmax1-j+jmin2
    elif trfk == 2: jp = k-kmin1+kmin2
    elif trfk == -2: jp = kmax1-k+kmin2

    # indice k
    if trfi == 3: kp = i-imin1+imin2
    elif trfi == -3:kp = imax1-i+imin2
    elif trfj == 3: kp = j-jmin1+jmin2
    elif trfj == -3: kp = jmax1-j+jmin2
    elif trfk == 3: kp = k-kmin1+kmin2
    elif trfk == -3: kp = kmax1-k+kmin2

    return (ip,jp,kp)

def reverseTrf(trf):
    trfi2 = 1; trfj2 = 2; trfk2 = 3
    (trfi,trfj,trfk) = trf
    if trfi == 1: trfi2 = 1
    if trfi == -1: trfi2 = -1
    if trfj == 1: trfi2 = 2
    if trfj == -1: trfi2 = -2
    if trfk == 1: trfi2 = 3
    if trfk == -1: trfi2 = -3

    if trfi == 2: trfj2 = 1
    if trfi == -2: trfj2 = -1
    if trfj == 2: trfj2 = 2
    if trfj == -2: trfj2 = -2
    if trfk == 2: trfj2 = 3
    if trfk == -2: trfj2 = -3

    if trfi == 3: trfk2 = 1
    elif trfi == -3: trfk2 = -1
    elif trfj == 3: trfk2 = 2
    elif trfj == -3: trfk2 = -2
    elif trfk == 3: trfk2 = 3
    elif trfk == -3: trfk2 = -3

    return (trfi2,trfj2,trfk2)

def trfDir(dir, trf):
    return trf[dir-1]


#====================================================================
# B0, dir0s : bloc courant et direction du traitement
# dir1 : direction du raccord sur B0
#====================================================================
def stackIt(B0, dir0s, dir1, graph, treated, stack):
    if graph[B0][dir1] is not None:
        (B1, rac, trf) = graph[B0][dir1]
        dirt = trf[dir0s-1] # direction du traitement sur B1
        if dirt == 1: diro = 1
        elif dirt == 2: diro = 3
        elif dirt == 3: diro = 5
        elif dirt == -1: diro = -1
        elif dirt == -2: diro = -3
        elif dirt == -3: diro = -5
        if (B1,abs(dirt)) not in treated: stack.append((B1, diro))


#============================================================================
# Construit le graph de connection
# Pour chaque bloc : 1=imin, 2=imax, 3=jmin, 4=jmax, 5=kmin, 6=kmax
# renvoit le bloc oppose, le raccord correspondant sur le bloc oppose
# et le trf
#============================================================================
def getGraph(t):
    graph = {}
    zones = Internal.getZones(t)
    for z in zones:
        neigh = [None, None, None, None, None, None, None]
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Unstructured': continue
        bcs = T.getBCMatchs__(z)
        for bc in bcs:
            (oppBlock, rnge, donor, trf, periodic) = T.getBCMatchData__(bc)
            racc = 0
            if rnge[0] == rnge[1] and rnge[0] == 1:
                racc = 1
                oppRacc = getOppRac(trf[0], donor)
            elif rnge[0] == rnge[1] and rnge[0] == dim[1]:
                racc = 2
                oppRacc = getOppRac(trf[0], donor)
            elif rnge[2] == rnge[3] and rnge[2] == 1:
                racc = 3
                oppRacc = getOppRac(trf[1], donor)
            elif rnge[2] == rnge[3] and rnge[2] == dim[2]:
                racc = 4
                oppRacc = getOppRac(trf[1], donor)
            elif rnge[4] == rnge[5] and rnge[4] == 1:
                racc = 5
                oppRacc = getOppRac(trf[2], donor)
            elif rnge[4] == rnge[5] and rnge[4] == dim[3]:
                racc = 6
                oppRacc = getOppRac(trf[2], donor)

            neigh[racc] = (oppBlock, oppRacc, trf)
        graph[z[0]] = neigh
    return graph

#==========================================================================
# recopie les BC Matchs de zi sur zo
# Les BCMatchs sont considerees full
# Seul les point range sont update, les pointRangeDonor sont faits a la fin
# zi : zone initiale
# zo : zone remappee
#==========================================================================
def _reportBCMatchs(zo, zi, t, dir0):
    dim = Internal.getZoneDim(zo)
    nio = dim[1]; njo = dim[2]; nko = dim[3]
    dim = Internal.getZoneDim(zi)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    gc = Internal.getNodeFromType1(zi, 'ZoneGridConnectivity_t')
    if gc is None: return
    bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
    for bc in bcs:
        (oppBlock, rnge, donor, trf, periodic) = T.getBCMatchData__(bc)
        donorWin = donor
        print("Adding bcmatch", getWinDir(rnge), oppBlock, donorWin)

        C._addBC2Zone(zo, 'match1', 'BCMatch', getWinDir(rnge),
                      zoneDonor=oppBlock, rangeDonor=donorWin, trirac=trf)
    #Internal.printTree(zo)


# The right one
def _adaptDonorRanges(t):
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        gc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if gc is None: continue
        bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        for bc in bcs:
            (oppBlock, rnge, donor, trf, periodic) = T.getBCMatchData__(bc)
            zopp = Internal.getNodeFromName2(t, oppBlock)
            dimOpp = Internal.getZoneDim(zopp)
            nio = dimOpp[1]; njo = dimOpp[2]; nko = dimOpp[3]
            donorDir = getWinDir(donor)
            if donorDir == 'imin':
                donorWin = [1, 1, 1, njo, 1, nko]
            elif donorDir == 'imax':
                donorWin = [nio, nio, 1, njo, 1, nko]
            elif donorDir == 'jmin':
                donorWin = [1, nio, 1, 1, 1, nko]
            elif donorDir == 'jmax':
                donorWin = [1, nio, njo, njo, 1, nko]
            elif donorDir == 'kmin':
                donorWin = [1, nio, 1, njo, 1, 1]
            elif donorDir == 'kmax':
                donorWin = [1, nio, 1, njo, nko, nko]
            pd = Internal.getNodeFromName1(bc, 'PointRangeDonor')
            donorWin = Internal.window2Range(donorWin)
            Internal.setValue(pd, donorWin)

#======================================================================
# Propage un traitement suivant le graph
# B0: bloc a remailler
# dir0: direction du remaillage sur B0 (1,2,3,4,5,6,-1,-2,-3,-4,-5,-6)
#======================================================================
def _propagate(t, graph, stack, treated, linelet, h1=None, h2=None, isAvg=False, nAvg=2):
    # dir0 est la direction du traitement a appliquer
    # peut-etre : 1,2,3,4,5,6,-1,-2,-3,-4,-5,-6
    (B0, dir0) = stack.pop(0)

    # traite B0 suivant dir0
    zB0 = Internal.getNodeFromName2(t, B0)

    if dir0 < 0:
        llt = T.reorder(linelet, (-1,2,3))
        C._initVars(llt, "{CoordinateX}=1.-{CoordinateX}")
    else:
        llt = linelet

    # traitement
    dir0s = getRac(dir0)
    print(">> Je traite %s (dir=%d, extDir=%d)\n"%(B0,dir0s,dir0))

    zp = G.map(zB0, llt, dir=dir0s,h1=h1,h2=h2,isAvg=isAvg,nAvg=nAvg)
    # zp a le meme nom, il faut reporter les matchs de zB0 full face sur zp
    _reportBCMatchs(zp, zB0, t, dir0s)
    # Adapte les donorRanges de ceux qui referencent zp
    #_adaptDonorRanges(t, zp)

    (p, c) = Internal.getParentOfNode(t, zB0)
    p[2][c] = zp

    dir0 = abs(dir0)
    treated.append((B0, dir0s))

    # propage dans les deux directions complementaires
    if dir0s == 1:
        stackIt(B0, dir0s, 3, graph, treated, stack)
        stackIt(B0, dir0s, 4, graph, treated, stack)
        stackIt(B0, dir0s, 5, graph, treated, stack)
        stackIt(B0, dir0s, 6, graph, treated, stack)

    elif dir0s == 2:
        stackIt(B0, dir0s, 1, graph, treated, stack)
        stackIt(B0, dir0s, 2, graph, treated, stack)
        stackIt(B0, dir0s, 5, graph, treated, stack)
        stackIt(B0, dir0s, 6, graph, treated, stack)

    elif dir0s == 3:
        stackIt(B0, dir0s, 1, graph, treated, stack)
        stackIt(B0, dir0s, 2, graph, treated, stack)
        stackIt(B0, dir0s, 3, graph, treated, stack)
        stackIt(B0, dir0s, 4, graph, treated, stack)

    if len(stack) > 0: _propagate(t, graph, stack, treated, linelet, h1=h1, h2=h2, isAvg=isAvg, nAvg=nAvg)

# Construction de la distribution a remapper a partir de dir, h1, h2 et N
# si h1,h2 < 0: conserve le h1,h2 original
def buildDistrib(t, block, dir, h1, h2, N):
    z = Internal.getNodeFromName2(t, block)
    dim = Internal.getZoneDim(z)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    if dir == 1:
        l = T.subzone(z, (1,1,1), (ni,1,1))
        ind1 = 0; ind2 = ni-1
    elif dir == 2:
        l = T.subzone(z, (1,1,1), (1,nj,1))
        l = T.reorder(l, (3,1,2))
        ind1 = 0; ind2 = nj-1
    else:
        l = T.subzone(z, (1,1,1), (1,1,nk))
        l = T.reorder(l, (2,3,1))
        ind1 = 0; ind2 = nk-1
    # Enforce
    P0 = C.getValue(l, 'GridCoordinates', ind1)
    P1 = C.getValue(l, 'GridCoordinates', ind1+1)
    h1s = Vector.norm(Vector.sub(P1, P0))
    P0 = C.getValue(l, 'GridCoordinates', ind2-1)
    P1 = C.getValue(l, 'GridCoordinates', ind2)
    h2s = Vector.norm(Vector.sub(P1, P0))
    if h1 > 0: D.setH(l, ind1, h1)
    else: D.setH(l, ind1, h1s)
    if h2 > 0: D.setH(l, ind2, h2)
    else: D.setH(l, ind2, h2s)
    l = D.enforceh(l, N=N)
    #C.convertPyTree2File(l, 'linelet.cgns')
    # Converti en abscisse curviligne
    l = D.getCurvilinearAbscissa(l)
    C._initVars(l, '{CoordinateX}={s}')
    C._initVars(l, '{CoordinateY}=0.')
    C._initVars(l, '{CoordinateZ}=0.')
    return l

# Construction de la distribution a remapper pour un maillage loi
# de paroi
# si h1 est negatif, wall law imin
# si h2 est negatif, wall law imax
def buildDistribWallLaw(t, block, dir, h1, h2, N):
    # subzone
    z = Internal.getNodeFromName2(t, block)
    dim = Internal.getZoneDim(z)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    if dir == 1:
        l = T.subzone(z, (1,1,1), (ni,1,1))
        ind1 = 0; ind2 = ni-1
    elif dir == 2:
        l = T.subzone(z, (1,1,1), (1,nj,1))
        l = T.reorder(l, (3,1,2))
        ind1 = 0; ind2 = nj-1
    else:
        l = T.subzone(z, (1,1,1), (1,1,nk))
        l = T.reorder(l, (2,3,1))
        ind1 = 0; ind2 = nk-1
    # subzone again
    if h1 < 0:
        l1 = T.subzone(l, (1,1,1), (N,1,1))
        l2 = T.subzone(l, (N,1,1), (-1,1,1))
        P0 = C.getValue(l1, 'GridCoordinates', 0)
        P1 = C.getValue(l1, 'GridCoordinates', N-1)
        ll = D.line(P0, P1, N=2) 
        l = T.join(ll, l2)
    else:
        npts = C.getNPts(l)
        l1 = T.subzone(l, (npts-N+1,1,1), (npts,1,1))
        l2 = T.subzone(l, (1,1,1), (npts-N+1,1,1))
        P0 = C.getValue(l1, 'GridCoordinates', 0)
        P1 = C.getValue(l1, 'GridCoordinates', N-1)
        ll = D.line(P0, P1, N=2) 
        l = T.join(l2, ll)
    C.convertPyTree2File(l, 'linelet.cgns')

    # Converti en abscisse curviligne
    l = D.getCurvilinearAbscissa(l)
    C._initVars(l, '{CoordinateX}={s}')
    C._initVars(l, '{CoordinateY}=0.')
    C._initVars(l, '{CoordinateZ}=0.')
    return l

# Return the first wall in t
# Return zoneName, window and window direction
def getAWallWindow(t):
    for z in Internal.getZones(t):
        BCs = Internal.getNodesFromName(z, 'BC_t')
        for b in BCs:
            if Internal.getValue(b) == 'BCWall':
                pr = Internal.getNodeFromName(b, 'PointRange')
                w = Internal.range2Window(pr[1])
                dir = T.getWinDir(w)
                return (z[0], w, dir)
    return None

#======================================================
# Main
#======================================================
if __name__ == "__main__":


    # FILE: nom du fichier contenant le maillage structure avec les raccords
    # block: block de depart sur lequel on definit la linelet
    # dir : direction for remeshing
    # 1: i = cte
    # 2: j = cte
    # 3: k = cte
    # h1: step for imin, si h1 negatif, loi de paroi sur N mailles en imin
    # h2: step for imax, si h2 negatif, loi de paroi sur N mailles en imax
    # N: number of points after remeshing

    CASE = 'asf2WL'

    # Remesh case1.cgns - OK
    if CASE == 'case1BL':
        FILE = 'case1.cgns'
        block = 'cart.0'
        dir = 2
        h1 = 0.001
        h2 = 1.
        N = 140

    # Remesh case2.cgns - OK
    elif CASE == 'case2BL':
        FILE = 'case2.cgns'
        block = 'cart.0'
        dir = 2
        h1 = 0.001
        h2 = 1.
        N = 140

    # Remesh case3.cgns - OK
    elif CASE == 'case3BL':
        FILE = 'case3.cgns'
        block = 'cart.1'
        dir = 2
        h1 = 0.001
        h2 = 1.
        N = 140

    # Remesh case4.cgns - OK
    elif CASE == 'case4BL':
        FILE = 'case4.cgns'
        block = 'StructZone0.32'
        dir = 2
        h1 = 0.001
        h2 = 0.0287
        N = 100

    # Remesh case7.cgns - OK
    elif CASE == 'case7BL':
        FILE = 'case7.cgns'
        block = 'cart.10'
        dir = 2
        h1 = 0.01
        h2 = 1.
        N = 60

    # Remesh case9.cgns - OK
    elif CASE == 'case9BL':
        FILE = 'case9.cgns'
        block = 'cart'
        dir = 2
        h1 = 0.01
        h2 = 1.
        N = 60

    # Remesh case10.cgns - OK
    elif CASE == 'case10BL':
        FILE = 'case10.cgns'
        block = 'cart'
        dir = 2
        h1 = 0.01
        h2 = 1.
        N = 60

    # Remesh case11.cgns - OK
    elif CASE == 'case11BL':
        FILE = 'case11.cgns'
        block = 'cart'
        dir = 2
        h1 = 0.01
        h2 = 1.
        N = 60

    # Remesh case12.cgns - OK
    elif CASE == 'case12BL':
        FILE = 'case12.cgns'
        block = 'cart'
        dir = 1
        h1 = 0.01
        h2 = 1.
        N = 60

    # Remesh naca.cgns - OK
    if CASE == 'nacaBL':
        FILE = 'naca.cgns'
        block = 'StructZone0.4'
        dir = 2
        h1 = 0.0001
        h2 = 2.
        N = 100

    # Remesh m6.cgns - OK
    elif CASE == 'm6BL':
        FILE = 'm6.cgns'
        block = 'm6wing.gfmt:3'
        dir = 2
        h1 = 0.004
        h2 = 3
        N = 40

    # Remesh asf2.cgns - KO
    elif CASE == 'asf2BL':
        FILE = 'asf2.cgns'
        block = 'blk-1-split-1.5'
        dir = 3
        h1 = 0.02
        h2 = 0.001
        N = 50

    # Remesh cylinder - OK
    elif CASE == 'cylBL':
        FILE = 'cyl.cgns'
        block = 'cylinder'
        dir = 2
        h1 = 0.01
        h2 = 0.01
        N = 32

    # Remesh cylinder2 - OK
    elif CASE == 'cyl2BL':
        FILE = 'cyl2.cgns'
        block = 'cylinder.1'
        dir = 2
        h1 = 0.01
        h2 = 0.1
        N = 32

    # Remesh cylinder2 with wall law - OK
    elif CASE == 'cyl2WL':
        FILE = 'cyl2.cgns'
        block = 'cylinder.1'
        dir = 2
        h1 = -0.01
        h2 = 0.01
        N = 5

    # Remesh asf2.cgns wall law - KO
    elif CASE == 'asf2WL':
        FILE = 'asf2.cgns'
        block = 'blk-1-split-1.5'
        dir = 3
        h1 = 0.02
        h2 = -0.001
        N = 10

    # Remesh dauphin - OK
    elif CASE == 'dauphinBL':
        FILE = 'dauphin.cgns'
        block = 'StructZone0.14'
        dir = 2
        h1 = 0.0001
        h2 = 1.
        N = 100

    # Remesh dauphin WL - OK
    elif CASE == 'dauphinWL':
        FILE = 'dauphin.cgns'
        block = 'StructZone0.14'
        dir = 2
        h1 = -0.1
        h2 = 1.
        N = 13

    t = C.convertFile2PyTree(FILE)

    # Break les conditions aux limites matchs
    print("== Exploding...")
    T._splitFullMatch(t)
    C.convertPyTree2File(t, 'split.cgns')

    # Construit le graph des raccords match
    print("== Computing match graph...")
    g = getGraph(t)
    print(g)

    # linelet reguliere
    #Npts = 10
    #linelet = G.cart((0,0,0), (1./(Npts-1),1,1), (Npts,1,1))

    # linelet irreguliere
    #Npts = 5
    #linelet1 = G.cart((0.3,0,0), (0.7/(Npts-1),1,1), (Npts,1,1))
    #linelet2 = G.cart((0,0,0), (0.3,1,1), (2,1,1))
    #linelet = T.join(linelet2, linelet1)
    #C.convertPyTree2File(linelet, 'line.plt')

    # linelet avec resserement
    #linelet = G.enforceMoinsX(linelet, 0.1, (5,5))

    # specification du resserement (dir=1,2,3)

    # Determine la distribution
    print("== Computing linelet...")

    if h1 < 0 or h2 < 0: linelet = buildDistribWallLaw(t, block, dir, h1, h2, N)
    else: linelet = buildDistrib(t, block, dir, h1, h2, N)

    if dir == 2: dir = 3
    elif dir == 3: dir = 5
    stack = [(block, dir)]

    # Run
    print("== Running propagate...")
    treated = []
    _propagate(t, g, stack, treated, linelet)

    # Adapte les donneurs a la fin
    _adaptDonorRanges(t)

    #Internal.printTree(t)
    C.convertPyTree2File(t, 'out.cgns')
