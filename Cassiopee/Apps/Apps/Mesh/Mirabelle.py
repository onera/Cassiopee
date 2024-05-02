# Remailleur structure ASM
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X
import Geom.PyTree as D
import numpy

#=====================================================================
# Helpers
#=====================================================================
# Retourne les raccords matchs d'une zone z
# Si donorName est specifie, retourne les raccords matchs correspondants
# uniquement a ce donneur
#=====================================================================
def getBCMatchs(z, donorName=None):
    bcs = []
    gc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')

    if donorName is None: # return all BCMatchs
        if gc is not None:
            bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        return bcs

    # Return BCMatchs with given donorName
    if gc is not None:
        l = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        for bc in l:
            oppBlock = Internal.getValue(bc)
            if oppBlock == donorName: bcs.append(bc)
    return bcs

#=====================================================================
# Retourne les datas d'un BCMatch : oppBlockName, range, donorRange, trf
#=====================================================================
def getBCMatchData(bc):
    oppBlock = Internal.getValue(bc)
    rnge = Internal.getNodeFromName1(bc, 'PointRange')
    rnge = Internal.range2Window(rnge[1])
    donor = Internal.getNodeFromName1(bc, 'PointRangeDonor')
    donor = Internal.range2Window(donor[1])
    trf = Internal.getNodeFromName1(bc, 'Transform')
    trf = Internal.getValue(trf)
    if trf.size == 2: # met toujours le transform en 3D
        trf2 = numpy.empty(3, dtype=numpy.int32)
        trf2[0:2] = trf[0:2]; trf2[2] = 3
        trf = trf2
    return (oppBlock, rnge, donor, trf)

#=====================================================================
# IN: rac = 1,2,3,4,5,6,-1,-2,-3,-4,-5,-6
# OUT: rac = 1,2,3
#=====================================================================
def getRac(rac):
    racl = abs(rac)
    if racl == 1: return 1
    elif racl == 2: return 1
    elif racl == 3: return 2
    elif racl == 4: return 2
    elif racl == 5: return 3
    else: return 3

#=====================================================================
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
# IN: window: [imin,imax,jmin,jmax,kmin,kmax]
# Retourne 1, 2, 3
#=====================================================================
def getWinDir2(win):
    if win[0] == win[1]: return 1
    if win[2] == win[3]: return 2
    if win[4] == win[5]: return 3
    return -1

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

#=====================================================================
# Casse les conditions aux limites pour que les raccords soient
# sur des faces pleines
#=====================================================================
def breakBCs(t):
    zones = Internal.getZones(t)
    stack = []; final = []
    for z in zones: stack.append(z)
    count = 0

    while len(stack) > 0:
        z = stack.pop(0)
        break__(z, stack, final, t)
        # DBG
        #count += 1
        #if count > 1: final += stack; break
        print(len(stack), len(final))

    t = C.newPyTree(['Base', final])
    return t

def break__(z, stack, final, t):
    allzones = stack+final
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Unstructured': final.append(z); return
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    bcs = getBCMatchs(z)
    for bc in bcs:
        (oppBlock, rnge, donor, trf) = getBCMatchData(bc)
        print(z[0], rnge, ni, nj, nk)
        dir = getWinDir2(rnge)

        # verifie si la fenetre est full
        if (dir == 2 or dir == 3) and rnge[0] > 1:
            z1 = T.subzone(z, (1,1,1), (rnge[0],nj,nk))
            z2 = T.subzone(z, (rnge[0],1,1), (ni,nj,nk))
            z1,z2 = T.split(z, dir, rnge[0], t)
            #adaptBCMatch(z, allzones, z1, z2, splitDir=1, splitIndex=rnge[0])
            stack.append(z1); stack.append(z2)
            return

        if (dir == 2 or dir == 3) and rnge[1] < ni:
            z1 = T.subzone(z, (1,1,1), (rnge[1],nj,nk))
            z2 = T.subzone(z, (rnge[1],1,1), (ni,nj,nk))
            #adaptBCMatch(z, allzones, z1, z2, splitDir=1, splitIndex=rnge[1])
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 3) and rnge[2] > 1: # need split, stack
            z1 = T.subzone(z, (1,1,1), (ni,rnge[2],nk))
            z2 = T.subzone(z, (1,rnge[2],1), (ni,nj,nk))
            #adaptBCMatch(z, allzones, z1, z2, splitDir=2, splitIndex=rnge[2])
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 3) and rnge[3] < nj: # need split, stack
            z1 = T.subzone(z, (1,1,1), (ni,rnge[3],nk))
            z2 = T.subzone(z, (1,rnge[3],1), (ni,nj,nk))
            #adaptBCMatch(z, allzones, z1, z2, splitDir=2, splitIndex=rnge[3])
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 2) and rnge[4] > 1: # need split, stack
            z1 = T.subzone(z, (1,1,1), (ni,nj,rnge[4]))
            z2 = T.subzone(z, (1,1,rnge[4]), (ni,nj,nk))
            #adaptBCMatch(z, allzones, z1, z2, splitDir=3, splitIndex=rnge[4])
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 2) and rnge[5] < nk: # need split, stack
            z1 = T.subzone(z, (1,1,1), (ni,nj,rnge[5]))
            z2 = T.subzone(z, (1,1,rnge[5]), (ni,nj,nk))
            #adaptBCMatch(z, allzones, z1, z2, splitDir=3, splitIndex=rnge[5])
            stack.append(z1); stack.append(z2)
            return

    final.append(z)

# adapte toutes les BCMatchs quand z est splitte en z1+z2
def adaptBCMatch(z, allzones, z1, z2, splitDir=1, splitIndex=1):
    print('Adapting BC of %s because of split'%z[0])
    dim = Internal.getZoneDim(z)
    ni = dim[1]; nj = dim[2]; nk = dim[3]

    # Suppression des raccords faisant intervenir z
    bcs = getBCMatchs(z)
    allOppBlocks = []
    for bc in bcs:
        oppBlock = Internal.getValue(bc)
        allOppBlocks.append(oppBlock)

    print(allOppBlocks)

    # Pour tous les blocs opposes au bloc splitte
    pool = []
    for a in allOppBlocks:
        zp = Internal.getNodeFromName(allzones, a)
        if zp is None:
            print('Warning: can not find %s in allzones'%a)
            bcs = []
        else:
            bcs = getBCMatchs(zp, z[0])
        for bc in bcs:
            (oppBlock, rnge, donor, trf) = getBCMatchData(bc)
            if oppBlock == z[0]: Internal._rmNode(zp, bc)
            if zp not in pool: pool.append(zp)

    pool += [z1,z2]
    print('pool')
    for i in pool: print(i[0])
    print('end pool')
    pool = X.connectMatch(pool, dim=dim[4])

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
        bcs = getBCMatchs(z)
        for bc in bcs:
            (oppBlock, rnge, donor, trf) = getBCMatchData(bc)
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
# zi : zone initiale
# zo : zone remappee
#==========================================================================
def _reportBCMatchs(zo, zi, t):
    dim = Internal.getZoneDim(zo)
    nio = dim[1]; njo = dim[2]; nko = dim[3]
    dim = Internal.getZoneDim(zi)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    gc = Internal.getNodeFromType1(zi, 'ZoneGridConnectivity_t')
    if gc is None: return
    bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
    for bc in bcs:
        (oppBlock, rnge, donor, trf) = getBCMatchData(bc)
        zopp = Internal.getNodeFromName2(t, oppBlock)
        dim = Internal.getZoneDim(zopp)
        niopp = dim[1]; njopp = dim[2]; nkopp = dim[3]
        donorDir = getWinDir(donor)
        if donorDir == 'imin':
            donorWin = [1, 1, 1, njo, 1, nko]
        elif donorDir == 'imax':
            donorWin = [niopp, niopp, 1, njo, 1, nko]
        elif donorDir == 'jmin':
            donorWin = [1, nio, 1, 1, 1, nko]
        elif donorDir == 'jmax':
            donorWin = [1, nio, njopp, njopp, 1, nko]
        elif donorDir == 'kmin':
            donorWin = [1, nio, 1, njo, 1, 1]
        elif donorDir == 'kmax':
            donorWin = [1, nio, 1, njo, nkopp, nkopp]

        print("Adding bcmatch", getWinDir(rnge), oppBlock, donorWin)

        C._addBC2Zone(zo, 'match1', 'BCMatch', getWinDir(rnge),
                      zoneDonor=oppBlock, rangeDonor=donorWin, trirac=trf)
    #Internal.printTree(zo)

# Adapt les BC range donor des Match referencant zo
def _adaptRanges(t, zo):
    dim = Internal.getZoneDim(zo)
    nio = dim[1]; njo = dim[2]; nko = dim[3]
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        gc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if gc is None: continue
        bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        for bc in bcs:
            (oppBlock, rnge, donor, trf) = getBCMatchData(bc)
            if oppBlock == zo[0]:
                #print("adapting range of %s"%z[0])
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
def _propagate(t, graph, stack, treated, linelet):

    # dir0 est la direction du traitement a appliquer
    # peut-etre : 1,2,3,4,5,6,-1,-2,-3,-4,-5,-6
    (B0, dir0) = stack.pop(0)

    # traite B0 suivant dir0
    zB0 = Internal.getNodeFromName2(t, B0)

    if dir0 < 0: l = T.reorder(linelet, (-1,2,3))
    else: l = linelet

    # traitement
    print(">> Je traite %s (dir=%d)\n"%(B0,dir0))
    dir0s = getRac(dir0)
    zp = G.map(zB0, linelet, dir=dir0s)
    # zp a le meme nom, il faut reporter les matchs de zB0 full face sur zp
    _reportBCMatchs(zp, zB0, t)
    # Adapte les donorRange de ceux qui referencent zp
    _adaptRanges(t, zp)

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

    if len(stack) > 0: _propagate(t, graph, stack, treated)


#=====================================================
#a = G.cart((0,0,0), (1,1,1), (10,10,10))
#b = G.cart((9,1,0), (1,1,1), (10,10,10))
#b = T.reorder(b, (2,1,3))
#t = C.newPyTree(['Base',a,b])
#t = X.connectMatch(t)
#zones = Internal.getZones(t)
#bcs = getBCMatchs(zones[0])
#(oppBlock, rnge, donor, trf) = getBCMatchData(bcs[0])

# suppose i,j,k index on donor zone
#print trfIndex((1,2,1), rnge, donor, trf)
#print reverseTrf(trf)
#trf = (1,3,2)
#print trfDir(2, trf)
#import sys; sys.exit()

# Construction de la distribution a remapper
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
    print(h1,h2,ind1,ind2)
    D.setH(l, ind1, h1)
    D.setH(l, ind2, h2)

    l = D.enforceh(l, N=N)
    C.convertPyTree2File(l, 'linelet.cgns')
    l = D.getCurvilinearAbscissa(l)
    C._initVars(l, '{CoordinateX}={s}')
    C._initVars(l, '{CoordinateY}=0.')
    C._initVars(l, '{CoordinateZ}=0.')
    #C._initVars(l,'{CoordinateX}')
    return l

#======================================================
# Main
#======================================================
if __name__ == "__main__":
    # Remesh case1.cgns - OK
    """
    FILE = 'case1.cgns'
    block = 'cart.0'
    dir = 2
    h1 = 0.001
    h2 = 1.
    N = 140
    """

    # Remesh case2.cgns - OK
    """
    FILE = 'case2.cgns'
    block = 'cart.0'
    dir = 2
    h1 = 0.001
    h2 = 1.
    N = 140
    """

    # Remesh naca2.cgns - OK

    FILE = 'naca.cgns'
    block = 'StructZone0.16'
    dir = 2
    h1 = 0.0001
    h2 = 2.
    N = 100


    # Remesh m6.cgns - OK
    """
    FILE = 'm6.cgns'
    block = 'm6wing.gfmt:3'
    dir = 2
    h1 = 0.004
    h2 = 3
    N = 40
    """

    # Remesh asf2.cgns - OK
    """
    FILE = 'asf2.cgns'
    block = 'blk-1-split-3.51'
    h1 = 0.02
    dir = 3
    h2 = 0.001
    N = 50
    """

    t = C.convertFile2PyTree(FILE)

    # Break les conditions aux limites matchs
    print("== Exploding...")
    t = breakBCs(t)
    C.convertPyTree2File(t, 'split.cgns')
    import sys; sys.exit()

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
    linelet = buildDistrib(t, block, dir, h1, h2, N)

    if dir == 2: dir = 3
    elif dir == 3: dir = 5
    stack = [(block, dir)]

    # Run
    print("== Running propagate...")
    treated = []
    _propagate(t, g, stack, treated)

    #Internal.printTree(t)
    C.convertPyTree2File(t, 'out.cgns')
