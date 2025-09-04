# - tkBlader -
"""App to create blade meshes."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Evalue la qualite du maillage m
# Plus le score est elevee pour le maillage est mauvais
#==============================================================================
def quality(meshes):
    import Converter as C
    import Generator as G
    score = 0.
    for m in meshes:
        ortho = G.getOrthogonalityMap(m)
        vol = G.getVolumeMap(m)
        min1 = C.getMinValue(ortho, 'orthogonality')
        max1 = C.getMaxValue(ortho, 'orthogonality')
        min2 = C.getMinValue(vol, 'vol')
        max2 = C.getMaxValue(vol, 'vol')
        score = max(score, min1); score = max(score, max1)
        if min2 < 1.e-12 and max2 > 0: score += 1000.
        elif max2 < 1.e-12 and min2 > 0: score += 1000.
    return score

#==============================================================================
# IN: a1,a2,a3: les 3 cotes du triangle
#==============================================================================
def trimesh(a1, a2, a3):
    import Converter as C
    import Generator as G
    import Transform as T
    import Geom as D
    N1 = a1[2]; N2 = a2[2]; N3 = a3[2]

    # Assure la continuite
    P0 = (a1[1][0,N1-1], a1[1][1,N1-1], a1[1][2,N1-1])
    P00 = (a2[1][0,0], a2[1][1,0], a2[1][2,0])
    P01 = (a2[1][0,N2-1], a2[1][1,N2-1], a2[1][2,N2-1])
    if (abs(P0[0]-P00[0]) + abs(P0[1]-P00[1]) + abs(P0[2]-P00[2]) > 1.e-6
            and abs(P0[0]-P01[0]) + abs(P0[1]-P01[1]) + abs(P0[2]-P01[2]) > 1.e-6):
        t = a2; a2 = a3; a3 = t
        N2 = a2[2]; N3 = a3[2]
        P00 = (a2[1][0,0], a2[1][1,0], a2[1][2,0])

    if abs(P0[0]-P00[0]) > 1.e-6: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[1]-P00[1]) > 1.e-6: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[2]-P00[2]) > 1.e-6: a2 = T.reorder(a2, (-1,2,3))
    P0 = (a2[1][0,N2-1], a2[1][1,N2-1], a2[1][2,N2-1])
    P00 = (a3[1][0,0], a3[1][1,0], a3[1][2,0])
    if abs(P0[0]-P00[0]) > 1.e-6: a3 = T.reorder(a3, (-1,2,3))
    elif abs(P0[1]-P00[1]) > 1.e-6: a3 = T.reorder(a3, (-1,2,3))
    elif abs(P0[2]-P00[2]) > 1.e-6: a3 = T.reorder(a3, (-1,2,3))
    #C.convertArrays2File([a1,a2,a3], 'order.plt')

    # Verif de N
    Nt = N3-N2+N1+1
    if Nt/2-Nt*0.5 != 0: return [0, 'N3-N2+N1 must be odd.', 0]
    N = Nt/2
    if N < 2: return [0, 'invalid number of points for this operation.', 0]
    if N > N1-1: return [0, 'invalid number of points for this operation.',0]
    if N > N2-1: return [0, 'invalid number of points for this operation.',0]

    # Center
    w1 = C.array('weight', N1, 1, 1)
    w1 = C.initVars(w1, 'weight', 1)
    w2 = C.array('weight', N2, 1, 1)
    w2 = C.initVars(w2, 'weight', 1)
    w3 = C.array('weight', N3, 1, 1)
    w3 = C.initVars(w3, 'weight', 1)
    CC = G.barycenter([a1,a2,a3], [w1,w2,w3])

    # Subzones
    s1 = T.subzone(a1, (1,1,1), (N1-N+1,1,1))
    s2 = T.subzone(a1, (N1-N+1,1,1), (N1,1,1))
    s3 = T.subzone(a2, (1,1,1), (N2-N1+N,1,1))
    s4 = T.subzone(a2, (N2-N1+N,1,1), (N2,1,1))
    s5 = T.subzone(a3, (1,1,1), (N,1,1))
    s6 = T.subzone(a3, (N,1,1), (N3,1,1))

    # Lines
    index = N1-N
    P01 = (a1[1][0,index], a1[1][1,index], a1[1][2,index])
    index = N2-N1+N-1
    P12 = (a2[1][0,index], a2[1][1,index], a2[1][2,index])
    index = N-1
    P23 = (a3[1][0,index], a3[1][1,index], a3[1][2,index])
    l1 = D.line(CC, P01, N=N2-N1+N)
    l2 = D.line(CC, P12, N=N)
    l3 = D.line(CC, P23, N=N1-N+1)

    #C.convertArrays2File([l1,l2,l3,s1,s2,s3,s4,s5,s6], 'borders.plt')

    # TFIs
    m1 = G.TFI([s1,l1,l3,s6])
    m2 = G.TFI([l1,s2,s3,l2])
    m3 = G.TFI([l2,s4,s5,l3])

    #return [l1,l2,l3,s1,s2,s3,s4,s5,s6]
    #return [s1,l1,l3,s6,m1]
    return [m1,m2,m3]

#==============================================================================
def mono2mesh(a1, a2):
    import Transform as T
    import Generator as G
    N1 = a1[2]; N2 = a2[2]
    diff = N2-N1
    if diff/2 != diff*0.5: return ['N1-N2 must be even.']
    if diff < 0: ap = a2; a2 = a1; a1 = ap; diff = -diff; N2 = N1
    Np = (diff+2)/2

    b1 = T.subzone(a2, (1,1,1), (Np,1,1))
    b2 = T.subzone(a2, (N2-Np+1,1,1), (N2,1,1))
    b3 = T.subzone(a2, (Np,1,1), (N2-Np+1,1,1))
    #C.convertArrays2File([a1,b1,b2,b3], 'borders.plt')

    m1 = G.TFI([a1,b1,b3,b2])
    return [m1]

#==============================================================================
def step1():
    if CTK.t == []: return

    # Recupere le profil
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    if len(nzs) > 2:
        CTK.TXT.insert('START', 'Input profile must be one curve or two curves (blunt profiles).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    if len(nzs) == 2: culot = 1
    else: culot = 0

    zones = []; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Unstructured':
            try: z = C.convertBAR2Struct(z)
            except Exception as e:
                #print('Error: blader: %s'%str(e))
                errors += [0,str(e)]
                CTK.TXT.insert('START', 'Input profile must be structured.\n')
                CTK.TXT.insert('START', 'Error: ', 'Error'); return
        zones.append(z)
    if len(errors)>0: Panels.displayErrors(errors, header='Error: blader')
    CTK.saveTree()

    # -- Go to array world!
    import Geom as D
    import Generator as G
    import Transform as T
    import Converter

    # repere le culot si 2 courbes sont fournies
    a = C.getAllFields(zones[0], 'nodes')[0]
    if culot == 1:
        ac = C.getAllFields(zones[1], 'nodes')[0]
        bb1 = G.bbox(a); bb2 = G.bbox(ac)
        if bb1[0] > bb2[0]: temp = a; a = ac; ac = temp

    # taille de maille trailing edge et culot
    h = float(VARS[1].get())
    # deraffinement par rapport a la distribution initiale
    h2 = float(VARS[2].get())
    # Point de split (pour la ligne medianne), en % de la longeur du profil
    Nsplit = float(VARS[0].get())
    # Creation delta (Ndelta est en nbre de pts sur le profil remaille)
    Ndelta = int(VARS[3].get())
    # Hauteur du maillage
    Dfar = float(VARS[4].get())

    #==========================================================================
    # Remaille uniforme du profil avec h2
    l = D.getLength(a)
    npts = int(l / h2)+1
    if npts/2 == npts*0.5: npts += 1
    distrib = G.cart( (0,0,0), (1./(npts-1.),1,1), (npts,1,1) )
    a = G.map(a, distrib)

    #===========================================================================
    # Split du profil en intrados/extrados
    N = a[2]; Ns = int(N*Nsplit)
    a1 = T.subzone( a, (1,1,1), (Ns+1,1,1) )
    a2 = T.subzone( a, (Ns+1,1,1), (N,1,1) )
    a2 = T.reorder(a2, (-1,2,3))

    #===========================================================================
    # Resserement bord de fuite et bord d'attaque
    l = D.getLength(a1)
    s = D.getCurvilinearAbscissa(a1); s[0] = 'x'
    s = Converter.initVars(s, 'y', 0.)
    s = Converter.initVars(s, 'z', 0.)
    N = a1[2]
    distrib = G.enforcePlusX(s, h/l, N/10, 1)
    distrib = G.enforceMoinsX(distrib, h/l, N/10, 1)
    a1 = G.map(a1, distrib)

    l = D.getLength(a2)
    s = D.getCurvilinearAbscissa(a2); s[0] = 'x'
    s = Converter.initVars(s, 'y', 0.)
    s = Converter.initVars(s, 'z', 0.)
    N = a2[2]
    distrib = G.enforcePlusX(s, h/l, N/10, 1)
    distrib = G.enforceMoinsX(distrib, h/l, N/10, 1)
    a2 = G.map(a2, distrib)

    #==========================================================================
    # ligne delta
    #===========================================================================
    ni = a1[2]
    b1 = T.subzone(a1, (1, 1,1), (Ndelta,1,1))
    c1 = T.subzone(a1, (Ndelta,1,1), (ni,1,1))
    b2 = T.subzone(a2, (1,1,1), (Ndelta,1,1))
    c2 = T.subzone(a2, (Ndelta,1,1), (ni,1,1))
    P1 = (c1[1][0,0], c1[1][1,0], c1[1][2,0])
    P2 = (c2[1][0,0], c2[1][1,0], c2[1][2,0])
    NdeltaDelta = Ndelta-4
    delta = D.line(P1, P2, 2*NdeltaDelta-1)

    #===========================================================================
    # ligne medianne
    #===========================================================================
    # ligne medianne droite
    P1 = Converter.getValue(c1,c1[2]-1)
    P2 = Converter.getValue(delta,delta[2]/2)
    median = D.line(P1, P2, N=10)

    # ligne medianne moyenne
    N1 = c1[2]
    # cree une ligne k2 remaillant c2 avec N1
    s = D.getCurvilinearAbscissa(c1); s[0] = 'x'
    s = Converter.initVars(s, 'y', 0.)
    s = Converter.initVars(s, 'z', 0.)
    k2 = G.map(c2, s)
    # cree la ligne moyenne
    median = Converter.copy(c1)
    median[1][0,:] = 0.5*(c1[1][0,:]+k2[1][0,:])
    median[1][1,:] = 0.5*(c1[1][1,:]+k2[1][1,:])
    median[1][2,:] = 0.5*(c1[1][2,:]+k2[1][2,:])

    # Remaillage ligne medianne
    s = D.getCurvilinearAbscissa(c1); s[0] = 'x'
    s = Converter.initVars(s, 'y', 0.)
    s = Converter.initVars(s, 'z', 0.)

    median = G.map(median, s)
    median = G.refine(median, 0.9, 1)
    N1 = c1[2]; N2 = median[2]; d = N1-N2
    if d/2 != d*0.5:
        factor = (N2+2.)/N2
        median = G.refine(median, factor, 1)

    #===========================================================================
    # Maillage TRI au bout
    if culot == 0:
        #Converter.convertArrays2File([b1,b2,delta], 'bout1.plt')
        m3 = trimesh(b1, b2, delta)
        if m3[0] == 0: raise ValueError(m3[1])
        #Converter.convertArrays2File([b1,b2,delta]+m3, 'bout1.plt')
    else:
        # Dans le cas avec culot, on remaille le culot comme delta
        l = D.getLength(ac)
        npts = delta[2]
        distrib = G.cart( (0,0,0), (1./(npts-1.),1,1), (npts,1,1) )
        ac = G.map(ac, distrib)
        m3 = [G.TFI([b1,delta,b2,ac])]

    #===========================================================================
    # Maillage du reste de l'interieur du profil
    ni = c1[2]
    d1 = T.subzone(c1, (1,1,1), (ni-NdeltaDelta+1,1,1))
    e1 = T.subzone(c1, (ni-NdeltaDelta+1,1,1), (ni,1,1))
    ni = c2[2]
    d2 = T.subzone(c2, (1,1,1), (ni-NdeltaDelta+1,1,1))
    e2 = T.subzone(c2, (ni-NdeltaDelta+1,1,1), (ni,1,1))

    # remaillage de l'axe median
    s = D.getCurvilinearAbscissa(d1); s[0] = 'x'
    s = Converter.initVars(s, 'y', 0.)
    s = Converter.initVars(s, 'z', 0.)
    median = G.map(median, s)

    npt = delta[2]
    delta1 = T.subzone(delta, (1,1,1), (npt/2+1,1,1))
    delta2 = T.subzone(delta, (npt/2+1,1,1), (npt,1,1))
    #Converter.convertArrays2File([d1,e1,delta1,delta2,median,d2,e2], 'curv1.plt')

    m1 = G.TFI([d1,e1,delta1,median])
    m2 = G.TFI([d2,e2,delta2,median])

    #===========================================================================
    # Ajout de la ligne arriere

    P1 = Converter.getValue(a1,0)
    P2 = Converter.getValue(a1,1)
    P3 = (P1[0]+Dfar,P1[1],P1[2])
    line2 = D.line(P1, P3, N=50)

    N = 50
    s = G.cart( (0,0,0), (1./(N-1.),1,1), (N,1,1))
    s = G.enforcePlusX(s, abs(P2[0]-P1[0])/Dfar, N, 1)
    line2 = G.map(line2, s)
    #Converter.convertArrays2File([a1,a2,line2], 'out.plt')

    if culot == 0:
        line2p = Converter.copy(line2)
    else:
        P1 = Converter.getValue(a2,0)
        P2 = Converter.getValue(a2,1)
        P3 = (P1[0]+Dfar,P1[1],P1[2])
        line2p = D.line(P1, P3, N=50)
        s = D.getCurvilinearAbscissa(line2); s[0] = 'x'
        s = Converter.initVars(s, 'y', 0.)
        s = Converter.initVars(s, 'z', 0.)
        line2p = G.map(line2p, s)

    #===========================================================================
    # Add normal layers
    N = 50; h = abs(P2[0]-P1[0])
    N = int(Dfar/(3*h)/2)+1
    d = G.cart( (0,0,0), (3*h,1,1), (N+1,1,1))

    # shift
    line2p[1][1,:] -= 1.e-9
    a2[1][1,0] -= 1.e-9
    #Converter.convertArrays2File([a2,line2,a1,line2p], 'curve0.plt')

    curve0 = T.join([line2, a1, a2, line2p])
    curve0 = T.reorder(curve0, (-1,2,3))
    #Converter.convertArrays2File([curve0], 'curve.plt')

    m4 = G.addNormalLayers([curve0], d, niter=1000, check=0)
    m4 = G.close(m4, 1.e-8)
    #Converter.convertArrays2File([m1,m2]+m3+m4, 'all.plt')

    # check volume + subzone
    vol = G.getVolumeMap(m4[0])
    nk = vol[4]; ni = vol[2]
    for k in range(nk-1):
        sub = T.subzone(vol, (1,1,k+1), (ni,1,k+2))
        volmin = Converter.getMinValue(sub, 'vol')
        if volmin < 0.:
            if k > 10: kc = k-4
            else: kc = k-1
            m4[0] = T.subzone(m4[0], (1,1,1), (m4[0][2],1,kc))
            break

    if culot == 1:
        r1p = T.translate(ac, (Dfar,0,0))
        Converter.convertArrays2File([ac,line2,r1p,line2p], 'sub.plt')
        bsup = G.TFI([ac,line2,r1p,line2p])

    # interieur
    M1 = [m1,m2]+m3
    # exterieur
    M2 = m4
    M2 = T.reorder(M2, (1,3,2))

    # split du bloc M2 comme les blocs interieurs (pour connectMatch)
    a = M2[0]
    mp = line2[2]+a1[2]-1
    i1 = T.subzone(a, (1,1,1), (mp,a[3],a[4]))
    i2 = T.subzone(a, (mp,1,1), (a[2],a[3],a[4]))
    M2 = [i1, i2]

    if culot == 1: M2 += [bsup]

    #==========================================================================
    # Zones resultantes dans l'arbre
    meshes = M1+M2
    zones = []
    for m in meshes:
        z = C.convertArrays2ZoneNode('zone', [m])
        zones.append(z)

    base = Internal.getNodesFromName1(CTK.t, 'STEP1')
    if base != []:
        (p, c) = Internal.getParentOfNode(CTK.t, base[0])
        del p[2][c]

    CTK.t = C.addBase2PyTree(CTK.t, 'STEP1', 3)
    base = Internal.getNodesFromName1(CTK.t, 'STEP1')[0]
    base[2] += zones

    CTK.TXT.insert('START', 'Step1 performed.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Deuxieme etape: passage au maillage volumique
#==============================================================================
def step2():

    import Transform as T
    import Generator as G

    CTK.saveTree()

    # taille de maille en envergure
    hp = float(VARS[5].get())
    # Hauteur du maillage
    Dfar = float(VARS[4].get())
    # Envergure
    span = float(VARS[6].get())

    # Recupere la base STEP1
    base = Internal.getNodesFromName1(CTK.t, 'STEP1')
    if base == []: return
    zones = base[0][2]
    l = len(zones)
    if l == 7: culot = 0
    else: culot = 1
    if culot == 0:
        # 2 zones exterieures, le reste interieur
        M2 = [C.getAllFields(zones[l-2], 'nodes')[0],
              C.getAllFields(zones[l-1], 'nodes')[0]]
        M1 = []
        for z in zones[0:l-2]:
            M1.append(C.getAllFields(z, 'nodes')[0])
    else:
        # 3 zones exterieures, le reste interieur
        M2 = [C.getAllFields(zones[l-2], 'nodes')[0],
              C.getAllFields(zones[l-1], 'nodes')[0],
              C.getAllFields(zones[l-3], 'nodes')[0]]
        M1 = []
        for z in zones[0:l-3]:
            M1.append(C.getAllFields(z, 'nodes')[0])

    #==========================================================================
    # stack + resserement vers les extremites
    #==========================================================================
    M1b = T.translate(M1, (0,0,Dfar))
    B1 = []
    for i in range(len(M1)):
        B1.append(G.stack(M1[i], M1b[i]))

    M1c = T.translate(M1, (0,0,-span))
    M1d = T.translate(M1, (0,0,-span-Dfar))
    B2 = []
    for i in range(len(M1c)):
        B2.append(G.stack(M1c[i], M1d[i]))

    #C.convertArrays2File(B1+B2, 'bouchon.plt')

    M2b = T.translate(M2, (0,0,Dfar))
    M2c = T.translate(M2, (0,0,-span-Dfar))
    I = []
    for i in range(len(M2b)):
        I.append(G.stack(M2c[i], M2b[i]))

    # B1, B2: les bouchons; I le reste
    #C.convertArrays2File(B1+B2+I, 'all.plt')

    #==========================================================================
    # Remaille la surface

    N = int(Dfar/hp)+1
    distrib = G.cart( (0,0,0), (1./(N-1),1,1), (N,1,1) )
    for i in range(len(B1)):
        B1[i] = G.map(B1[i], distrib, 3)
    for i in range(len(B2)):
        B2[i] = G.map(B2[i], distrib, 3)
    N = int((2*Dfar+span)/hp)+1
    distrib = G.cart( (0,0,0), (1./(N-1),1,1), (N,1,1) )
    for i in range(len(I)):
        I[i] = G.map(I[i], distrib, 3)

    # Back to zones
    zones = []
    for b in B1+B2+I:
        zones.append(C.convertArrays2ZoneNode('zone', [b]))

    base = Internal.getNodesFromName1(CTK.t, 'STEP2')
    if base != []:
        (p, c) = Internal.getParentOfNode(CTK.t, base[0])
        del p[2][c]

    CTK.t = C.addBase2PyTree(CTK.t, 'STEP2', 3)
    base = Internal.getNodesFromName1(CTK.t, 'STEP2')[0]
    (p, c) = Internal.getParentOfNode(CTK.t, base)
    base[2] += zones

    # Add BCs
    base = X.connectMatch(base, tol=1.e-6)

    # # Blocs exterieurs
    if culot == 0:
        z = base[2][10]
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'jmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmin')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'imin')
        base[2][10] = z

        z = base[2][11]
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'jmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmin')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'imax')
        base[2][11] = z

        for i in range(5):
            z = base[2][i]
            z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
            base[2][i] = z

        for i in range(5):
            z = base[2][5+i]
            z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
            base[2][5+i] = z
    else:
        z = base[2][6]
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'jmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmin')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'imax')
        base[2][6] = z

        z = base[2][8]
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'jmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmin')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'imin')
        base[2][8] = z

        z = base[2][7]
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'imax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmin')
        base[2][7] = z

        for i in range(3):
            z = base[2][i]
            z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
            base[2][i] = z

        for i in range(3):
            z = base[2][3+i]
            z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
            base[2][3+i] = z

    base = C.fillEmptyBCWith(base, 'wall', 'BCWall')
    CTK.t[2][c] = base

    CTK.TXT.insert('START', 'Step2 performed.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkBlader  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Automatic mesher for\nblades.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkBlader')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- front split % -
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    if 'tkBladerFrontSplit' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerFrontSplit'])
    # -1- front step -
    V = TK.StringVar(win); V.set('0.001'); VARS.append(V)
    if 'tkBladerFrontStep' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerFrontStep'])
    # -2- other step -
    V = TK.StringVar(win); V.set('0.01'); VARS.append(V)
    if 'tkBladerStep' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerStep'])
    # -3- delta line index -
    V = TK.StringVar(win); V.set('15'); VARS.append(V)
    # -4- Dfar. Mesh height -
    V = TK.StringVar(win); V.set('0.3'); VARS.append(V)
    if 'tkBladerHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerHeight'])
    # -5- hp: step en envergure -
    V = TK.StringVar(win); V.set('0.02'); VARS.append(V)
    if 'tkBladerSpanStep' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerSpanStep'])
    # -6- span: longeur de l'envergure -
    V = TK.StringVar(win); V.set('5.'); VARS.append(V)
    if 'tkBladerSpan' in CTK.PREFS:
        V.set(CTK.PREFS['tkBladerSpan'])

    # - Step1 -
    B = TTK.Label(Frame, text='Front split %:')
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split profile at %.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Index for front split.')

    B = TTK.Label(Frame, text='Front step:')
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Front step.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Step for front and back of profile.')

    B = TTK.Label(Frame, text='Profile step:')
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Profile step.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Step on profile.')

    B = TTK.Label(Frame, text='Delta line index:')
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Index for delta line.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Index for delta line.')

    B = TTK.Label(Frame, text='Mesh height:')
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Mesh height.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Mesh height.')

    B = TTK.Button(Frame, text="Step 1", command=step1)
    B.grid(row=5, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Run step 1. Plane mesh.')

    B = TTK.Label(Frame, text='Span length:')
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Span length.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=5)
    B.grid(row=6, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Span length.')

    B = TTK.Label(Frame, text='Span step:')
    B.grid(row=7, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Step in span direction.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=5)
    B.grid(row=7, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Step in span direction.')

    B = TTK.Button(Frame, text="Step 2", command=step2)
    B.grid(row=8, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Run step 2. Volume mesh.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkBlader')
    except: pass
    CTK.WIDGETS['MeshNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['MeshNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkBladerFrontSplit'] = VARS[0].get()
    CTK.PREFS['tkBladerFrontStep'] = VARS[1].get()
    CTK.PREFS['tkBladerStep'] = VARS[2].get()
    CTK.PREFS['tkBladerHeight'] = VARS[4].get()
    CTK.PREFS['tkBladerSpanStep'] = VARS[5].get()
    CTK.PREFS['tkBladerSpan'] = VARS[6].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('0.5')
    VARS[1].set('0.001')
    VARS[2].set('0.01')
    VARS[4].set('0.3')
    VARS[5].set('0.02')
    VARS[6].set('5.')
    CTK.PREFS['tkBladerFrontSplit'] = VARS[0].get()
    CTK.PREFS['tkBladerFrontStep'] = VARS[1].get()
    CTK.PREFS['tkBladerStep'] = VARS[2].get()
    CTK.PREFS['tkBladerHeight'] = VARS[4].get()
    CTK.PREFS['tkBladerSpanStep'] = VARS[5].get()
    CTK.PREFS['tkBladerSpan'] = VARS[6].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkBlader '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
