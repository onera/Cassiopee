# Cartesian generators
from . import PyTree as G
try:
    import Converter.Mpi as Cmpi
    import Converter.Internal as Internal
    import Converter.PyTree as C
    import Distributor2.PyTree as D2
    import Transform.PyTree as T
except ImportError:
    def cartRx(X0, H, N, Nb, depth=0, addCellN=False, addBCMatch=False, rank=None, size=None): return None
    def _cartRxRefit(a): return None
    def cartRx2(XC0, XC1, HC, XF0, XF1, R, rank=None, size=None): return None
    def cartRx3(XC0, XC1, HC, XF0, XF1, R, rank=None, size=None): return None
    def cartRxHollow(XC0, XC1, HC, XH0, XH1, XF0, XF1, R, rank=None, size=None): return None


def cartRx(X0, H, N, Nb, depth=0, addCellN=False, addBCMatch=False,
           rank=None, size=None):
    """Create a set of regular cartesian grids."""
    out = []
    for k in range(Nb[2]):
        for j in range(Nb[1]):
            for i in range(Nb[0]):
                if rank is None or size is None or rank == (i+j*Nb[0]+k*Nb[0]*Nb[1])%size:
                    Xp = [X0[0]+H[0]*(N[0]-1)*i,X0[1]+H[1]*(N[1]-1)*j,X0[2]+H[2]*(N[2]-1)*k]
                    Np = [N[0],N[1],N[2]]
                    if i > 0: Xp[0] -= depth*H[0]; Np[0] += depth
                    if i < Nb[0]-1: Xp[0] += depth*H[0]; Np[0] += depth
                    if j > 0: Xp[1] -= depth*H[1]; Np[1] += depth
                    if j < Nb[1]-1: Xp[1] += depth*H[1]; Np[1] += depth
                    if k > 0: Xp[2] -= depth*H[2]; Np[2] += depth
                    if k < Nb[2]-1: Xp[2] += depth*H[2]; Np[2] += depth
                    z = G.cart(Xp, H, Np); z[0] = 'cart%d-%d-%d'%(i,j,k)
                    if rank is not None: Cmpi._setProc(z, rank)
                    if addCellN:
                        C._initVars(z, 'centers:cellN', 1)
                        cellN = Internal.getNodeFromName2(z, 'cellN')[1]
                        if i > 0: cellN[0:depth,:,:] = 2
                        if i < Nb[0]-1: cellN[Np[0]-depth-1:Np[0]-1,:,:] = 2
                        if j > 0: cellN[:,0:depth,:] = 2
                        if j < Nb[1]-1: cellN[:,Np[1]-depth-1:Np[1]-1,:] = 2
                        if k > 0: cellN[:,:,0:depth] = 2
                        if k < Nb[2]-1: cellN[:,:,Np[2]-depth-1:Np[2]-1] = 2
                    if addBCMatch and depth == 0:
                        if i > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i-1,j,k))
                        if i < Nb[0]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i+1,j,k))
                        if j > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z, 'jmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j-1,k))
                        if j < Nb[1]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j+1,k))
                        if k > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k-1))
                        if k < Nb[2]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k+1))
                    out.append(z)
    return out

def _cartRxRefit(a):
    """Refit a cartRx to cartesian after addGhostCells."""
    zones = Internal.getZones(a)
    for z in zones:
        dim = Internal.getZoneDim(z)
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        xp = Internal.getNodeFromName2(z, 'CoordinateX')[1]
        if ni > 3:
            dx1 = xp[1,0,0]-xp[0,0,0]
            dx2 = xp[2,0,0]-xp[1,0,0]
            dx3 = xp[3,0,0]-xp[2,0,0]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[0,:,:] = xp[2,:,:]-2*dx3
                xp[1,:,:] = xp[2,:,:]-dx3
            dx1 = xp[-1,0,0]-xp[-2,0,0]
            dx2 = xp[-2,0,0]-xp[-3,0,0]
            dx3 = xp[-3,0,0]-xp[-4,0,0]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[-1,:,:] = xp[-3,:,:]+2*dx3
                xp[-2,:,:] = xp[-3,:,:]+dx3
        xp = Internal.getNodeFromName2(z, 'CoordinateY')[1]
        if nj > 3:
            dx1 = xp[0,1,0]-xp[0,0,0]
            dx2 = xp[0,2,0]-xp[0,1,0]
            dx3 = xp[0,3,0]-xp[0,2,0]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[:,0,:] = xp[:,2,:]-2*dx3
                xp[:,1,:] = xp[:,2,:]-dx3
            dx1 = xp[0,-1,0]-xp[0,-2,0]
            dx2 = xp[0,-2,0]-xp[0,-3,0]
            dx3 = xp[0,-3,0]-xp[0,-4,0]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[:,-1,:] = xp[:,-3,:]+2*dx3
                xp[:,-2,:] = xp[:,-3,:]+dx3
        xp = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
        if nk > 3:
            dx1 = xp[0,0,1]-xp[0,0,0]
            dx2 = xp[0,0,2]-xp[0,0,1]
            dx3 = xp[0,0,3]-xp[0,0,2]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[:,:,0] = xp[:,:,2]-2*dx3
                xp[:,:,1] = xp[:,:,2]-dx3
            dx1 = xp[0,0,-1]-xp[0,0,-2]
            dx2 = xp[0,0,-2]-xp[0,0,-3]
            dx3 = xp[0,0,-3]-xp[0,0,-4]
            if abs(dx1) < 1.e-12 and abs(dx2) < 1.e-12:
                xp[:,:,-1] = xp[:,:,-3]+2*dx3
                xp[:,:,-2] = xp[:,:,-3]+dx3
    return None

def cartRx2(XC0, XC1, HC, XF0, XF1, R, dim=3, rank=None, size=None):
    """Create a set of regular and geometric cartesian grids."""

    L0x = XC0[0]-XF0[0]
    L1x = XC1[0]-XC0[0]
    L2x = XF1[0]-XC1[0]
    L0y = XC0[1]-XF0[1]
    L1y = XC1[1]-XC0[1]
    L2y = XF1[1]-XC1[1]
    L0z = XC0[2]-XF0[2]
    L1z = XC1[2]-XC0[2]
    L2z = XF1[2]-XC1[2]

    X0x = [XC0[0], XC0[0], XC0[0]+L1x]
    X0y = [XC0[1], XC0[1], XC0[1]+L1y]
    X0z = [XC0[2], XC0[2], XC0[2]+L1z]
    X1x = [XC0[0]-L0x, XC0[0]+L1x, XC0[0]+L1x+L2x]
    X1y = [XC0[1]-L0y, XC0[1]+L1y, XC0[1]+L1y+L2y]
    X1z = [XC0[2]-L0z, XC0[2]+L1z, XC0[2]+L1z+L2z]
    Rx = [R[0],1.,R[0]]
    Ry = [R[1],1.,R[1]]
    Rz = [R[2],1.,R[2]]

    a = [0]* (3*3*3)
    dimj = 3; dimk1 = 0; dimk2 = 3
    if dim == 2: dimk1 = 1; dimk2 = 2

    # squelette
    data = {}
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(dimk1, dimk2):
                Px = X0x[i]; Py = X0y[j]; Pz = X0z[k]
                Qx = X1x[i]; Qy = X1y[j]; Qz = X1z[k]
                (ni,nj,nk,rio,rjo,rko,hio,hjo,hko) = G.cartr2((Px,Py,Pz), HC, (Rx[i],Ry[j],Rz[k]), (Qx,Qy,Qz), skeleton=True)
                z = Internal.newZone('Zone', zsize=[[ni,ni-1,0], [nj,nj-1,0], [nk,nk-1,0]], ztype='Structured')
                n = Internal.newGridCoordinates(parent=z)
                Internal.newDataArray('CoordinateX', value=None, parent=n)
                z[0] = 'cart%d-%d-%d'%(i,j,k)

                if i > 0:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i-1,j,k))
                if i < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i+1,j,k))
                if j > 0:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z,'jmax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j-1,k))
                if j < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j+1,k))
                if k > 0 and dim == 3:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k-1))
                if k < 2 and dim == 3:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k+1))

                if i == 0: hio=hio*rio**(ni-2); rio=1./rio; Px=Qx
                if j == 0: hjo=hjo*rjo**(nj-2); rjo=1./rjo; Py=Qy
                if k == 0: hko=hko*rko**(nk-2); rko=1./rko; Pz=Qz

                data[z[0]] = [(Px,Py,Pz), (hio,hjo,hko), (rio,rjo,rko), (ni,nj,nk)]

                a[i+3*j+9*k] = z

    if dim == 2: # clean list for 2D case
        out = []
        for i in a:
            if i != 0: out.append(i)
        a = out

    t = C.newPyTree(['CARTESIAN','FLEX'])
    for z in a:
        if z[0] == 'cart1-1-1': core = z; break
    t[2][1][2].append(core)
    a.remove(core)
    t[2][2][2] += a

    # correction des fenetres max des BCs
    for z in Internal.getZones(t):
        bcs = Internal.getNodesByType(z, 'GridConnectivity1to1_t')
        for BC in bcs:
            PtRangeDonor = Internal.getNodeFromName1(BC, 'PointRangeDonor')[1]
            donorName = Internal.getValue(BC)
            zd = Internal.getNodeFromName2(t, donorName)
            dimz = Internal.getZoneDim(zd)
            imaxDonor = dimz[1]; jmaxDonor = dimz[2]; kmaxDonor = dimz[3]
            if PtRangeDonor[0,0] > 1 and PtRangeDonor[0,0] == PtRangeDonor[0,1]:
                PtRangeDonor[0,0] = imaxDonor
                PtRangeDonor[0,1] = imaxDonor
            if PtRangeDonor[1,0] > 1 and PtRangeDonor[1,0] == PtRangeDonor[1,1]:
                PtRangeDonor[1,0] = jmaxDonor
                PtRangeDonor[1,1] = jmaxDonor
            if PtRangeDonor[2,0] > 1 and PtRangeDonor[2,0] == PtRangeDonor[2,1]:
                PtRangeDonor[2,0] = kmaxDonor
                PtRangeDonor[2,1] = kmaxDonor

    # SplitNParts on core
    if size is None: size = 1
    if rank is None: rank = 0    
    b = Internal.getNodeFromName(t, 'CARTESIAN')
    T._splitNParts(b, N=size, topTree=t)
    D2._distribute(b, NProc=size, algorithm='fast')

    # SplitSize + ressource : distribue en meme temps
    b = Internal.getNodeFromName(t, 'FLEX')
    T._splitSize(b, R=size, topTree=t)
    #D2._distribute(t2, NProc=size, algorithm='fast') # deja fait par splitSize
    D2.printStats(b)

    # Generation reelle
    bases = Internal.getBases(t)
    for b in bases:
        for c in range(len(b[2])):
            z = b[2][c]

            if z[3] == 'Zone_t' and Cmpi.getProc(z) == rank:
                if z[0] in data: # bloc non splitte
                    #print(z[0],'bloc non splite', flush=True)
                    d = data[z[0]]
                    zn = G.cartr1(d[0], d[1], d[2], d[3])
                else:
                    #print(z[0],'bloc splitte', flush=True)
                    source, dest = Internal.getLoc2Glob(z)
                    d = data[source]
                    #print('source', source, flush=True)
                    #print('dest', dest, flush=True)
                    P = d[0]; H = d[1]; R = d[2] ; N = d[3]
                    i1 = dest[0]-1; j1 = dest[2]-1; k1 = dest[4]-1
                    i2 = dest[1]-1; j2 = dest[3]-1; k2 = dest[5]-1

                    if R[0] == 1.: 
                        ratiox = i1; Hx = H[0]
                    else:
                        ratiox = (R[0]**i1-1.)/(R[0]-1.); Hx = H[0]*R[0]**i1
                    if R[1] == 1.: 
                        ratioy = j1; Hy = H[1]
                    else:
                        ratioy = (R[1]**j1-1.)/(R[1]-1.); Hy = H[1]*R[1]**j1
                    if R[2] == 1.: 
                        ratioz = k1; Hz = H[2]
                    else:
                        ratioz = (R[2]**k1-1.)/(R[2]-1.); Hz = H[2]*R[2]**k1
                    Px = P[0] + ratiox*H[0]
                    Py = P[1] + ratioy*H[1]
                    Pz = P[2] + ratioz*H[2]
                    Rx = R[0]; Ry = R[1]; Rz = R[2]
                    N = (i2-i1+1,j2-j1+1,k2-k1+1)
                    zn = G.cartr1((Px,Py,Pz), (Hx,Hy,Hz), (Rx,Ry,Rz), N)

                zn[0] = z[0]
                D2._addProcNode(zn, rank)
                n = Internal.getNodesFromName(z, 'ZoneBC')
                zn[2] += n
                n = Internal.getNodesFromName(z, 'ZoneGridConnectivity')
                zn[2] += n
                b[2][c] = zn

    Cmpi._convert2PartialTree(t)
    return t

def cartRx3(XC0, XC1, HC, XF0, XF1, R, dim=3, rank=None, size=None):
    """Create a set of regular and geometric cartesian grids with double steps."""

    L0x = XC0[0]-XF0[0]
    L1x = XC1[0]-XC0[0]
    L2x = XF1[0]-XC1[0]
    L0y = XC0[1]-XF0[1]
    L1y = XC1[1]-XC0[1]
    L2y = XF1[1]-XC1[1]
    L0z = XC0[2]-XF0[2]
    L1z = XC1[2]-XC0[2]
    L2z = XF1[2]-XC1[2]

    X0x = [XC0[0], XC0[0], XC0[0]+L1x]
    X0y = [XC0[1], XC0[1], XC0[1]+L1y]
    X0z = [XC0[2], XC0[2], XC0[2]+L1z]
    X1x = [XC0[0]-L0x, XC0[0]+L1x, XC0[0]+L1x+L2x]
    X1y = [XC0[1]-L0y, XC0[1]+L1y, XC0[1]+L1y+L2y]
    X1z = [XC0[2]-L0z, XC0[2]+L1z, XC0[2]+L1z+L2z]
    Rx = [R[0],1.,R[0]]
    Ry = [R[1],1.,R[1]]
    Rz = [R[2],1.,R[2]]

    a = [0]* (3*3*3)
    dimj = 3; dimk1 = 0; dimk2 = 3
    if dim == 2: dimk1 = 1; dimk2 = 2

    # squelette
    data = {}
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(dimk1, dimk2):
                Px = X0x[i]; Py = X0y[j]; Pz = X0z[k]
                Qx = X1x[i]; Qy = X1y[j]; Qz = X1z[k]
                if i == 0: dli = 1 ; dri = 0
                elif i == 1: dli = 1 ; dri = 1
                elif i == 2: dli = 1 ; dri = 0
                if j == 0: dlj = 1 ; drj = 0
                elif j == 1: dlj = 1 ; drj = 1
                elif j == 2: dlj = 1 ; drj = 0
                if k == 0: dlk = 1 ; drk = 0
                elif k == 1: dlk = 1 ; drk = 1
                elif k == 2: dlk = 1 ; drk = 0
                doubleLeft = [dli,dlj,dlk]
                doubleRight = [dri,drj,drk]
                (ni,nj,nk,rio,rjo,rko,hio,hjo,hko) = G.cartr2((Px,Py,Pz), HC, (Rx[i],Ry[j],Rz[k]), (Qx,Qy,Qz), doubleLeft, doubleRight, skeleton=True)
                z = Internal.newZone('Zone', zsize=[[ni,ni-1,0], [nj,nj-1,0], [nk,nk-1,0]], ztype='Structured')
                n = Internal.newGridCoordinates(parent=z)
                Internal.newDataArray('CoordinateX', value=None, parent=n)
                z[0] = 'cart%d-%d-%d'%(i,j,k)

                if i > 0:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i-1,j,k))
                if i < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i+1,j,k))
                if j > 0:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z,'jmax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j-1,k))
                if j < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j+1,k))
                if k > 0 and dim == 3:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k-1))
                if k < 2 and dim == 3:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k+1))

                if i == 0: hio=hio*rio**(ni-3); rio=1./rio; Px=Qx; doubleLeft[0] = 0; doubleRight[0] = 1
                if j == 0: hjo=hjo*rjo**(nj-3); rjo=1./rjo; Py=Qy; doubleLeft[1] = 0; doubleRight[1] = 1
                if k == 0: hko=hko*rko**(nk-3); rko=1./rko; Pz=Qz; doubleLeft[2] = 0; doubleRight[2] = 1

                data[z[0]] = [(Px,Py,Pz), (hio,hjo,hko), (rio,rjo,rko), (ni,nj,nk), doubleLeft, doubleRight]

                a[i+3*j+9*k] = z

    if dim == 2: # clean list for 2D case
        out = []
        for i in a:
            if i != 0: out.append(i)
        a = out

    t = C.newPyTree(['CARTESIAN','FLEX'])
    for z in a:
        if z[0] == 'cart1-1-1': core = z; break
    t[2][1][2].append(core)
    a.remove(core)
    t[2][2][2] += a

    # correction des fenetres max des BCs
    for z in Internal.getZones(t):
        bcs = Internal.getNodesByType(z, 'GridConnectivity1to1_t')
        for BC in bcs:
            PtRangeDonor = Internal.getNodeFromName1(BC, 'PointRangeDonor')[1]
            donorName = Internal.getValue(BC)
            zd = Internal.getNodeFromName2(t, donorName)
            dimz = Internal.getZoneDim(zd)
            imaxDonor = dimz[1]; jmaxDonor = dimz[2]; kmaxDonor = dimz[3]
            if PtRangeDonor[0,0] > 1 and PtRangeDonor[0,0] == PtRangeDonor[0,1]:
                PtRangeDonor[0,0] = imaxDonor
                PtRangeDonor[0,1] = imaxDonor
            if PtRangeDonor[1,0] > 1 and PtRangeDonor[1,0] == PtRangeDonor[1,1]:
                PtRangeDonor[1,0] = jmaxDonor
                PtRangeDonor[1,1] = jmaxDonor
            if PtRangeDonor[2,0] > 1 and PtRangeDonor[2,0] == PtRangeDonor[2,1]:
                PtRangeDonor[2,0] = kmaxDonor
                PtRangeDonor[2,1] = kmaxDonor

    # SplitNParts on core
    if size is None: size = 1
    if rank is None: rank = 0    
    b = Internal.getNodeFromName(t, 'CARTESIAN')
    T._splitNParts(b, N=size, topTree=t)
    D2._distribute(b, NProc=size, algorithm='fast')

    # SplitSize + ressource : distribue en meme temps
    b = Internal.getNodeFromName(t, 'FLEX')
    T._splitSize(b, R=size, topTree=t)
    #D2._distribute(t2, NProc=size, algorithm='fast') # deja fait par splitSize
    D2.printStats(b)

    # Generation reelle
    bases = Internal.getBases(t)
    for b in bases:
        for c in range(len(b[2])):
            z = b[2][c]

            if z[3] == 'Zone_t' and Cmpi.getProc(z) == rank:
                if z[0] in data: # bloc non splitte
                    #print(z[0], 'bloc non splite', flush=True)
                    d = data[z[0]]
                    zn = G.cartr1(d[0], d[1], d[2], d[3], d[4], d[5])
                else:
                    #print(z[0],'bloc splitte', flush=True)
                    source, dest = Internal.getLoc2Glob(z)
                    d = data[source]
                    #print('source', source, flush=True)
                    # print('dest', dest, flush=True)
                    P = d[0]; H = d[1]; R = d[2] ; N = d[3] ; dL = d[4] ; dR = d[5]
                    i1 = dest[0]-1; j1 = dest[2]-1; k1 = dest[4]-1
                    i2 = dest[1]-1; j2 = dest[3]-1; k2 = dest[5]-1
                    doubleLeft=[0,0,0]; doubleRight=[0,0,0]
                    if dL[0] == 1 and i1 == 0: doubleLeft[0] = 1
                    if dR[0] == 1 and i2 == N[0]-1: doubleRight[0] = 1 
                    if dL[1] == 1 and j1 == 0: doubleLeft[1] = 1
                    if dR[1] == 1 and j2 == N[1]-1: doubleRight[1] = 1 
                    if dL[2] == 1 and k1 == 0: doubleLeft[2] = 1
                    if dR[2] == 1 and k2 == N[2]-1: doubleRight[2] = 1 

                    if doubleLeft[0] == 1 and doubleRight[0] == 0:
                        if R[0] == 1.: ratiox = i1; Hx = H[0]
                        else: ratiox = (R[0]**(i1) - 1. )/(R[0]-1.); Hx= H[0]
                    elif doubleLeft[0] == 0 and doubleRight[0] ==0:
                        if R[0] == 1.: ratiox = i1; Hx = H[0]
                        else:
                            if i1 == 0: ratiox = (R[0]**i1-1.)/(R[0]-1.); Hx = H[0]*R[0]**i1
                            else: ratiox = 1. + (R[0]**(i1-1)-1.)/(R[0]-1.); Hx = H[0]*R[0]**(i1-1)
                    elif doubleLeft[0] == 0 and doubleRight[0] == 1:
                        if R[0] == 1: ratiox=i1; Hx = H[0]
                        else: 
                            if i1 == 0: ratiox = (R[0]**i1 - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**i1
                            else: ratiox = (R[0]**(i1) - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**(i1)                                        
                    elif doubleLeft[0] == 1 and doubleRight[0] == 1:
                        if R[0] == 1: ratiox=i1; Hx = H[0]
                        else: 
                            if i1 == 0: ratiox = (R[0]**i1 - 1.)/(R[0] - 1.); Hx = H[0]
                            else: ratiox =  1.+(R[0]**(i1-1) - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**(i1-1)

                    if doubleLeft[1] == 1 and doubleRight[1] == 0:
                        if R[1] == 1.: ratioy = j1; Hy = H[1]
                        else: ratioy = (R[1]**(j1) - 1.)/(R[1]-1.); Hy= H[1] 
                    elif doubleLeft[1] == 0 and doubleRight[1] == 0:
                        if R[1] == 1.: ratioy = j1; Hy = H[1]
                        else:
                            if j1 == 0: ratioy = (R[1]**j1-1.)/(R[1]-1.); Hy = H[1]*R[1]**j1 
                            else: ratioy = 1. + (R[1]**(j1-1.)-1.)/(R[1]-1.); Hy = H[1]*R[1]**(j1 -1)
                    elif doubleLeft[1] == 0 and doubleRight[1] == 1:
                        if R[1] == 1: ratioy=j1; Hy = H[1]
                        else: 
                            if j1 == 0: ratioy = (R[1]**j1 - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**j1
                            else: ratioy = (R[1]**(j1) - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**(j1) 
                    elif doubleLeft[1] == 1 and doubleRight[1] ==1:
                        if R[1] == 1: ratioy=j1 ; Hy = H[1]
                        else: 
                            if j1 == 0: ratioy = (R[1]**j1 - 1.)/(R[1] - 1.); Hy = H[1]
                            else: ratioy =  1.+(R[1]**(j1-1) - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**(j1-1)

                    if doubleLeft[2] == 1 and doubleRight[2] == 0:
                        if R[2] == 1.: ratioz = k1 ; Hz = H[2]
                        else:
                            if k1 == 0: ratioz = (R[2]**(k1) - 1.)/(R[2]-1.); Hz= H[2]

                    elif doubleLeft[2] == 0 and doubleRight[2] == 0:
                        if R[2] == 1.: ratioz = k1; Hz = H[2]
                        else:
                            if k1 == 0: ratioz = (R[2]**k1-1.)/(R[2]-1.); Hz = H[2]*R[2]**k1                    
                            else: ratioz = 1. + (R[2]**(k1-1) -1.)/(R[2]-1.); Hz = H[2]*R[2]**(k1 - 1) 
                    elif doubleLeft[2] == 0 and doubleRight[2] ==1:
                        if R[2] == 1: ratioz=k1; Hz = H[2]
                        else: 
                            if k1 == 0: ratioz = (R[2]**k1 - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**k1
                            else: ratioz =  (R[2]**(k1) - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**(k1)     

                    elif doubleLeft[2] == 1 and doubleRight[2] == 1:
                        if R[2] == 1: ratioz=k1; Hz = H[2]
                        else: 
                            if k1 == 0: ratioz = (R[2]**k1 - 1.)/(R[2] - 1.); Hz = H[2]
                            else: ratioz =  1.+(R[2]**(k1-1) - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**(k1-1)

                    Px = P[0] + ratiox*H[0]
                    Py = P[1] + ratioy*H[1]
                    Pz = P[2] + ratioz*H[2]
                    Rx = R[0]; Ry = R[1]; Rz = R[2]
                    N = (i2-i1+1,j2-j1+1,k2-k1+1)
                    zn = G.cartr1((Px,Py,Pz), (Hx,Hy,Hz), (Rx,Ry,Rz), N, doubleLeft, doubleRight)

                zn[0] = z[0]
                D2._addProcNode(zn, rank)
                n = Internal.getNodesFromName(z, 'ZoneBC')
                zn[2] += n
                n = Internal.getNodesFromName(z, 'ZoneGridConnectivity')
                zn[2] += n
                b[2][c] = zn

    Cmpi._convert2PartialTree(t)
    return t

def cartRxHollow(XC0, XC1, HC, XH0, XH1, XF0, XF1, R, dim=3, rank=None, size=None):
    """Create a set of regular and geometric cartesian grids with double steps."""

    if abs(XH0[0]) >= abs(XC0[0]) or abs(XH1[0]) >= abs(XC1[0]): raise ValueError ("cartRxHollow: Hollow is bigger than cart on x-axis")
    if abs(XH0[1]) >= abs(XC0[1]) or abs(XH1[1]) >= abs(XC1[1]): raise ValueError ("cartRxHollow: Hollow is bigger than cart on y-axis")
    if abs(XH0[2]) >= abs(XC0[2]) or abs(XH1[2]) >= abs(XC1[2]): raise ValueError ("cartRxHollow: Hollow is bigger than cart on z-axis")

    L0x = XC0[0]-XF0[0]
    L1x = XC1[0]-XC0[0]
    L2x = XF1[0]-XC1[0]
    L0y = XC0[1]-XF0[1]
    L1y = XC1[1]-XC0[1]
    L2y = XF1[1]-XC1[1]
    L0z = XC0[2]-XF0[2]
    L1z = XC1[2]-XC0[2]
    L2z = XF1[2]-XC1[2]

    H0x = XH0[0] - XC0[0]
    H1x = XH1[0] - XH0[0]
    H2x = XC1[0] - XH1[0]
    H0y = XH0[1] - XC0[1]
    H1y = XH1[1] - XH0[1]
    H2y = XC1[1] - XH1[1]
    H0z = XH0[2] - XC0[2]
    H1z = XH1[2] - XH0[2]
    H2z = XC1[2] - XH1[2]

    X0x = [XC0[0], XC0[0], XC0[0]+H0x, XH0[0]+H1x, XC0[0]+L1x]
    X0y = [XC0[1], XC0[1], XC0[1]+H0y, XH0[1]+H1y, XC0[1]+L1y]
    X0z = [XC0[2], XC0[2], XC0[2]+H0z, XH0[2]+H1z, XC0[2]+L1z]
    X1x = [XC0[0]-L0x, XC0[0]+H0x, XH0[0]+H1x, XC0[0]+L1x, XC0[0]+L1x+L2x]
    X1y = [XC0[1]-L0y, XC0[1]+H0y, XH0[1]+H1y, XC0[1]+L1y, XC0[1]+L1y+L2y]
    X1z = [XC0[2]-L0z, XC0[2]+H0z, XH0[2]+H1z, XC0[2]+L1z, XC0[2]+L1z+L2z]

    Rx = [R[0],1.,1.,1.,R[0]]
    Ry = [R[1],1.,1.,1.,R[1]]
    Rz = [R[2],1.,1.,1.,R[2]]

    a = [0]* (5*5*5)
    dimj = 5; dimk1 = 0; dimk2 = 5
    if dim == 2: dimk1 = 2; dimk2 = 3
    coreL = []
    # squelette
    data = {}
    for i in range(0, 5):
        for j in range(0, 5):
            for k in range(dimk1, dimk2):
                Px = X0x[i]; Py = X0y[j]; Pz = X0z[k]
                Qx = X1x[i]; Qy = X1y[j]; Qz = X1z[k]
                if i == 0: dli = 1 ; dri = 0
                elif 1 <= i <= 3: dli = 1 ; dri = 1
                elif i == 4: dli = 1 ; dri = 0
                if j == 0: dlj = 1 ; drj = 0
                elif 1 <= j <= 3: dlj = 1 ; drj = 1
                elif j == 4: dlj = 1 ; drj = 0
                if k == 0: dlk = 1 ; drk = 0
                elif 1 <= k <= 3: dlk = 1 ; drk = 1
                elif k == 4: dlk = 1 ; drk = 0
                doubleLeft = [dli,dlj,dlk]
                doubleRight = [dri,drj,drk]
                if (i,j,k)!=(2,2,2) :
                    (ni,nj,nk,rio,rjo,rko,hio,hjo,hko) = G.cartr2((Px,Py,Pz), HC, (Rx[i],Ry[j],Rz[k]), (Qx,Qy,Qz), doubleLeft, doubleRight, skeleton=True)
                    z = Internal.newZone('Zone', zsize=[[ni,ni-1,0], [nj,nj-1,0], [nk,nk-1,0]], ztype='Structured')
                    n = Internal.newGridCoordinates(parent=z)
                    Internal.newDataArray('CoordinateX', value=None, parent=n)
                    z[0] = 'cart%d-%d-%d'%(i,j,k)
                    if 1<=i<=3 and 1<=j<=3 and 1<=k<=3:
                        coreL.append('cart%d-%d-%d'%(i,j,k))
                    if i > 0 and (i,j,k)!=(3,2,2) :
                        C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i-1,j,k))
                    if i < 4 and (i,j,k)!=(1,2,2):
                        C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i+1,j,k))
                    if j > 0 and (i,j,k)!=(2,3,2):
                        C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z,'jmax', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j-1,k))
                    if j < 4 and (i,j,k)!=(2,1,2):
                        C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j+1,k))
                    if k > 0 and dim == 3 and (i,j,k)!=(2,2,3):
                        C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k-1))
                    if k < 4 and dim == 3 and (i,j,k)!=(2,2,1):
                        C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                        bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                        Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k+1))

                    if i == 0: hio=hio*rio**(ni-3); rio=1./rio; Px=Qx; doubleLeft[0] = 0; doubleRight[0] = 1 
                    if j == 0: hjo=hjo*rjo**(nj-3); rjo=1./rjo; Py=Qy; doubleLeft[1] = 0; doubleRight[1] = 1
                    if k == 0: hko=hko*rko**(nk-3); rko=1./rko; Pz=Qz; doubleLeft[2] = 0; doubleRight[2] = 1
                    data[z[0]] = [(Px,Py,Pz), (hio,hjo,hko), (rio,rjo,rko), (ni,nj,nk), doubleLeft, doubleRight]

                a[i+5*j+25*k] = z

    if dim == 2: # clean list for 2D case
        out = []
        for i in a:
            if i != 0: out.append(i)
        a = out

    t = C.newPyTree(['CARTESIAN','FLEX'])
    core = []
    for z in a:
        if z[0] in coreL: 
            core.append(z)
    for zone in core : 
        t[2][1][2].append(zone)
        a.remove(zone)
    t[2][2][2] += a
    # correction des fenetres max des BCs
    for z in Internal.getZones(t):
        bcs = Internal.getNodesByType(z, 'GridConnectivity1to1_t')
        for BC in bcs:
            PtRangeDonor = Internal.getNodeFromName1(BC, 'PointRangeDonor')[1]
            donorName = Internal.getValue(BC)
            zd = Internal.getNodeFromName2(t, donorName)
            dimz = Internal.getZoneDim(zd)
            imaxDonor = dimz[1]; jmaxDonor = dimz[2]; kmaxDonor = dimz[3]
            if PtRangeDonor[0,0] > 1 and PtRangeDonor[0,0] == PtRangeDonor[0,1]:
                PtRangeDonor[0,0] = imaxDonor
                PtRangeDonor[0,1] = imaxDonor
            if PtRangeDonor[1,0] > 1 and PtRangeDonor[1,0] == PtRangeDonor[1,1]:
                PtRangeDonor[1,0] = jmaxDonor
                PtRangeDonor[1,1] = jmaxDonor
            if PtRangeDonor[2,0] > 1 and PtRangeDonor[2,0] == PtRangeDonor[2,1]:
                PtRangeDonor[2,0] = kmaxDonor
                PtRangeDonor[2,1] = kmaxDonor

    # SplitNParts on core
    if size is None: size = 1
    if rank is None: rank = 0    
    b = Internal.getNodeFromName(t, 'CARTESIAN')
    T._splitNParts(b, N=size, topTree=t)
    D2._distribute(b, NProc=size, algorithm='fast')

    # SplitSize + ressource : distribue en meme temps
    b = Internal.getNodeFromName(t, 'FLEX')
    T._splitSize(b, R=size, topTree=t)
    #D2._distribute(t2, NProc=size, algorithm='fast') # deja fait par splitSize
    D2.printStats(b)
    # Generation reelle
    bases = Internal.getBases(t)
    for b in bases:
        for c in range(len(b[2])):
            z = b[2][c]

            if z[3] == 'Zone_t' and Cmpi.getProc(z) == rank:
                if z[0] in data: # bloc non splitte
                    # print(z[0], 'bloc non splite', flush=True)
                    d = data[z[0]]
                    zn = G.cartr1(d[0], d[1], d[2], d[3], d[4], d[5])
                else:
                    # print(z[0],'bloc splitte', flush=True)
                    source, dest = Internal.getLoc2Glob(z)
                    # print('source', source, flush=True)
                    # print('dest', dest, flush=True)
                    d = data[source]

                    P = d[0]; H = d[1]; R = d[2] ; N = d[3] ; dL = d[4] ; dR = d[5]
                    i1 = dest[0]-1; j1 = dest[2]-1; k1 = dest[4]-1
                    i2 = dest[1]-1; j2 = dest[3]-1; k2 = dest[5]-1
                    doubleLeft=[0,0,0]; doubleRight=[0,0,0]
                    if dL[0] == 1 and i1 == 0: doubleLeft[0] = 1
                    if dR[0] == 1 and i2 == N[0]-1: doubleRight[0] = 1 
                    if dL[1] == 1 and j1 == 0: doubleLeft[1] = 1
                    if dR[1] == 1 and j2 == N[1]-1: doubleRight[1] = 1 
                    if dL[2] == 1 and k1 == 0: doubleLeft[2] = 1
                    if dR[2] == 1 and k2 == N[2]-1: doubleRight[2] = 1 

                    if doubleLeft[0] == 1 and doubleRight[0] == 0:
                        if R[0] == 1.: ratiox = i1; Hx = H[0]
                        else: ratiox = (R[0]**(i1) - 1. )/(R[0]-1.); Hx= H[0]
                    elif doubleLeft[0] == 0 and doubleRight[0] ==0:
                        if R[0] == 1.: ratiox = i1; Hx = H[0]
                        else:
                            if i1 == 0: ratiox = (R[0]**i1-1.)/(R[0]-1.); Hx = H[0]*R[0]**i1
                            else: ratiox = 1. + (R[0]**(i1-1)-1.)/(R[0]-1.); Hx = H[0]*R[0]**(i1-1)
                    elif doubleLeft[0] == 0 and doubleRight[0] == 1:
                        if R[0] == 1: ratiox=i1; Hx = H[0]
                        else: 
                            if i1 == 0: ratiox = (R[0]**i1 - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**i1
                            else: ratiox = (R[0]**(i1) - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**(i1)                                        
                    elif doubleLeft[0] == 1 and doubleRight[0] == 1:
                        if R[0] == 1: ratiox=i1; Hx = H[0]
                        else: 
                            if i1 == 0: ratiox = (R[0]**i1 - 1.)/(R[0] - 1.); Hx = H[0]
                            else: ratiox =  1.+(R[0]**(i1-1) - 1.)/(R[0] - 1.); Hx = H[0]*R[0]**(i1-1)

                    if doubleLeft[1] == 1 and doubleRight[1] == 0:
                        if R[1] == 1.: ratioy = j1; Hy = H[1]
                        else: ratioy = (R[1]**(j1) - 1.)/(R[1]-1.); Hy= H[1] 
                    elif doubleLeft[1] == 0 and doubleRight[1] == 0:
                        if R[1] == 1.: ratioy = j1; Hy = H[1]
                        else:
                            if j1 == 0: ratioy = (R[1]**j1-1.)/(R[1]-1.); Hy = H[1]*R[1]**j1 
                            else: ratioy = 1. + (R[1]**(j1-1.)-1.)/(R[1]-1.); Hy = H[1]*R[1]**(j1 -1)
                    elif doubleLeft[1] == 0 and doubleRight[1] == 1:
                        if R[1] == 1: ratioy=j1; Hy = H[1]
                        else: 
                            if j1 == 0: ratioy = (R[1]**j1 - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**j1
                            else: ratioy = (R[1]**(j1) - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**(j1) 
                    elif doubleLeft[1] == 1 and doubleRight[1] ==1:
                        if R[1] == 1: ratioy=j1 ; Hy = H[1]
                        else: 
                            if j1 == 0: ratioy = (R[1]**j1 - 1.)/(R[1] - 1.); Hy = H[1]
                            else: ratioy =  1.+(R[1]**(j1-1) - 1.)/(R[1] - 1.); Hy = H[1]*R[1]**(j1-1)

                    if doubleLeft[2] == 1 and doubleRight[2] == 0:
                        if R[2] == 1.: ratioz = k1 ; Hz = H[2]
                        else:
                            if k1 == 0: ratioz = (R[2]**(k1) - 1.)/(R[2]-1.); Hz= H[2]

                    elif doubleLeft[2] == 0 and doubleRight[2] == 0:
                        if R[2] == 1.: ratioz = k1; Hz = H[2]
                        else:
                            if k1 == 0: ratioz = (R[2]**k1-1.)/(R[2]-1.); Hz = H[2]*R[2]**k1                    
                            else: ratioz = 1. + (R[2]**(k1-1) -1.)/(R[2]-1.); Hz = H[2]*R[2]**(k1 - 1) 
                    elif doubleLeft[2] == 0 and doubleRight[2] ==1:
                        if R[2] == 1: ratioz=k1; Hz = H[2]
                        else: 
                            if k1 == 0: ratioz = (R[2]**k1 - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**k1
                            else: ratioz =  (R[2]**(k1) - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**(k1)     

                    elif doubleLeft[2] == 1 and doubleRight[2] == 1:
                        if R[2] == 1: ratioz=k1; Hz = H[2]
                        else: 
                            if k1 == 0: ratioz = (R[2]**k1 - 1.)/(R[2] - 1.); Hz = H[2]
                            else: ratioz =  1.+(R[2]**(k1-1) - 1.)/(R[2] - 1.); Hz = H[2]*R[2]**(k1-1)

                    Px = P[0] + ratiox*H[0]
                    Py = P[1] + ratioy*H[1]
                    Pz = P[2] + ratioz*H[2]
                    Rx = R[0]; Ry = R[1]; Rz = R[2]
                    N = (i2-i1+1,j2-j1+1,k2-k1+1)
                    zn = G.cartr1((Px,Py,Pz), (Hx,Hy,Hz), (Rx,Ry,Rz), N, doubleLeft, doubleRight)

                zn[0] = z[0]
                D2._addProcNode(zn, rank)
                n = Internal.getNodesFromName(z, 'ZoneBC')
                zn[2] += n
                n = Internal.getNodesFromName(z, 'ZoneGridConnectivity')
                zn[2] += n
                b[2][c] = zn

    Cmpi._convert2PartialTree(t)
    return t