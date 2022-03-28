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
                    z = G.cart(Xp, H, Np); z[0] = 'cart%d.%d.%d'%(i,j,k)
                    if rank is not None:
                        Cmpi._setProc(z, rank)
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
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i-1,j,k))
                        if i < Nb[0]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i+1,j,k))
                        if j > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z, 'jmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j-1,k))
                        if j < Nb[1]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j+1,k))
                        if k > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j,k-1))
                        if k < Nb[2]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j,k+1))
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

def cartRx2(XC0, XC1, HC, XF0, XF1, R, rank=None, size=None):
    """Create a set of regular and geometric cartesian grids."""

    # Ajuste HC pour core match
    Nx = int((XC1[0]-XC0[0])/HC[0])+1
    Ny = int((XC1[1]-XC0[1])/HC[1])+1
    Nz = int((XC1[2]-XC0[2])/HC[2])+1
    HCx = (XC1[0]-XC0[0])/(Nx-1.)
    HCy = (XC1[1]-XC0[1])/(Ny-1.)
    HCz = (XC1[2]-XC0[2])/(Nz-1.)
    HC = (HCx,HCy,HCz)

    L0x = XC0[0]-XF0[0]
    L1x = XC1[0]-XC0[0]
    L2x = XF1[0]-XC1[0]
    L0y = XC0[1]-XF0[1]
    L1y = XC1[1]-XC0[1]
    L2y = XF1[1]-XC1[1]
    L0z = XC0[2]-XF0[2]
    L1z = XC1[2]-XC0[2]
    L2z = XF1[2]-XC1[2]

    a = [0]* (3*3*3)
    X0x = [XC0[0], XC0[0], XC0[0]+L1x]
    X0y = [XC0[1], XC0[1], XC0[1]+L1y]
    X0z = [XC0[2], XC0[2], XC0[2]+L1z]
    X1x = [XC0[0]-L0x, XC0[0]+L1x, XC0[0]+L1x+L2x]
    X1y = [XC0[1]-L0y, XC0[1]+L1y, XC0[1]+L1y+L2y]
    X1z = [XC0[2]-L0z, XC0[2]+L1z, XC0[2]+L1z+L2z]
    Rx = [R[0],1.,R[0]]
    Ry = [R[1],1.,R[1]]
    Rz = [R[2],1.,R[2]]

    # squelette
    data = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Px = X0x[i]; Py = X0y[j]; Pz = X0z[k]
                Qx = X1x[i]; Qy = X1y[j]; Qz = X1z[k]
                (ni,nj,nk,rio,rjo,rko,hio,hjo,hko) = G.cartr2((Px,Py,Pz), HC, (Rx[i],Ry[j],Rz[k]), (Qx,Qy,Qz), skeleton=True)
                z = Internal.newZone('Zone', zsize=[[ni,ni-1,0], [nj,nj-1,0], [nk,nk-1,0]], ztype='Structured')
                n = Internal.newGridCoordinates(parent=z)
                Internal.newDataArray('CoordinateX', value=None, parent=n)
                z[0] = 'cart%d-%d-%d'%(i,j,k)
                
                if i > 0:
                    C._addBC2Zone(z,'match','BCMatch','imin',z,'imin',[1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i-1,j,k))
                if i < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i+1,j,k))
                if j > 0:
                    C._addBC2Zone(z,'match','BCMatch','jmin',z,'jmin',[1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j-1,k))
                if j < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j+1,k))
                if k > 0:
                    C._addBC2Zone(z,'match','BCMatch','kmin',z,'kmin',[1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k-1))
                if k < 2:
                    C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                    bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                    Internal._setValue(bcs[-1], 'cart%d-%d-%d'%(i,j,k+1))
           
                if i == 0: iv = 1; hio = -hio
                else: iv = 0
                if j == 0: jv = 1; hjo = -hjo
                else: jv = 0
                if k == 0: kv = 1; hko = -hko
                else: kv = 0

                data[z[0]] = [(Px,Py,Pz), (hio,hjo,hko), (rio,rjo,rko), (ni,nj,nk), (iv,jv,kv)]

                a[i+3*j+9*k] = z

    t = C.newPyTree(['Base',a])

    # SplitNParts on skeleton tree
    core = Internal.getNodeFromName2(t, 'cart1-1-1')
    t1 = C.newPyTree(['CARTESIAN',core])
    t1 = T.splitNParts(t1, N=size)
    D2._distribute(t1, NProc=size, algorithm='fast')

    a.remove(core)
    t2 = C.newPyTree(['FLEX',a])
    t2 = T.splitSize(t2, R=size)
    D2._distribute(t2, NProc=size, algorithm='fast')

    t = Internal.merge([t1,t2])
    
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
                    P = d[0]
                    H = d[1]
                    R = d[2]
                    N = d[3]
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
            
                (iv,jv,kv) = d[4]
                if iv == 1: T._reorder(zn, (-1,2,3))
                if jv == 1: T._reorder(zn, (1,-2,3))
                if kv == 1: T._reorder(zn, (1,2,-3))
                zn[0] = z[0]
                D2._addProcNode(zn, rank)
                n = Internal.getNodesFromName(z, 'ZoneBC')
                zn[2] += n
                n = Internal.getNodesFromName(z, 'ZoneGridConnectivity')
                zn[2] += n
                b[2][c] = zn
            
    Cmpi._convert2PartialTree(t)
    return t
