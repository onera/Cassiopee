# - Rotor/propeller specific post-processing -
import Converter.PyTree as C
from . import ExtraVariables2 as PE
import Converter.Internal as Internal
import Transform.PyTree as T
import Generator.PyTree as G
from . import PyTree as P
from . import Mpi as Pmpi
import Converter.Mpi as Cmpi
import RigidMotion.PyTree as R
import numpy
import math

# Helpers
# Echange deux noeuds
def switch(node1, node2):
    temp = node1[1]
    node1[1] = node2[1]
    node2[1] = temp
    return node1, node2

def ProdVect(U,V):  
    Wx = U[1]*V[2]-U[2]*V[1]
    Wy = U[2]*V[0]-U[0]*V[2]
    Wz = U[0]*V[1]-U[1]*V[0]
    return Wx,Wy,Wz

#=========================================================================
# Detection des bords de fuite, attaque et quart de corde 
# Creation du repere local et calcul de la matrice de passage
# IN: a: slice
#     r: rayon de la slice
# OUT: PLE: coord du bord d'attaque repere fixe
#      PTE: coord du bord de fuite repere fixe
#      PF: coord du foyer repere fixe
#      Glob2Loc: Matrice de passage du repere global/fixe au repere local
#=========================================================================
def detect_BA_BF_QC(a, r):

    # Obtention du Leading Edge et Trailing Edge avec xc sur pale non deformee
    xcp = Internal.getNodeFromName2(a, 'xc')[1]
    #iLE = list(xcp).index(C.getMinValue(a,'xc'))
    #iTE = list(xcp).index(C.getMaxValue(a,'xc'))
    iLE = numpy.argmin(xcp)
    iTE = numpy.argmax(xcp)

    # Replace la pale avec son mouvement 
    gridx = Internal.getNodeFromName2(a, 'CoordinateX')
    gridy = Internal.getNodeFromName2(a, 'CoordinateY')
    gridz = Internal.getNodeFromName2(a, 'CoordinateZ')
    xabs = Internal.getNodeFromName2(a, 'Xabs')
    yabs = Internal.getNodeFromName2(a, 'Yabs')
    zabs = Internal.getNodeFromName2(a, 'Zabs')
    switch(gridx, xabs)
    switch(gridy, yabs)
    switch(gridz, zabs)

    # Coordonnees repere absolu/fixe de la pale deformee
    PTE = C.getValue(a, 'GridCoordinates', iTE)
    PLE = C.getValue(a, 'GridCoordinates', iLE)
    PF = [(PTE[i]-PLE[i])*0.25 + PLE[i] for i in range(0,3)]

    #===== Repere local
    # Vx selon la corde, positif de LE vers TE
    Vx = [PTE[i]-PLE[i] for i in range(3)]
    norme = numpy.sqrt(Vx[0]**2+Vx[1]**2+Vx[2]**2)
    Vx = [Vx[i]/norme for i in range(3)]

    # Vy normal au plan de la slice et pointe vers le bout de la pale

    #mesh = G.tetraMesher(a) # optimise par les lignes suivantes
    p = G.fittingPlaster(a)
    mesh = G.gapfixer(a, p, refine=0)

    G._getNormalMap(mesh)
    #C.convertPyTree2File(mesh, 'mesh.cgns')
    (sxmean, symean, szmean) = (-C.getMeanValue(mesh, 'centers:sx'), -C.getMeanValue(mesh, 'centers:sy'), -C.getMeanValue(mesh, 'centers:sz'))
    norme = numpy.sqrt(sxmean**2+symean**2+szmean**2)
    Vy = (sxmean/norme, symean/norme, szmean/norme)

    # Vz pointe vers le haut 
    Vz = ProdVect(Vx, Vy)

    #===== Matrices de passage 
    MatPass = numpy.array([[Vx[i],Vy[i],Vz[i]] for i in range(3)]).real
    #print(MatPass)
    Glob2Loc = numpy.linalg.inv(MatPass)
    #print(Glob2Loc)

    return PLE,PTE,PF,Glob2Loc

#=========================================================================
# IN: t: volume computation tree
# Extrait la vitesse et calcul la vorticite, le critere Q, ... pour les isos
#=========================================================================
def computeVelocityRelated(t, VorticityMagnitude=False, QCriterion=False, VelocityMagnitude=False):
    tp = C.extractVars(t, ['centers:VelocityX','centers:VelocityY','centers:VelocityZ','centers:cellN'], keepOldNodes=False)
    # calcul de la vorticite
    if VorticityMagnitude: PE._computeVorticityMagnitude2(tp, ghostCells=True)
    # calcul du critere Q
    if QCriterion: PE._computeQCriterion2(tp, ghostCells=True)
    # calcul du module de la vitess
    if VelocityMagnitude: P._computeVariables2(t, ['centers:VelocityMagnitude'])
    tp = Internal.rmGhostCells(tp, tp, 2)
    # extra overlap cleaning for blades
    #for name in ['Blade7A_00', 'Blade7A_01', 'Blade7A_02', 'Blade7A_03']:
    #    b = Internal.getNodeFromName1(tp, name)
    #    for c in range(len(b[2])):
    #        if b[2][c][3] == 'Zone_t':
    #            b[2][c] = T.subzone(b[2][c], (1,1,1), (-1,-1,-3))
    tp = C.node2Center(tp)
    return tp

#====================================================================
# Export a (psi,rad) accumulator storing vars to a zone
def exportAccumulatorPerPsi(accumulator, psi=0., vars=['F1','F2']):
    """Export accumulator (psi,rad) in a zone for a given psi."""
    accumulator = Cmpi.allgatherDict(accumulator)

    radius = []
    for k in accumulator.keys(): radius.append(k[1])

    radius = sorted(set(radius))
    # Create zone
    z = G.cart((0,0,0), (1,0,0), (len(radius),1,1))
    C._addVars(z, 'radius')
    C._addVars(z, vars)
    qx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    qr = Internal.getNodeFromName2(z, 'radius')[1]
        
    for r, rad in enumerate(radius):
        acu = accumulator[(psi,rad)]
        qx[r] = rad
        qr[r] = rad
        for c, v in enumerate(vars):
            q = Internal.getNodeFromName2(z, v)[1]
            q[r] = acu[c]
    return z

#====================================================================
# Export a (psi,rad) accumulator storing vars to a zone
def exportAccumulatorPerRadius(accumulator, rad=0., vars=['F1','F2']):
    """Export accumulator (psi,rad) in a zone for a given radius."""

    accumulator = Cmpi.allgatherDict(accumulator)

    psis = []
    for k in accumulator.keys(): psis.append(k[0])

    psis = sorted(set(psis))
    # Create zone
    z = G.cart((0,0,0), (1,0,0), (len(psis),1,1))
    C._addVars(z, 'psi')
    C._addVars(z, vars)
    qx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    qp = Internal.getNodeFromName2(z, 'psi')[1]
        
    for p, psi in enumerate(psis):
        acu = accumulator[(psi,rad)]
        qx[p] = psi
        qp[p] = psi
        for c, v in enumerate(vars):
            q = Internal.getNodeFromName2(z, v)[1]
            q[p] = acu[c]
    return z

#====================================================================
# export un accumulateur (psi,rad) sous forme d'une map psi/radius
# IN: vars: the variables stored in accumulator
#====================================================================
def exportAccumulatorMap(accumulator, vars=['Fx','Fy','Fz']):
    """Export accumulator (psi,rad) in a map zone."""
    accumulator = Cmpi.allgatherDict(accumulator)

    psis = []; radius = []
    for k in accumulator.keys():
        psis.append(k[0]) # psi
        radius.append(k[1]) # rad
    psis = sorted(set(psis))
    radius = sorted(set(radius))

    # Create zone
    z = G.cart((0,0,0), (1,0,0), (len(psis),len(radius),1))
    qx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    qy = Internal.getNodeFromName2(z, 'CoordinateY')[1]

    C._addVars(z, vars)
    if len(psis) == 1:
        psirad = psis[0]*math.pi/180.
        for r, rad in enumerate(radius):
                acu = accumulator[(psis[0],rad)]
                qx[r] = rad * math.cos(psirad)
                qy[r] = rad * math.sin(psirad)
                for c, v in enumerate(vars):
                    q = Internal.getNodeFromName2(z, v)[1]
                    q[r] = acu[c]
    elif len(radius) == 1:
        for p, psi in enumerate(psis):
            psirad = psi*math.pi/180.
            acu = accumulator[(psi,radius[0])]
            qx[p] = radius[0] * math.cos(psirad)
            qy[p] = radius[0] * math.sin(psirad)
            for c, v in enumerate(vars):
                q = Internal.getNodeFromName2(z, v)[1]
                q[p] = acu[c]
    else:
        for p, psi in enumerate(psis):
            psirad = psi*math.pi/180.
            for r, rad in enumerate(radius):
                acu = accumulator[(psi,rad)]
                qx[p,r] = rad * math.cos(psirad)
                qy[p,r] = rad * math.sin(psirad)
                for c, v in enumerate(vars):
                    q = Internal.getNodeFromName2(z, v)[1]
                    q[p,r] = acu[c]
    return z

#========================================================================
# Zb: dans le repere vent
# IN: teff
# IN: RoInf: far field density (kg/m3)
# IN: ASOUND: sound speed (OMEGA * R = Mtip * ASOUND)
# IN: Mtip: tip mach number
# IN: AR: rotor radius (m)
# IN: SIGMA: rotor solidity = (Nb*c) / (pi*AR)
# IN: relativeShaft: si le repere du maillage n'est pas le repere vent
# OUT: rotor traction
#========================================================================
def computeZb(teff, psi, RoInf, ASOUND, Mtip, AR, SIGMA, 
              relativeShaft=0., accumulatorZb=None):
    """Compute Zb."""
    PE._extractShearStress(teff)
    PE._extractForce(teff, withPInf=RoInf*ASOUND**2/1.4)
    if relativeShaft != 0.:
        teff = T.rotate(teff, (0,0,0), (0,1,0), relativeShaft, vectors=[['centers:sx','centers:sy','centers:sz'],['centers:Fx','centers:Fy','centers:Fz']])
    xb = Pmpi.integ(teff, 'centers:Fx')[0]
    yb = Pmpi.integ(teff, 'centers:Fy')[0]
    zb = Pmpi.integ(teff, 'centers:Fz')[0]
    S = math.pi * AR**2
    adim = 0.5*RoInf*S*SIGMA*(Mtip*ASOUND)**2
    xb = 100. * xb / adim
    yb = 100. * yb / adim
    zb = 100. * zb / adim
    if accumulatorZb is not None: accumulatorZb[(psi,0)] = [xb,yb,zb]
    return [xb,yb,zb]

#========================================================================
# Thrust: dans le repere shaft (ortho au rotor)
# IN: teff
# IN: PINF: infinite pressure
# IN: relativeShaft: si le repere du maillage n'est pas le repere shaft (ortho rotor)
# OUT: Thrust en Newton et Torque en N.m
#========================================================================
def computeThrustAndTorque(teff, psi, PInf, center=(0,0,0), relativeShaft=0., accumulatorThrust=None):
    """Compute thrust."""
    PE._extractShearStress(teff)
    PE._extractForce(teff, withPInf=PInf)
    if relativeShaft != 0.:
        teff = T.rotate(teff, (0,0,0), (0,1,0), relativeShaft, vectors=[['centers:sx','centers:sy','centers:sz'],['centers:Fx','centers:Fy','centers:Fz']])
    thrustx = Pmpi.integ(teff, 'centers:Fx')[0]
    thrusty = Pmpi.integ(teff, 'centers:Fy')[0]
    thrustz = Pmpi.integ(teff, 'centers:Fz')[0]
    torque = Pmpi.integMoment(teff, center=center, vector=['centers:Fx','centers:Fy','centers:Fz'])
    if accumulatorThrust is not None: accumulatorThrust[(psi,0)] = [thrustx,thrusty,thrustz,torque[0],torque[1],torque[2]]
    return [thrustx,thrusty,thrustz],[torque[0],torque[1],torque[2]]

#========================================================================
# IN: teff avec frictionXYZ
# OUT: streamlines sur la surface
#========================================================================
def frictionLines(teff):
    """Compute friction lines on a surface."""
    b = Internal.getNodeFromName1(teff, 'Blade7A_00')

    # Points de depart des streamLines
    points = []
    z = Internal.getNodeFromName1(b, 'Blade7A_00_15_BCWall0')
    for i in range(20):
        Pt = C.getValue(z, 'GridCoordinates', (i,35,1))
        Pt = (Pt[0]+1.e-6,Pt[1],Pt[2])
        points.append(Pt)
    z = Internal.getNodeFromName1(b, 'Blade7A_00_12_BCWall.50')
    for i in range(20):
        Pt = C.getValue(z, 'GridCoordinates', (i,35,1))
        Pt = (Pt[0]+1.e-6,Pt[1],Pt[2])
        points.append(Pt)    
    s1 = P.streamLine2(b, points, vector=['centers:frictionX','centers:frictionY','centers:frictionZ'])
    return s1

#========================================================================
# extrait les slices en parallele et calcule les variables de post traitement 
# IN: teff
# IN: bladeName: name of base of blade
# IN: radius: liste des radius a extraire
# IN: delta: interne, largeur des bandelettes
# IN: relativeShaft: si localFrame=False et le repere maillage n'est pas le repere vent
# OUT: dictionnaire des slices, slices[rad] est une zone
# OUT: accumule les valeurs de CnM2 dans un dictionnaire (psi,rad)
#========================================================================
def extractSlices(teff, bladeName, psi, radius, 
                  RoInf, PInf, ASOUND, Mtip, AR, CHORD, MU,
                  accumulatorSlices=None,
                  accumulatorCnM2=None, accumulatorCmM2=None, 
                  relativeShaft=0., localFrame=True, delta=0.05):
    """Extract slices on blade and compute Kp,Cf,CnM2,CmM2."""
    b = Internal.getNodeFromName1(teff, bladeName)
    if b is not None:
        # Complete blade
        PE._extractShearStress(b)
        PE._extractFrictionVector(b)
        PE._extractFrictionMagnitude(b)
        PE._extractForce(b, withPInf=PInf)
        G._getNormalMap(b)
        C._normalize(b, ['centers:sx','centers:sy','centers:sz'])
        if relativeShaft != 0.:
            b = T.rotate(b, (0,0,0), (0,1,0), relativeShaft, vectors=[['centers:sx','centers:sy','centers:sz'],['centers:frictionX','centers:frictionY','centers:frictionZ'],['centers:Fx','centers:Fy','centers:Fz']]) 

        # switch GridInit and GridCoordinates and save gridCoordinates in field
        bp = Internal.copyRef(b)
        zones = Internal.getZones(bp)
        for z in zones:
            ga = Internal.getNodeFromName1(z, 'GridCoordinates')
            gi = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
            fs = Internal.createUniqueChild(z, 'FlowSolution', 'FlowSolution_t')
            xc = Internal.getNodeFromName1(ga, 'CoordinateX')
            xc = Internal.copyRef(xc)
            xc[0] = 'Xabs'
            fs[2].append(xc)
            yc = Internal.getNodeFromName1(ga, 'CoordinateY')
            yc = Internal.copyRef(yc)
            yc[0] = 'Yabs'
            fs[2].append(yc)
            zc = Internal.getNodeFromName1(ga, 'CoordinateZ')
            zc = Internal.copyRef(zc)
            zc[0] = 'Zabs'
            fs[2].append(zc)
            R._switchGridAndGridInit(z)
    else: bp = []

    # Bandeletage sur radius et gather sur des procs differents
    zones = Internal.getZones(bp)
    proc = 0
    allbandes = {}
    for rad in radius:
        bandes = []
        for z in zones:
            C._addVars(z, 'tag')
            ga = Internal.getNodeFromName1(z, 'GridCoordinates')
            xc = Internal.getNodeFromName1(ga, 'CoordinateX')
            xc[1] = xc[1].ravel('k')
            tag = Internal.getNodeFromName2(z, 'tag')
            tag[1] = tag[1].ravel('k')
            tag[1][:] = (xc[1][:] < rad+delta) & (xc[1][:] > rad-delta)
            sel = P.selectCells2(z, 'tag')
            Internal._rmNodesFromName(sel, 'GridCoordinates#Init')
            Internal._rmNodesFromName(sel, '.Solver#ownData')
            Internal._rmNodesFromName(sel, 'tag')
            if C.getNPts(sel) > 0: bandes.append(sel)
        #print('sending to ',proc,bandes)
        allbandes[rad] = Cmpi.gatherZones(bandes, root=proc)
        proc += 1
        proc = proc%Cmpi.size
    
    # Chaque proc a sa ou ses bandelettes (liste de zones)
    #print(Cmpi.rank, allbandes.keys())
    #for rad in allbandes:
    #    b = allbandes[rad]
    #    if b != []:
    #        print('Proc ',Cmpi.rank, 'writes bandes ',rad)
    #        C.convertPyTree2File(allbandes[rad], 'bandes%f.cgns'%rad)
    
    slices = {}
    for rad in allbandes:
        #Internal.printTree(allbandes[rad])
        b = allbandes[rad] # liste de zones
        #print(rad,b)
        if b != []:
            bandes = C.convertArray2Hexa(b)
            bandes = T.join(bandes)
            bandes = C.center2Node(bandes, 'FlowSolution#Centers')
            C._rmVars(bandes, 'FlowSolution#Centers')
            slice = P.isoSurfMC(bandes, 'CoordinateX', rad)
            slices[rad] = slice[0]

    CnM2All = []; CmM2All = []
    for rad in slices:
        iso = slices[rad]

        #no = radius.index(rad) # position du rayon de la slice actuelle dans la liste complete des rayons
        #print('Proc ',Cmpi.rank, 'writes slice ',rad)
        #C.convertPyTree2File(iso, 'slice%f.cgns'%rad)

        # x/c
        ymin = C.getMinValue(iso, 'CoordinateY')
        ymax = C.getMaxValue(iso, 'CoordinateY')
        C._initVars(iso, '{xc}= 1.-({CoordinateY}-%20.16g)/(%20.16g-%20.16g)'%(ymin,ymax,ymin))
        
        # Kp
        adimKp = 0.5*RoInf*(rad*Mtip*ASOUND/AR+MU*Mtip*ASOUND*math.sin(psi/180.*math.pi))**2
        C._initVars(iso, '{Kp}= ({Pressure}-%20.16g)/ %20.16g'%(PInf,adimKp))
        # Cf
        C._initVars(iso, '{Cf}= {frictionMagnitude}/ %20.16g'%adimKp)

        # Passage struct
        iso = C.convertBAR2Struct(iso)
        iso[0] = '%f'%rad        
        #C.convertPyTree2File(iso, 'slice%f.cgns'%rad)

        # Detection et repere local
        (PLE,PTE,PF,Glob2Loc) = detect_BA_BF_QC(iso, rad)

        # CnM2
        adimCnM2 = 0.5*RoInf*ASOUND**2*CHORD
        CnM2x = P.integ(iso, var='Fx')[0]
        CnM2y = P.integ(iso, var='Fy')[0]
        CnM2z = P.integ(iso, var='Fz')[0]
        CnM2x = CnM2x/adimCnM2
        CnM2y = CnM2y/adimCnM2
        CnM2z = CnM2z/adimCnM2
        Cn = (CnM2x, CnM2y, CnM2z)
        if localFrame: Cn = numpy.dot(Glob2Loc, Cn)

        if accumulatorCnM2 is not None: accumulatorCnM2[(psi,rad)] = [Cn[0],Cn[1],Cn[2]]
        CnM2All.append([Cn[0],Cn[1],Cn[2]])

        # CmM2
        adimCmM2 = 0.5*RoInf*ASOUND**2*CHORD
        (CmM2x,CmM2y,CmM2z) = P.integMoment(iso, center=PF, vector=['Fx','Fy','Fz'])
        CmM2x = CmM2x/adimCmM2
        CmM2y = CmM2y/adimCmM2
        CmM2z = CmM2z/adimCmM2
        Cm = (CmM2x, CmM2y, CmM2z)
        if localFrame: Cm = numpy.dot(Glob2Loc, Cm)
        
        if accumulatorCmM2 is not None: accumulatorCmM2[(psi,rad)] = [Cm[0],Cm[1],Cm[2]]
        CmM2All.append([Cm[0],Cm[1],Cm[2]])

    if accumulatorSlices is not None:
        for rad in slices: accumulatorSlices[(psi,rad)] = slices[rad]

    slicesAll = []
    for rad in slices: slicesAll += [slices[rad]]
    return slicesAll, CnM2All, CmM2All

#============================================================
# extract the files for H2T RotorLoads
# IN: teff: stress tree 
# IN: bladeName: nom de la base blade a extraire
# IN: nblade: numero de la pale a extraire (pour les noms de fichiers de sortie)
# IN: it: iteration correspondant a l'extraction de teff (pour les noms de fichiers de sortie)
# IN: relativeShaft: decalage de shaft (si on fait un calcul a plat
# avec le shaft dans la vitesse infinie = -alpha_shaft)
#============================================================
def extractForH2TRotorLoads(teff, bladeName, nblade, it, relativeShaft=0.):
    b = Internal.getNodeFromName1(teff, bladeName)
    bp = Internal.copyRef(b)
    
    # extrait le maillage de reference 
    R._switchGridAndGridInit(bp)
    for c, z in enumerate(Internal.getZones(bp)):
        Internal._rmNodesFromName(z, 'FlowSolution#Centers')
        Internal._rmNodesFromName(z, 'GridCoordinates#Init')
        C.convertPyTree2File(z, 'BLADESURF_MESHREF/%s_SurfRef%03d.tp'%(bladeName,c+1), 'fmt_tp')

    # sortie psta
    for c, z in enumerate(Internal.getZones(bp)):
        #C._initVars(z, '{centers:Pressure}={centers:Pressure}*%16.12g'%adimP)
        z = C.extractVars(z, ['centers:Pressure'], keepOldNodes=False)
        P._renameVars(z, ['centers:Pressure'], ['centers:p'])
        C.convertPyTree2File(z, 'SURFACES/ExtBlade%04dpsta_%04d.tp%04d'%(nblade,c+1,it), 'bin_tp')
        
    # sortie friction vector
    bp = Internal.copyRef(b)
    PE._extractShearStress(bp)
    PE._extractFrictionVector(bp)
    PE._extractFrictionMagnitude(bp)
    G._getNormalMap(bp)

    # Pour retrouver le repere absolu
    T._rotate(bp, (0,0,0), (0,1,0), relativeShaft, vectors=[['centers:sx','centers:sy','centers:sz'],['centers:frictionX','centers:frictionY','centers:frictionZ']]) 

    for c, z in enumerate(Internal.getZones(bp)):
        z = C.extractVars(z, ['centers:sx', 'centers:sy', 'centers:sz', 'centers:frictionX','centers:frictionY','centers:frictionZ','centers:frictionMagnitude'], keepOldNodes=False)
        #C._initVars(z, '{centers:frictionX}={centers:frictionX}*%16.12g'%adimF)
        #C._initVars(z, '{centers:frictionY}={centers:frictionY}*%16.12g'%adimF)
        #C._initVars(z, '{centers:frictionZ}={centers:frictionZ}*%16.12g'%adimF)
        #C._initVars(z, '{centers:frictionMagnitude}={centers:frictionMagnitude}*%16.12g'%adimF)
        P._renameVars(z, ['centers:sx','centers:sy','centers:sz','centers:frictionX','centers:frictionY', 'centers:frictionZ','centers:frictionMagnitude'], ['centers:nx','centers:ny','centers:nz','centers:frictionvectorx','centers:frictionvectory','centers:frictionvectorz','centers:frictionmodulus'])
        C.convertPyTree2File(z, 'SURFACES/ExtBlade%04dvect_%04d.tp%04d'%(nblade,c+1,it), 'bin_tp')

    return None

#=========================================================================
# Extrait des slices a certains radius
# Calcul Kp, Cf, CnM2, CmM2
# Ancienne version
#=========================================================================
""" def extractSlices(b, psib, radius, ROINF, PINF, ASOUND, MTIP, MU, AR, CHORD):
    global CnM2dict, CmM2dict
    slices = []
    CnM2dict[psib] = numpy.empty((3,len(radius)), dtype=numpy.float64)
    CmM2dict[psib] = numpy.empty((3,len(radius)), dtype=numpy.float64)
    R._switchGridAndGridInit(b)
    for no, r in enumerate(radius):
        iso = P.isoSurfMC(b, 'CoordinateX', r)
        if iso == []:
            print('Warning: empty slice.',r)
            continue
        iso = T.join(iso)
        iso[0] = '%20.16g'%r
        # x/c
        ymin = C.getMinValue(iso, 'CoordinateY')
        ymax = C.getMaxValue(iso, 'CoordinateY')
        C._initVars(iso, '{xc}= 1.-({CoordinateY}-%20.16g)/(%20.16g-%20.16g)'%(ymin,ymax,ymin))
        # Kp
        adimKp = 0.5*ROINF*(r*MTIP*ASOUND/AR+MU*MTIP*ASOUND*math.sin(psib/180.*math.pi))**2
        C._initVars(iso, '{Kp}= ({Pressure}-%20.16g)/ %20.16g'%(PINF,adimKp))
        # Cf
        C._initVars(iso, '{Cf}= {frictionMagnitude}/ %20.16g'%adimKp)

        # passage struct
        iso = C.convertBAR2Struct(iso)
        iso[0] = '%20.16g'%r
        
        # CnM2
        adimCnM2 = 0.5*ROINF*ASOUND**2*CHORD
        CnM2x = P.integ(iso, var='Fx')[0]
        CnM2y = P.integ(iso, var='Fy')[0]
        CnM2z = P.integ(iso, var='Fz')[0]
        CnM2x = CnM2x/adimCnM2
        CnM2y = CnM2y/adimCnM2
        CnM2z = CnM2z/adimCnM2
        CnM2dict[psib][0,no] = CnM2x
        CnM2dict[psib][1,no] = CnM2y
        CnM2dict[psib][2,no] = CnM2z
        
        # CmM2
        adimCmM2 = 0.5*ROINF*ASOUND**2*CHORD**2
        #(PLE,PTE,PF) = detect_BA_BF_QC(iso,r)
        (CmM2x,CmM2y,CmM2z) = P.integMoment(iso, center=(r,0.,0.), vector=['Fx','Fy','Fz'])
        CmM2x = CmM2x/adimCmM2
        CmM2y = CmM2y/adimCmM2
        CmM2z = CmM2z/adimCmM2
        CmM2dict[psib][0,no] = CmM2x
        CmM2dict[psib][1,no] = CmM2y
        CmM2dict[psib][2,no] = CmM2z
        # retour
        slices.append(iso)
    R._switchGridAndGridInit(b)
    return slices
 """

#=========================================================================
# Calcul les champs derives sur teff
# IN: b: blade base of teff
#=========================================================================
""" def completeBlade(teff):
    # Calcul des grandeurs supplementaires
    PE._extractShearStress(teff)
    PE._extractFrictionVector(teff)
    PE._extractFrictionMagnitude(teff)
    PE._extractForce(teff)
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx','centers:sy','centers:sz'])
    # Passe de centres en noeuds proprement
    teff = C.center2Node(teff, 'FlowSolution#Centers')
    Internal._rmNodesFromName(teff, 'FlowSolution#Centers')
    return teff
 """