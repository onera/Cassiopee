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
def detectBA_BF_QC(a, r):

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

    radii = []
    for k in accumulator.keys(): radii.append(k[1])

    radii = sorted(set(radii))
    # Create zone
    z = G.cart((0,0,0), (1,0,0), (len(radii),1,1))
    C._addVars(z, 'radius')
    C._addVars(z, vars)
    qx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    qr = Internal.getNodeFromName2(z, 'radius')[1]

    for r, rad in enumerate(radii):
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

    psis = []; radii = []
    for k in accumulator.keys():
        psis.append(k[0]) # psi
        radii.append(k[1]) # rad
    psis = sorted(set(psis))
    radii = sorted(set(radii))

    # Create zone
    z = G.cart((0,0,0), (1,0,0), (len(psis),len(radii),1))
    qx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    qy = Internal.getNodeFromName2(z, 'CoordinateY')[1]

    # C._initVars(z, 'Radius', 0)
    # radii = Internal.getNodeFromName2(z, 'Radius')[1]
    # for r, rad in enumerate(fadii):
    #     radii[r] = rad
    # Internal.getNodeFromName2(z, 'Radius')[1] = radii

    C._addVars(z, vars)
    if len(psis) == 1:
        psirad = psis[0]*math.pi/180.
        for r, rad in enumerate(radii):
            acu = accumulator[(psis[0],rad)]
            qx[r] = rad * math.cos(psirad)
            qy[r] = rad * math.sin(psirad)
            for c, v in enumerate(vars):
                q = Internal.getNodeFromName2(z, v)[1]
                q[r] = acu[c]
    elif len(radii) == 1:
        for p, psi in enumerate(psis):
            psirad = psi*math.pi/180.
            acu = accumulator[(psi,radii[0])]
            qx[p] = radii[0] * math.cos(psirad)
            qy[p] = radii[0] * math.sin(psirad)
            for c, v in enumerate(vars):
                q = Internal.getNodeFromName2(z, v)[1]
                q[p] = acu[c]
    else:
        for p, psi in enumerate(psis):
            psirad = psi*math.pi/180.
            for r, rad in enumerate(radii):
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
    """Compute thrust and torque."""
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
# extract the radius field on a surface, using the rotation axis and center
# IN: teff
# OUT: teff with radius field
#========================================================================
def extractRadius(teff, axis_pnt, axis_vct, loc='node'):
    """Extract the radius field, using the rotation axis and center."""
    def function(x,y,z):
        ux,uy,uz = axis_vct
        cx,cy,cz = axis_pnt
        ax = (y-cy)*uz - (z-cz)*uy
        ay = (x-cx)*uz - (z-cz)*ux
        az = (x-cx)*uy - (y-cy)*ux
        return math.sqrt(ax**2 + ay**2 + az**2)/math.sqrt(ux**2 + uy**2 + uz**2)

    if loc == 'center':
        teff = C.initVars(teff, 'centers:Radius', function, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    else:
        teff = C.initVars(teff, 'Radius', function, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
    return teff

#========================================================================
# extract the theta field on a surface, using the rotation axis and center
# IN: teff
# OUT: teff with theta field
#========================================================================
def extractTheta(teff, axis_pnt, axis_vct, loc='node'):
    """Extract the theta field, using the rotation axis and center."""
    # local functions for 3x3 matrices
    def determinant(m):
        det = 0

        det += m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
        det -= m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
        det += m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0])

        return det

    def comatrix(m):
        cm = [[0,0,0], [0,0,0], [0,0,0]]

        for i in range(3):
            i1 = 0 if i != 0 else 1
            i2 = 2 if i != 2 else 1
            for j in range(3):
                j1 = 0 if j != 0 else 1
                j2 = 2 if j != 2 else 1
                cm[i][j] = (-1)**(i+j)*(m[i1][j1]*m[i2][j2] - m[i1][j2]*m[i2][j1])
        return cm

    def transpose(m):
        tm = [[0,0,0], [0,0,0], [0,0,0]]
        for i in range(3):
            for j in range(3):
                tm[i][j] = m[j][i]
        return tm

    def inverse(m):
        det = determinant(m)
        tm = transpose(comatrix(m))
        for i in range(3):
            for j in range(3):
                tm[i][j] /= det
        return tm
    ###################################

    # create a base from axis_vct
    a1,b1,c1 = axis_vct
    vec1 = (a1,b1,c1)

    # get ortho vec
    if abs(a1) > abs(b1):
        a2,b2,c2 = c1,0,-a1
    else:
        a2,b2,c2 = 0,c1,-b1
    vec2 = (a2,b2,c2)

    # complete the base with cross product vec1 x vec2
    a3 = b1*c2 - b2*c1
    b3 = a1*c2 - a2*c1
    c3 = a1*b2 - a2*b1
    vec3 = (a3,b3,c3)

    # normalize the three vectors
    base = []
    for vec in [vec1, vec2, vec3]:
        a,b,c = vec
        norm = math.sqrt(a*a + b*b + c*c)
        base.append((a/norm, b/norm, c/norm))

    # matrice de passage de R dans base
    P = [[base[0][0], base[1][0], base[2][0]],
         [base[0][1], base[1][1], base[2][1]],
         [base[0][2], base[1][2], base[2][2]]]

    # inverse de P
    invP = inverse(P)

    # C (canonical base) -> Cprime (new base) using
    # C = P.Cprime <=> Cprime = invP.C
    # C = [axis_pnt[0], axis_pnt[1], axis_pnt[2]]

    Cprime = [axis_pnt[0]*invP[0][0] + axis_pnt[1]*invP[0][1] + axis_pnt[2]*invP[0][2],
              axis_pnt[0]*invP[1][0] + axis_pnt[1]*invP[1][1] + axis_pnt[2]*invP[1][2],
              axis_pnt[0]*invP[2][0] + axis_pnt[1]*invP[2][1] + axis_pnt[2]*invP[2][2]]

    def function(x,y,z):
        # X' = invP*X
        # first direction is the rotation axis
        dist2 = x*invP[1][0] + y*invP[1][1] + z*invP[1][2] - Cprime[1]
        dist3 = x*invP[2][0] + y*invP[2][1] + z*invP[2][2] - Cprime[2]

        theta = numpy.arctan2(-dist2,dist3)
        if theta < 0: theta += 2*numpy.pi

        return theta

    if loc == 'center':
        teff = C.initVars(teff, 'centers:Theta', function, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    else:
        teff = C.initVars(teff, 'Theta', function, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
    return teff

#========================================================================
# extrait les slices en parallele et calcule les variables de post traitement
# IN: teff
# IN: bladeName: name of base of blade
# IN: psi: current angle of rotation, used as key in all accumulator dict
# IN: radii: list of radii at which the solution on the blades must be extracted
# IN: RoInf, PInf, ASOUND, Mtip, AR, CHORD, MU: flow information used for adimensioning
# IN: accumulatorSlices: dict of slices used to track their evolution. Keys are (psi,rad)
# IN: accumulatorCnM2: dict of CnM2 used to track their evolution. Keys are (psi,rad)
# IN: accumulatorCmM2: dict of CmM2 used to track their evolution. Keys are (psi,rad)
# IN: adimCnM2: scaling value for CnM2. If adimCnM2=0: computes the value with adimCnM2=0.5*RoInf*ASOUND**2*CHORD
# IN: adimCmM2: scaling value for CmM2. If adimCmM2=0: computes the value with adimCmM2=0.5*RoInf*ASOUND**2*CHORD
# IN: adimKp: scaling value for Kp. If Kp=0: computes the value with Kp=0.5*RoInf*(abs(rad)*Mtip*ASOUND/AR+MU*Mtip*ASOUND*math.sin(psi/180.*math.pi))**2
# IN: relativeShaft: relative shaft angle if the mesh is not in the wind frame
# IN: localFrame: if True, returns CnM2 and CmM2 in relative (blade section) frame
# IN: delta: mean mesh step on blade in the span wise direction
# IN: rotationCenter: center of rotation
# IN: coordDir: axis of rotation ('CoordinateX', 'CooridnateY' or 'CoordinateZ')
# IN: coordSlice: slicing direction ('CoordinateX', 'CoordinateY' or 'CoordinateZ')
# IN: sliceNature: if 'straight', slices the blade(s) in the slicing direction coordSlice. If 'curved', slices at constant radii.

# OUT: dictionnaire des slices, slices[rad] est une zone
# OUT: accumule les valeurs de CnM2 dans un dictionnaire (psi,rad)
#========================================================================
def extractSlices(teff, bladeName, psi, radii,
                  RoInf, PInf, ASOUND, Mtip, AR, CHORD, MU,
                  accumulatorSlices=None,
                  accumulatorCnM2=None, accumulatorCmM2=None,
                  adimCnM2=0, adimCmM2=0, adimKp=0,
                  relativeShaft=0., localFrame=True, delta=0.05, rotationCenter=[0.,0.,0.],
                  coordDir='CoordinateZ', coordSlice='CoordinateX', sliceNature='straight'):
    """Extract slices on blade and compute Kp,Cf,CnM2,CmM2."""
    if coordDir == coordSlice:
        raise ValueError('extractSlices: coordDir and coordSlice are identical.')

    if sliceNature not in ['curved', 'straight']:
        print('Warning: extractSlices: invalid sliceNature name. Default value is used')
        sliceNature = 'straight'

    if 'CoordinateX' not in [coordDir, coordSlice]: coordXSC = 'CoordinateX'
    elif 'CoordinateY' not in [coordDir, coordSlice]: coordXSC = 'CoordinateY'
    else: coordXSC = 'CoordinateX'

    if sliceNature == 'curved':
        if coordDir == 'CoordinateX':
            cy,cz = rotationCenter[1],rotationCenter[2]
            C._initVars(teff,'{Radius}=sqrt(({CoordinateY}-%f)*({CoordinateY}-%f) + ({CoordinateZ}-%f)*({CoordinateZ}-%f))'%(cy,cy,cz,cz))
        elif coordDir == 'CoordinateY':
            cx,cz = rotationCenter[0],rotationCenter[2]
            C._initVars(teff,'{Radius}=sqrt(({CoordinateX}-%f)*({CoordinateX}-%f) + ({CoordinateZ}-%f)*({CoordinateZ}-%f))'%(cx,cx,cz,cz))
        else:
            cx,cy = rotationCenter[0],rotationCenter[1]
            C._initVars(teff,'{Radius}=sqrt(({CoordinateX}-%f)*({CoordinateX}-%f) + ({CoordinateY}-%f)*({CoordinateY}-%f))'%(cx,cx,cy,cy))

    # def arctan3(y,z):
    #     theta = numpy.arctan2(-y,z)
    #     if theta < 0: theta += 2*numpy.pi
    #     return theta
    # C._initVars(teff,'Theta', arctan3, [coordXSC, coordSlice])
    # C._initVars(teff,'{Theta2}={Theta}*%f/%f'%(180., math.pi))

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

        # switch GridInit and GridCoordinates and save GridCoordinates in field
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
    for rad in radii:
        bandes = []
        for z in zones:
            C._addVars(z, 'tag')
            ga = Internal.getNodeFromName1(z, 'GridCoordinates')
            if sliceNature == 'straight':
                xc = Internal.getNodeFromName1(ga, coordSlice)
            else: # sliceNature curved
                xc = Internal.getNodeFromName(z, 'Radius')
            xc[1] = xc[1].ravel('k')
            tag = Internal.getNodeFromName2(z, 'tag')
            tag[1] = tag[1].ravel('k')
            tag[1][:] = (xc[1][:] < rad+delta) & (xc[1][:] > rad-delta)
            if sliceNature == 'curved':
                sel = P.selectCells2(z, 'tag', strict=1)
            else:
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
            if sliceNature == 'straight':
                slice = P.isoSurfMC(bandes, coordSlice, rad)
            else: # sliceNature curved
                slice = P.isoSurfMC(bandes, 'Radius', rad)
            #
            if slice == []: print('Warning: extractSlices: no slice found at position %f'%rad)
            else: slices[rad] = slice[0]

    CnM2All = []; CmM2All = []
    for rad in slices:
        if rad not in slices: continue

        iso = slices[rad]

        #no = radius.index(rad) # position du rayon de la slice actuelle dans la liste complete des rayons
        #print('Proc ',Cmpi.rank, 'writes slice ',rad)
        #C.convertPyTree2File(iso, 'slice%f.cgns'%rad)

        # x/c
        xmin = C.getMinValue(iso, coordXSC)
        xmax = C.getMaxValue(iso, coordXSC)
        if coordXSC == 'CoordinateX':
            C._initVars(iso, '{xc}= 1.-({CoordinateX}-%20.16g)/(%20.16g-%20.16g)'%(xmin,xmax,xmin))
        elif coordXSC == 'CoordinateY':
            C._initVars(iso, '{xc}= 1.-({CoordinateY}-%20.16g)/(%20.16g-%20.16g)'%(xmin,xmax,xmin))
        else:
            C._initVars(iso, '{xc}= 1.-({CoordinateZ}-%20.16g)/(%20.16g-%20.16g)'%(xmin,xmax,xmin))


        # Kp
        if adimKp == 0: adimKp_loc = 0.5*RoInf*(abs(rad)*Mtip*ASOUND/AR + MU*Mtip*ASOUND*math.sin(psi/180.*math.pi))**2
        else: adimKp_loc = adimKp
        C._initVars(iso, '{Kp}= ({Pressure}-%20.16g)/ %20.16g'%(PInf,adimKp_loc))
        # Cf
        C._initVars(iso, '{Cf}= {frictionMagnitude}/ %20.16g'%adimKp_loc)

        # Passage struct
        iso[0] = '%f'%rad
        iso = C.convertBAR2Struct(iso)
        #C.convertPyTree2File(iso, 'slice%f.cgns'%rad)

        # Detection et repere local
        (PLE,PTE,PF,Glob2Loc) = detectBA_BF_QC(iso, rad)

        # CnM2
        if adimCnM2 == 0: adimCnM2_loc = 0.5*RoInf*ASOUND**2*CHORD
        else: adimCnM2_loc = adimCnM2
        CnM2x = P.integ(iso, var='Fx')[0]
        CnM2y = P.integ(iso, var='Fy')[0]
        CnM2z = P.integ(iso, var='Fz')[0]
        CnM2x = CnM2x/adimCnM2_loc
        CnM2y = CnM2y/adimCnM2_loc
        CnM2z = CnM2z/adimCnM2_loc
        Cn = (CnM2x, CnM2y, CnM2z)
        if localFrame: Cn = numpy.dot(Glob2Loc, Cn)

        if accumulatorCnM2 is not None: accumulatorCnM2[(psi,rad)] = [Cn[0],Cn[1],Cn[2]]
        CnM2All.append([Cn[0],Cn[1],Cn[2]])

        # CmM2
        if adimCmM2 == 0: adimCmM2_loc = 0.5*RoInf*ASOUND**2*CHORD
        else: adimCmM2_loc = adimCmM2
        (CmM2x,CmM2y,CmM2z) = P.integMoment(iso, center=PF, vector=['Fx','Fy','Fz'])
        CmM2x = CmM2x/adimCmM2_loc
        CmM2y = CmM2y/adimCmM2_loc
        CmM2z = CmM2z/adimCmM2_loc
        Cm = (CmM2x, CmM2y, CmM2z)
        if localFrame: Cm = numpy.dot(Glob2Loc, Cm)

        if accumulatorCmM2 is not None: accumulatorCmM2[(psi,rad)] = [Cm[0],Cm[1],Cm[2]]
        CmM2All.append([Cm[0],Cm[1],Cm[2]])

        # Rename Xabs en Xref
        nx = Internal.getNodeFromName2(iso, 'Xabs')
        nx[0] = 'Xref'
        ny = Internal.getNodeFromName2(iso, 'Yabs')
        ny[0] = 'Yref'
        nz = Internal.getNodeFromName2(iso, 'Zabs')
        nz[0] = 'Zref'
        slices[rad] = iso

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

#============================================================
# Clean stress for KIM post-processing
# IN: teff: stress tree
# IN: cpt: flowsol and grodcoordinate container number
#============================================================
def cleanStress(teff, cpt):
    tclean = Internal.copyRef(teff)

    flowSol = Internal.getNodeFromName(tclean, Internal.__FlowSolutionCenters__)
    varnames = [z[0] for z in flowSol[2][1:]]
    for v in varnames:
        if v != 'Pressure': Internal._rmNodesByName(tclean, v)

    Internal._rmNodesByName(tclean, 'CARTESIAN')
    Internal._rmNodesByName(tclean, 'FLEX')
    Internal._rmNodesByName(tclean, '.Solver#ownData')
    Internal._rmNodesByName(tclean, 'ReferenceState')
    Internal._rmNodesByName(tclean, 'FlowEquationSet')
    Internal._rmNodesByName(tclean, 'TimeMotion')
    Internal._rmNodesByName(tclean, Internal.__GridCoordinates__+"#Init")

    Internal._renameNode(tclean, Internal.__GridCoordinates__, Internal.__GridCoordinates__+"#%d"%cpt)
    Internal._renameNode(tclean, Internal.__FlowSolutionCenters__, Internal.__FlowSolutionCenters__+"#%d"%cpt)

    listOfZones = []
    for b in Internal.getBases(tclean):
        listOfZones = listOfZones+Internal.getZones(b)

    tclean = C.newPyTree(['Base', listOfZones])

    return tclean

#============================================================
# Concatenate stress files (F-order) for KIM post-processing
# IN: t1, t2: stress trees
# IN: cpt: flowsol and grodcoordinate container number
#============================================================
def concatenateStress(t1, t2, cpt):
    for z1, z2 in zip(Internal.getZones(t1), Internal.getZones(t2)):
        coordinates1 = Internal.getNodeFromName(z1, Internal.__GridCoordinates__+"#%d"%cpt)
        coordinates2 = Internal.getNodeFromName(z2, Internal.__GridCoordinates__+"#%d"%cpt)

        flowSol1 = Internal.getNodeFromName(z1, Internal.__FlowSolutionCenters__+"#%d"%cpt)
        flowSol2 = Internal.getNodeFromName(z2, Internal.__FlowSolutionCenters__+"#%d"%cpt)

        for cname in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
            c1 = Internal.getNodeFromName(coordinates1, cname)[1]
            c2 = Internal.getNodeFromName(coordinates2, cname)[1]
            if numpy.shape(c1) == numpy.shape(c2):
                c1 = numpy.stack((c1,c2), axis=2)
            else:
                c1 = numpy.concatenate((c1,c2[:,:,numpy.newaxis]), axis=2)

            Internal.getNodeFromName(coordinates1, cname)[1] = numpy.asfortranarray(c1)

        for vname in ['Pressure']:
            v1 = Internal.getNodeFromName(flowSol1, vname)[1]
            v2 = Internal.getNodeFromName(flowSol2, vname)[1]

            if numpy.shape(v1) == numpy.shape(v2):
                v1 = numpy.stack((v1,v2), axis=2)
            else:
                v1 = numpy.concatenate((v1,v2[:,:,numpy.newaxis]), axis=2)

            Internal.getNodeFromName(flowSol1, vname)[1] = numpy.asfortranarray(v1)

    return t1

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