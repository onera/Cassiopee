import Converter.PyTree as C
import Converter.Internal as Internal
import Post.PyTree as P
import Transform.PyTree as T
import Converter.Filter as Filter
import Connector.PyTree as X
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Geom.PyTree as D

IMETHOD = 'AABB'
# Extraction de la solution sur une surface dans le plan (X,R) a partir d un maillage recouvrant
# et d une courbe ou un nuage de points definis dans le plan(X,R)
# t_sol et t_pts sont dans le systeme de coordonnees cartesiennes
# t_sol et t_pts sont des fichiers
# aX + bR + c Theta +d = 0
def extractSurface(t_sol, t_pts=None, eq=(0.,0.,0.,0.), XC=(0.,0.,0.), AXIS=(1.,0.,0.), loc='centers', cellNName='cellN',variables=[], mode='robust'):
    if loc == 'centers': cellNName2 = 'centers:'+cellNName
    else: cellNName2 = cellNName
    a = eq[0]; b = eq[1]; c = eq[2]; d=eq[3]
    if t_pts is not None:
        if isinstance(t_pts, str): Pts = C.convertFile2PyTree(t_pts)
        else: Pts = t_pts

        Pts = C.convertArray2Tetra(Pts)# si nuage de pts issu du structure -> shape 2D
        Pts = Internal.getZones(Pts)[0]
        PtsBB = G.BB(Pts,method=IMETHOD)

        PtsXR = T.cart2Cyl(Pts,XC,AXIS)
        ThetaMean = C.getMeanValue(PtsXR,"CoordinateZ")
        C.convertPyTree2File(PtsXR,"Pts_XR.cgns")
        DTheta=C.getMaxValue(PtsXR,"CoordinateZ")-C.getMinValue(PtsXR,"CoordinateZ")

    if isinstance(t_sol,str):
        h = Filter.Handle(t_sol)
        t = h.loadSkeleton()
        h._loadZonesWoVars(t)
        h._loadVariables(t, var=['centers:'+cellNName])
    else:
        t = Internal.copyRef(t_sol)# pour ne pas modifier ensuite

    T._cart2Cyl(t,XC,AXIS)

    if t_pts is None:
        varc = '{'+cellNName2+'}'
        res = P.selectCells(t,"%s>0.2"%varc)
        C._initVars(res,'{eq}=%g*{CoordinateX}+%g*{CoordinateY}+%g*{CoordinateZ}+%g'%(a,b,c,d))
        res = P.isoSurfMC(res, 'eq',0.)
    else:
        if DTheta>1e-3:
            Npts = Internal.getZoneDim(PtsXR)[1]
            CoordX=[];CoordY=[];CoordZ=[]
            for i in [0, (Npts-1)//2, Npts-1]:
                CoordX.append(C.getValue(PtsXR,'CoordinateX',i))
                CoordY.append(C.getValue(PtsXR,'CoordinateY',i))
                CoordZ.append(C.getValue(PtsXR,'CoordinateZ',i))
            DX1 = CoordX[1]-CoordX[0];DX2 = CoordX[2]-CoordX[0]
            DY1 = CoordY[1]-CoordY[0];DY2 = CoordY[2]-CoordY[0]
            DZ1 = CoordZ[1]-CoordZ[0];DZ2 = CoordZ[2]-CoordZ[0]
            V1 = DY1*DZ2-DY2*DZ1
            V2 = DZ1*DX2-DZ2*DX1
            V3 = DX1*DY2-DX2*DY1
            # (X-X0)*V1+(Y-Y0)*V2+(Z-Z0)*V3=0
            a = V1;b=V2;c=V3;d=-(CoordX[0]*V1+CoordY[0]*V2+CoordZ[0]*V3)
            C._initVars(t,'{eq}=%g*{CoordinateX}+%g*{CoordinateY}+%g*{CoordinateZ}+%g'%(a,b,c,d))
            varc = '{'+cellNName2+'}'
            res = P.selectCells(t,"%s>0.2"%varc)
            res = P.isoSurfMC(res, 'eq',0.)
        else:
            varc = '{'+cellNName2+'}'
            res = P.selectCells(t,"%s>0.2"%varc)
            res = P.isoSurfMC(res, 'CoordinateZ',ThetaMean)

    Internal._rmNodesFromType(res,"FlowSolution_t")
    T._cyl2Cart(res,XC,AXIS)
    res = C.convertArray2Tetra(res)
    res = T.join(res); res = G.close(res,1.e-6)
    # res = T.splitConnexity(res)
    # resBB = G.BB(res,method=IMETHOD)
    res = XOR.conformUnstr(res,tol=0.,itermax=1)

    if isinstance(t_sol,str):
        if len(variables)==0:# on interpole tout
            t = C.convertFile2PyTree(t_sol)
        else:
            h = Filter.Handle(t_sol)
            t = h.loadSkeleton()
            h._loadZonesWoVars(t)
            h._loadVariables(t, var=variables+[cellNName2])
    else:
        T._cyl2Cart(t, XC, AXIS)

    tBB = G.BB(t,method=IMETHOD)
    resBB = G.BB(res,method=IMETHOD)
    dnrZones = []
    for zbbd in Internal.getZones(tBB):
        if G.bboxIntersection(resBB,zbbd,method=IMETHOD,isBB=True)==1:
            zd = Internal.getNodeFromName(t,zbbd[0])
            dnrZones.append(zd)
    P._extractMesh(dnrZones,res,order=2,constraint=10.,mode=mode)
    C._rmVars(res,[cellNName])
    return res

def extractIJSurface(t_sol, t_pts, XC=(0.,0.,0.), AXIS=(1.,0.,0.), loc='centers', cellNName='cellN',variables=[]):
    # distrib en i (le long de la ligne) et j
    NI = 101; NJ = 201
    dhi = G.cart((0.,0.,0.),(1./(NI-1),1,1),(NI,1,1))
    # on peut mettre un resserrement :
    # dhi = G.enforcePlusX(dhi, ...
    dhj = G.cart((0.,0.,0.),(1./(NJ-1),1,1),(NJ,1,1))

    if loc == 'centers': cellNName2 = 'centers:'+cellNName
    else: cellNName2 = cellNName

    if isinstance(t_sol,str):
        h = Filter.Handle(t_sol)
        t = h.loadSkeleton()
        h._loadZonesWoVars(t)
        h._loadVariables(t, var=['centers:'+cellNName])
    else:
        t = Internal.copyRef(t_sol)# pour ne pas modifier ensuite

    T._cart2Cyl(t,XC,AXIS)


    if isinstance(t_pts, str): Pts = C.convertFile2PyTree(t_pts)
    else: Pts = t_pts

    Pts = C.convertArray2Tetra(Pts)# si nuage de pts issu du structure -> shape 2D
    T._cart2Cyl(Pts,XC,AXIS)
    Pts = C.convertArray2Node(Pts)
    Pts = Internal.getZones(Pts)[0]
    NPts = Internal.getZoneDim(Pts)[1]
    Pts0 = []
    for i in range(NPts):
        xi = C.getValue(Pts,'CoordinateX',i)
        yi = C.getValue(Pts,'CoordinateY',i)
        zi = C.getValue(Pts,'CoordinateZ',i)
        Pts0.append((xi,yi,zi))
    # creation de la spline cubique
    res = D.polyline(Pts0)
    #res = G.map(res,dhi,dir=1)

    ThetaMean = C.getMeanValue(res,'CoordinateZ')
    #print('ThetaMean = ', ThetaMean, C.getMeanValue(Pts,'CoordinateZ'), C.getMaxValue(Pts,'CoordinateZ')-C.getMinValue(Pts,'CoordinateZ'))
    # PLAN MERIDIEN Theta=ThetaMean
    PM = P.isoSurfMC(t,"CoordinateZ",ThetaMean)
    PM = T.join(PM); PM = T.splitConnexity(PM)
    # on affine
    PM = [z0 for z0 in PM if G.bboxIntersection(z0,res)==1]

   # on essaie de trouver les angles optimaux - pour les cas centrifuges
    RMin = C.getMinValue(res,'CoordinateY')# %R
    RMax = C.getMaxValue(res,'CoordinateY')
    XMin = C.getMinValue(res,'CoordinateX')
    XMax = C.getMaxValue(res,'CoordinateX')
    C.convertPyTree2File(res,"toto1.cgns")
    C.convertPyTree2File(PM,"toto2.cgns")

    T._projectOrtho(res,PM)
    del PM

    intersectedZones = []
    for z in Internal.getZones(t):
        if G.bboxIntersection(res,z)==1: intersectedZones.append(z)

    C._initVars(intersectedZones,'{tag}=logical_and({CoordinateX}>%g,{CoordinateX}<%g)'%(XMin,XMax))
    C._initVars(intersectedZones,'{tag}={tag}*logical_and({CoordinateY}>%g,{CoordinateY}<%g)'%(RMin,RMax))
    intersectedZones = P.selectCells2(intersectedZones,'tag')

    ThetaMax = C.getMaxValue(intersectedZones,'CoordinateZ')
    ThetaMin = C.getMinValue(intersectedZones,'CoordinateZ')
    del intersectedZones
    print(' RMin = %g %g , XMin = %g %g '%(RMin,RMax,XMin,XMax))
    print(' ThetaMin/Mean/Max : ', ThetaMin, ThetaMean, ThetaMax)

    # EXTRUSION DE LA LIGNE
    res1 = T.translate(res,(0,0,-ThetaMean))
    NP = int(NJ/2)+1
    T._addkplane(res1,N=NP)
    DZP = ThetaMax-ThetaMean; DZ = DZP/NP
    C._initVars(res1,"{CoordinateZ}={CoordinateZ}*%g+%g"%(DZ,ThetaMean))
    #
    res2 = T.translate(res,(0,0,-ThetaMean))
    NP2 = NJ-NP+1
    T._addkplane(res2,N=NP2)
    DZP = ThetaMean-ThetaMin; DZ = DZP/NP
    C._initVars(res2,"{CoordinateZ}={CoordinateZ}*%g+%g"%(DZ,ThetaMin)); res2[0]='toto'
    T._reorder(res2,(-1,2,3))
    res = T.join(res1,res2)
    del res1, res2
    res = G.map(res,dhj,dir=1)
    Internal._rmNodesFromType(res,'FlowSolution_t')
    T._cyl2Cart(res,XC,AXIS)

    #C._cpVars(t,"centers:cellnf",t,"centers:cellN")
    C._initVars(t,'{centers:cellN}={%s}'%cellNName2)
    C.convertPyTree2File(res,"meshij.cgns")
    P._extractMesh(t,res, mode='robust',constraint=0.)
    #C._initVars(res,"{cellN}=({cellN}>0.5)")
    C.convertPyTree2File(res,"out.cgns")
    return res
