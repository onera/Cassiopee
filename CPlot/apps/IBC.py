#============================================================================================
# Creation des premiers pts IBC
# IN : a : zone avec cellN=2 pour les pts IBC
# OUT : retourne les pts IBC sous forme de 'NODE'
#============================================================================================
def getIBCFrontForZone__(a):
    f0 =  P.selectCells(a, '{centers:cellN} == 2.')
    f0 = T.join(f0); f0 = G.close(f0)
    #
    # recuperation des champs en centres perdus par selectCells
    a2 = C.initVars(a,'centers:cellN',1.)
    ta = C.newPyTree(['Base']); ta[2][1][2] = [a2]
    tf0 = C.newPyTree(['Base']); tf0[2][1][2] = [f0]
    tf0 = P.extractMesh(ta,tf0)
    tf0 = C.rmVars(tf0,['centers:cellN','cellN'])
    coords = C.getFields(Internal.__GridCoordinates__,tf0)
    solc =  C.getFields(Internal.__FlowSolutionCenters__,tf0)
    coords = Converter.node2Center(coords) # passage en centres pour recuperer les coordonnees des centres
    coords = Converter.addVars([coords,solc])
    # passage en 'NODE'
    coords = Converter.convertArray2Node(coords)[0]
    # Sortie : conversion en zone node CGNS
    return C.convertArrays2ZoneNode('Front1_'+a[0], [coords])

#==================================================================================================================
# Calcul de la normale et de delta pour les pts du front
# IN : fc1 : front IBC sous forme de zone
# IN : donorZone : zone dont provient la zone IBC, localisee comme souhaitee pour les interpolations (ex en centres avec cell fict)
# IN : numZone numero de la zone de donorZones contenant le front, demarre a 1
# IN : listPts : si != [] alors delta = delta+eps pour ces points
# IN : eps : deplacement des points d indices dans listPts
# OUT : front fc1 avec delta (eventuellement corrige de eps)
#                      la normale = delta * grad TurbulentDistance
#                      le numero du bloc d origine 'noblk', demarrant a 1
#                      l indice en centres avec ghost cells dans le maillage d origine correspondant au pt IBC 
#==================================================================================================================
def getIBCFrontInfo(fc1, donorZone, numZone, dhloc, listPts =[], eps=0.):
    toldist = 1.e-10
    # Ajout du numero du bloc d origine
    fc1 = C.initVars(fc1,'noblk',numZone)

    # Determination de l indice du pt dans le maillage donneur
    hook = C.createHook(donorZone, function='nodes')
    indices,distances = C.nearestNodes(hook,fc1)
    C.freeHook(hook)
    coords1 = C.getFields(Internal.__GridCoordinates__,fc1)[0]
    coords1 = Converter.initVars(coords1,'ind',-1.)
    coords1 = Converter.initVars(coords1,'indI',-1.)
    coords1 = Converter.initVars(coords1,'indJ',-1.)
    coords1 = Converter.initVars(coords1,'indK',-1.)
    coords1 = Converter.extractVars(coords1,['ind','indI','indJ','indK'])

    npts = coords1[1].shape[1]
    dimGC = Internal.getZoneDim(donorZone)
    nigc = dimGC[1]; njgc = dimGC[2]; nkgc = dimGC[3]; nigcnjgc = nigc*njgc
    dhLoc = [dhloc]*npts # tableau dhLoc
    xt = C.getField('CoordinateX',donorZone)[0][1]
    yt = C.getField('CoordinateY',donorZone)[0][1]
    zt = C.getField('CoordinateZ',donorZone)[0][1]
    
    for ind in xrange(npts):
        index = indices[ind]-1
        dist  = distances[ind]
        if dist < toldist:
            indk = index/nigcnjgc
            indj = (index-indk*nigcnjgc)/nigc
            indi = index - indj*nigc - indk*nigcnjgc
            coords1[1][0,ind] = index
            coords1[1][1,ind] = indi
            coords1[1][2,ind] = indj
            coords1[1][3,ind] = indk
            if indi == nigc-1: dxloc = abs(xt[0,index]-xt[0,index-1])
            else: dxloc = abs(xt[0,index]-xt[0,index+1])
            if indj == njgc-1: dyloc = abs(yt[0,index]-yt[0,index-nigc])
            else: dyloc = abs(yt[0,index]-yt[0,index+nigc])
            if indk == nkgc-1: dzloc = abs(zt[0,index]-zt[0,index-nigcnjgc])
            else: dzloc = abs(zt[0,index]-zt[0,index+nigcnjgc])
            
            if dxloc < toldist: dxloc = 1e10
            if dyloc < toldist: dyloc = 1e10
            if dzloc < toldist: dzloc = 1e10            
            dhLoc[ind] = min(dxloc,dyloc,dzloc)

    C.setFields([coords1], fc1, loc='nodes')

    # delta * normale
    varnx = 'gradxTurbulentDistance'
    varny = 'gradyTurbulentDistance'
    varnz = 'gradzTurbulentDistance'
    fc1 = C.normalize(fc1, [varnx, varny, varnz])
    #
    if listPts == []:
        fc1 = C.initVars(fc1,'delta',0.)
        deltaa = C.getField('delta',fc1)[0]
        distance = C.getField('TurbulentDistance',fc1)[0][1]
        # formule d obtention du frontC2
        for ind in xrange(deltaa[1].shape[1]):
            dist = distance[0,ind]
            # NOUVELLE VERSION
            # # cas 1 : le centre est proche paroi, le point interpole est alors positionne a dhloc+eps de la paroi
	    # if abs(dist) < dhLoc[ind]: deltaa[1][0,ind] =  2*dhLoc[ind] + eps
            # # cas 2 : le centre est loin de la paroi, le point interpole est alors positionne a dist+eps de la paroi
	    # else: deltaa[1][0,ind] = 2.*abs(dist) + eps
            # FIN NOUVELLE VERSION

            # cas 1 : le centre est proche paroi, le point interpole est alors positionne a dhloc+eps de la paroi
	    if abs(dist) < dhloc: deltaa[1][0,ind] = abs(dist) + dhloc + eps
            # cas 2 : le centre est loin de la paroi, le point interpole est alors positionne a dist+eps de la paroi
	    else: deltaa[1][0,ind] = 2.*abs(dist) + eps
        C.setFields([deltaa], fc1, loc='nodes')

    else:
        # modification locale de delta : a modifier ?
        deltaa = C.getField('delta',fc1)[0]
        for ind in listPts: deltaa[1][0,ind] += eps
        C.setFields([deltaa], fc1, loc='nodes')
    #
    fc1 = C.initVars(fc1, '{nx} = {gradxTurbulentDistance} * {delta}')
    fc1 = C.initVars(fc1, '{ny} = {gradyTurbulentDistance} * {delta}')
    fc1 = C.initVars(fc1, '{nz} = {gradzTurbulentDistance} * {delta}')
    return fc1
