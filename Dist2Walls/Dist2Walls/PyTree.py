#
# Python Interface to compute distance to walls from PyTrees
#
from . import Dist2Walls
__version__ = Dist2Walls.__version__

try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except:
    raise ImportError("Dist2Walls: requires Converter module.")

PHIMAX = 1.e12
# Les differents algorithmes qu'on peut choisir pour la resolution de l'equation Eikonale
fmm = 0; fim = 1
# temporaire
fim_old = 2

#==============================================================================
# Calcul de la distance a la paroi pour a (tree, base, zone)
#==============================================================================
def distance2Walls(t, bodies, type='ortho', loc='centers', signed=0, dim=3):
    """Compute distance field.
    Usage: distance2Walls(a, bodies, type, loc, signed, dim)"""
    tp = Internal.copyRef(t)
    _distance2Walls(tp, bodies, type, loc, signed, dim)
    return tp

#==============================================================================
def _distance2Walls(t, bodies, type='ortho', loc='centers', signed=0, dim=3):
    """Compute distance field.
    Usage: distance2Walls(a, bodies, type, loc, signed, dim)"""
    bodyZones = Internal.getZones(bodies)
    bodiesa = C.getFields(Internal.__GridCoordinates__, bodyZones)
    cellnba = [] # cellN localise au meme endroit que bodies
    varn = 'cellN'; varn2 = 'cellnf'
    for zb in bodyZones:
        posn = C.isNamePresent(zb, varn)
        if posn != -1: cellnba += C.getField(varn, zb)
        else:
            posn = C.isNamePresent(zb, varn2)
            if posn != -1: cellnba += C.getField(varn2, zb)
            else: cellnba.append([])
    
    coords = C.getFields(Internal.__GridCoordinates__, t)
    if loc == 'centers': flag = C.getField('centers:flag',t)
    else: flag = C.getField('flag',t)

    distances = Dist2Walls.distance2Walls(
        coords, bodiesa, flags=flag, cellnbodies=cellnba, type=type,
        loc=loc, signed=signed, dim=dim)
    if loc == 'centers': return C.setFields(distances, t, 'centers')
    else: return C.setFields(distances, t, 'nodes')

#==============================================================================
# Eikonal equation starting from spring points 
# Multidomain not taken into account, no transfer is done !  
#==============================================================================
def eikonal(t,tc=None,loc='nodes',nitmax=10, err=0.01,algo=fim_old):
    """Solve the eikonal equation.
    Usage: eikonal(t,loc='nodes')"""
    tp = Internal.copyRef(t)
    _eikonal(tp,tc=tc,loc=loc, nitmax=nitmax, err=err, algo=algo)
    return tp

def _eikonal(t,tc=None,loc='nodes', nitmax=10, err=0.01,algo=fmm):
    MB = 1
    try: import Connector.PyTree as X
    except: MB = 0

    if tc is None or MB==0:
        for z in Internal.getNodesFromType2(t, 'Zone_t'):
            _eikonalForZone(z,loc,algo)
        return None

    # Eikonal
    nzones = len(Internal.getNodesFromType2(t, 'Zone_t'))
    isConverged=[0]*nzones
    it = 0
    nocv = 0 # nb de zones convergees
    # calcul de la taille de maille sur le niveau le plus fin
    dhmin = PHIMAX
    for z in Internal.getNodesFromType2(t, 'Zone_t'):
        dhmin = min(dhmin,C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
    while nocv < nzones and it < nitmax+1:
        print('Iteration %d'%it)
        # Eikonal sur les zones non convergees et sources
        if loc == 'nodes': C._initVars(t,'{PhiM}={Phi}')
        else: C._initVars(t,'{centers:PhiM}={centers:Phi}')
        no = 0 
        for z in Internal.getNodesFromType2(t,"Zone_t"): 
            if isConverged[no] < 1 and C.getMaxValue(z,loc+':flag')>0.:
                if loc == 'nodes': C._initVars(z,'{Phi}=({flag}<0.1)*%f+({flag}>0.)*{Phi}'%PHIMAX)
                else: C._initVars(z,'{centers:Phi}=({centers:flag}<0.1)*%f+({centers:flag}>0.)*{Phi}'%PHIMAX)
                _eikonalForZone(z,loc=loc,algo=algo)
                isConverged[no] = -1
            no+=1
        # Synchro
        print('Synchronization/transfers')
        if tc is not None:
            no = 0
            for z in Internal.getNodesFromType2(t,"Zone_t"):
                if isConverged[no] == -1: # a ete eikonalise: transferts
                    z2 = Internal.getNodeFromName(tc,z[0])
                    C._cpVars(z,loc+':Phi',z2,'Phi')
                    C._initVars(z2,'flag',1.)
                    # PAS THREADE ?????
                    X._setInterpTransfers(t,z2, variables=['Phi','flag'],variablesIBC=None)
                no += 1

        # Convergence 
        if it > 0:
            no = 0 ; nocv = 0
            for z in Internal.getNodesFromType2(t,"Zone_t"):
                if loc=='nodes':C._initVars(z,'{DPhi}=abs({Phi}-{PhiM})/(maximum(1.e-12,{Phi}))')
                else:C._initVars(z,'{centers:DPhi}=abs({centers:Phi}-{centers:PhiM})/(maximum(1.e-12,{centers:Phi}))')
                valmax = C.getMaxValue(z,loc+':DPhi')
                dhloc = C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0)
                errloc = err*dhloc/dhmin # on augmente la tolerance au fur et a mesure qu on diminue la resolution
                if valmax< errloc and isConverged[no]==-1:
                    isConverged[no] = 1
                    nocv += 1
                elif isConverged[no] == 1: nocv+=1
                no += 1

        # Iteration 
        it += 1
    #-----------------------------------------------------------------------------
    if it < nitmax+1: print('Distance by Eikonal converged after %d subiterations.'%it)
    else: 
        print('Warning: distance by Eikonal did not converged after %d subiterations.'%nitmax)
        noi = 0
        for i in isConverged: 
            if i != 1: print('%d, %d'%(i, noi))
            noi += 1
    return None

#==============================================================================
# Eikonal equation starting from spring points  for a zone
# IN: Phi > 1.e10 for non-BC points, initial value elsewhere
# In place version
#==============================================================================
def _eikonalForZone(z,loc='nodes',algo=fim_old):
    """Solve the eikonal equation.
    Usage: _eikonalForZone(z,loc='nodes')"""
    locNode = False
    if loc == 'nodes': locNode = True
    
    nodes = C.getFields(Internal.__GridCoordinates__, z)[0]
    if locNode: 
        phi = C.getField('Phi',z)[0]
        speed = C.getField('speed',z)[0]
    else:
        phi = C.getField('centers:Phi',z)[0]
        speed = C.getField("centers:speed",z)[0]
        nodes = Converter.node2Center(nodes)
    if phi == []: 
        raise ValueError("Dist2Walls:_eikonalForZone: Phi variable not defined in zone.")
    if speed == []:
        raise ValueError("Dist2Walls:_eikonalForZone: speed variable not defined in zone.")

    nodes = Converter.addVars([nodes,phi,speed])
    nodes = Dist2Walls.eikonal(nodes,algo)
    phi = Converter.extractVars([nodes],['Phi'])
    C.setFields(phi,z,loc)
    return None

#------------------------------------------------------------------
# transfert du cellN aux raccords
# min/max cellN
#------------------------------------------------------------------
def transferCellN__(t,tc,DEPTH,loc):
    if tc is None: return t
    try:
        import Connector.PyTree as X
    except:
        raise ImportError("Dist2Walls: Eikonal version requires Connector module.")
    # POINTS EXTERIEURS
    # Marquage des pts de front entre du 0 et du 1
    t = X.setHoleInterpolatedPoints(t,depth=DEPTH,loc=loc)
    # transfert du cellN aux raccords
    if tc is not None:
        C._cpVars(t,loc+':cellN',tc,'cellN')
        for zc in Internal.getNodesFromType2(tc,"Zone_t"):
            if C.getMaxValue(zc,'cellN')==2.:
                X._setInterpTransfers(t,zc,variables=["cellN"])
    if loc == 'nodes':
        C._initVars(t,"{cellN}=({cellN}>1.5)*2.+({cellN}>0.)*({cellN}<1.5)")
        C._initVars(t,'{cellN} = 1-{cellN}+({cellN}>1.5)*3')
    else:
        C._initVars(t,"{centers:cellN}=({centers:cellN}>1.5)*2.+({centers:cellN}>0.)*({centers:cellN}<1.5)")
        C._initVars(t,'{centers:cellN} = 1-{centers:cellN}+({centers:cellN}>1.5)*3')
    t = X.setHoleInterpolatedPoints(t,depth=DEPTH,loc=loc)
    if loc == 'nodes':
        C._initVars(t,'{cellN} = 1-{cellN}+({cellN}>1.5)*3')
        C._initVars(t,'{cellN2}={cellN}')
    else:
        C._initVars(t,'{centers:cellN} = 1-{centers:cellN}+({centers:cellN}>1.5)*3')
        C._initVars(t,'{centers:cellN2}={centers:cellN}')

    # POINTS INTERIEURS
    if loc == 'nodes': C._initVars(t,'{cellN}=minimum(1.,{cellN})')
    else: C._initVars(t,'{centers:cellN}=minimum(1.,{centers:cellN})')
        
    t = X.setHoleInterpolatedPoints(t,depth=-DEPTH,loc=loc)
    # transfert du cellN aux raccords
    if tc is not None:
        C._cpVars(t,loc+':cellN',tc,'cellN')
        for zc in Internal.getNodesFromType2(tc,"Zone_t"):
            if C.getMaxValue(zc,'cellN')==2.:
                X._setInterpTransfers(t,zc,variables=["cellN"])
    if loc == 'nodes':
        C._initVars(t,"{cellN}=({cellN}>1.5)*2.+({cellN}>0.)*({cellN}<1.5)")
        C._initVars(t,'{cellN} = 1-{cellN}+({cellN}>1.5)*3')
    else:
        C._initVars(t,"{centers:cellN}=({centers:cellN}>1.5)*2.+({centers:cellN}>0.)*({centers:cellN}<1.5)")
        C._initVars(t,'{centers:cellN} = 1-{centers:cellN}+({centers:cellN}>1.5)*3')
    t = X.setHoleInterpolatedPoints(t,depth=-DEPTH,loc=loc)
    if loc == 'nodes':
        C._initVars(t,'{cellN} = 1-{cellN}+({cellN}>1.5)*3')
        C._initVars(t,'{cellN}=maximum({cellN},{cellN2})')
        C._rmVars(t,['cellN2'])
    else:
        C._initVars(t,'{centers:cellN} = 1-{centers:cellN}+({centers:cellN}>1.5)*3')
        C._initVars(t,'{centers:cellN}=maximum({centers:cellN},{centers:cellN2})')
        C._rmVars(t,['centers:cellN2'])
    return t
    
#=============================================================================
# Distance field computation using FIM method
# IN: err: relative error on distance at convergence
# IN: nitmax: nb of iterative loops for multidomain
#=============================================================================
def distance2WallsEikonal(t, body, tc=None, DEPTH=2, loc='nodes', err=0.01, nitmax=10, type=0,algo=fim_old):
    #import time
    #beg = time.time()
    flagName='flag'; distName='TurbulentDistance'
    if loc == 'centers':
        flagName='centers:'+flagName; distName='centers:'+distName

    if loc == 'nodes':
        C._initVars(t,'{sign}=({cellN}>0.)-1.*({cellN}<1.)')
        # Marquage des points 
        C._initVars(t,'{Phi}=%g*{flag}'%PHIMAX)#'Phi=%g*({flag}<1.)+({flag}>0.)'%PHIMAX)

    else:
        C._initVars(t,'{centers:sign}=({centers:cellN}>0.)-1.*({centers:cellN}<1.)')
        # Marquage des points 
        C._initVars(t,'{centers:Phi}=%g*({centers:flag}<1.)+({centers:flag}>0.)'%PHIMAX)

    #----------------------------------------------
    # Marquage des pts de front entre du 0 et du 1
    #----------------------------------------------
    #print('transfer cellN : a passer par la fonction recente de Connector')
    t = transferCellN__(t,tc,DEPTH,loc)

    # Initialisation du front
    #print('initDistance')
    #beg2 = time.time()
    for z in Internal.getZones(t):
        dims = Internal.getZoneDim(z)
        if dims[0] != 'Structured':
            raise ValueError('dist2WallsEikonal works only on structured grids currently.')
        #beg4 = time.time()
        C._initVars(t,distName,PHIMAX)
        #end4 = time.time()
        #print("Temps init vars phi : {} secondes".format(end4-beg4))
        # calcul de la distance a la paroi reelle
        if C.getMaxValue(z,flagName) == 1.:
            #beg3 = time.time()
            _distance2Walls(z,body,type='ortho',loc=loc)
            #end3 = time.time()
            #print("Calcul distance initiale : {} secondes".format(end3-beg3))

        #beg5 = time.time()
        if loc == 'nodes': 
            C._initVars(z,'{Phi}={TurbulentDistance}*({flag}>0.)+%g*({flag}<1.)'%PHIMAX)
        else:
            C._initVars(z,'{centers:Phi}={centers:TurbulentDistance}*({centers:flag}>0.)+%g*({centers:flag}<1.)'%PHIMAX)
        #end5 = time.time()
        #print("Temps passe init var turbulentDistance : {} secondes".format(end5-beg5))
            
        if type == 0: 
            ni = dims[1]; nj = dims[2]; nk = dims[3]
            i = max(1,ni//2); j = max(1,nj//2); k = max(1,nk//2)
            ind1 = i+j*ni+k*ni*nj
            ind2 = ind1+1
            dh = C.getValue(z,'CoordinateX',ind2)-C.getValue(z,'CoordinateX',ind1)
            C._initVars(z,loc+':speed',1./dh)
             
        else: C._initVars(z,loc+':speed',1.)
    #end2 =time.time()
    #print("Temps initialisation au front : {} secondes".format(end2-beg2))
    #end = time.time()
    #print("Temps initialisation champs pour l'Eikonal : {} secondes".format(end-beg))
    # Eikonal
    #print('eikonal')
    _eikonal(t,tc,loc=loc, nitmax=nitmax, err=err,algo=algo)

    #-----------------------------------------------------------------------------
    if loc =='nodes':
        C._initVars(t,'{TurbulentDistance}={sign}*{Phi}')
        C._rmVars(t,['flag','PhiM','DPhi','speed','Phi']) # pour l instant on detruit tout
    else:
        C._initVars(t,'{centers:TurbulentDistance}={centers:sign}*{centers:Phi}')
        C._rmVars(t,['centers:flag','centers:PhiM','centers:DPhi','centers:speed','centers:Phi']) # pour l instant on detruit tout
    return t
