#Class for coupled FastLBM - FastS simulations
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast
import FastLBM.PyTree as FastLBM
import FastS.PyTree as FastS
import Connector.PyTree as X
import Transform.PyTree as T
from Apps.Fast.Common import Common
import numpy

try: range = xrange
except: pass

#===============================================================================
#   Multiblock prepare
#   NP : nombre de processeurs
#===============================================================================

def prepare(t_case, t_out, tc_out, solvertype=None, translation=[0.,0.,0.], NP=0, format='signle', NG=2):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    #sequential prep
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, solvertype, translation, NP, format, NG=NG)
    return ret

def prepare0(t_case, t_out, tc_out, solvertype=None, translation=[0.,0.,0.], NP=0, format='single',NG=2):
    #Cas ou l'on specifie un fichier .cgns pour le t_case
    if isinstance(t_case, str): t = C.convertFile2PyTree(t_case)
    else: t = t_case

    TX = translation[0]; TY = translation[1]; TZ = translation[2]

    #Recupere la dimension. Si non specifiee : dim = 3
    dim = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)
    else: C._addState(t, 'EquationDimension', 3)

    #Regarde si on a defini des equations
    eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
    model = 'CouplageNSLBM'
    if eqs is not None: model = Internal.getValue(eqs)
    else: C._addState(t, 'GoverningEquations','CouplageNSLBM')

    #On affecte a chaque zone son modele : LBM, Euler, NS,...
    #Attention !! Dans solvertype on doit mettre les zones dans le meme ordre que dans l'arbre !!
    if solvertype==None: raise ValueError('prepare0 : solvertype has not been specified.')
    infos_couplage_zones = {}
    zones_LBM = []; zones_NS = []
    for i,z in enumerate(Internal.getZones(t)):
        #On initialise le dico qui contiendra les infos pour les BC de couplage
        infos_couplage_zones[z[0]] = []
        #Affectation des equations a chaque zone
        if solvertype[i]==0:
            C._addState(z,'GoverningEquations','LBMLaminar')
            zones_LBM.append(z[0])
        elif solvertype[i]==1:
            C._addState(z,'GoverningEquations','Euler')
            zones_NS.append(z[0])
        elif solvertype[i]==2:
            C._addState(z,'GoverningEquations','NSLaminar')
            zones_NS.append(z[0])
        elif solvertype[i]==3:
            C._addState(z,'GoverningEquations','NSTurbulent')
            zones_NS.append(z[0])
        else:
            raise ValueError('prepare0 : solvertype only takes values 0,1,2,3. ')

    # 1. PREMIERE RECHERCHE DES INTERFACES NS_LBM
    #    La fonction connectNSLBM cherche les interfaces 1-to-1, c'est-a-dire
    #    les interfaces qui correspondent a un bord complet entre domaines.
    #===========================================================================
    nb_couplage_bc = 0
    t = X.connectNSLBM(t)
    nb_vert = []
    for z in Internal.getZones(t):
        dim = Internal.getZoneDim(z)
        nb_vert.append(dim[1:4])
        zoneBC = Internal.getNodesFromType2(z, 'BC_t')
        if zoneBC !=[]:
            for bc in zoneBC:
                v = Internal.getValue(bc)
                if v=='BCReconsLBM' or v=='BCdimNS':
                    infos_couplage_zones[z[0]].append(Internal.copyNode(bc))
                    nb_couplage_bc += 1

    print('Found {} internal NSLBM boundary conditions corresponding to {} NSLBM interfaces.'.format(nb_couplage_bc,int(nb_couplage_bc/2)))
    #Internal._rmNodesByType(t,'ZoneBC_t')
    C._rmBCOfType(t,'BCdimNS')
    C._rmBCOfType(t,'BCReconsLBM')
    Internal._rmNodesByType(t,'ZoneGridConnectivity_t')

    isPerio = 0
    if NP>1: t = T.splitSize(t, R=NP, type=2, minPtsPerDir=9)

    # 2. MISE EN PLACE DES CONNECTIONS POUR LES TRANSFERTS
    #===========================================================================
    present_match = []

    # 2.1 On commence par connecter les bords qui se touchent en 1-to-1
    #===========================================================================
    t = X.connectMatch(t)
    for z in Internal.getZones(t):
        match = Internal.getNodesFromType(z,'GridConnectivity1to1_t')
        for m in match:
            present_match.append(m[0])

    # 2.2 Gestion de la periodicite du domaine complet selon X et Y
    #     Si on a une pêriodicie entre un domaine NS et LBM alors on stocke les
    #     informations pour creer la BC plus tard.
    #===========================================================================
    count=0
    if TX != 0.:
        isPerio=1
        t = X.connectMatchPeriodic(t,translation=[TX,0.,0.])
        for z in Internal.getZones(t):
            match = Internal.getNodesFromType(z,'GridConnectivity1to1_t')
            for m in match:
                #On ne s'interesse qu'aux nouvelles connections qui sont forcement
                #issues de la periodicite selon X.
                if m[0] not in present_match:
                    present_match.append(m[0])
                    name_other = Internal.getValue(m)
                    if (z[0] in zones_NS) and (name_other in zones_LBM):
                        #Dans ce cas on a une CL periodique avec une frontiere NSLBM
                        node = Internal.getNodeFromName(m, 'PointRange')
                        per_range = Internal.getValue(node)
                        per_range = sum(per_range.tolist(),[])
                        C._addBC2Zone(z,'match_PerX'+z[0],'BCdimNS',per_range)
                        node = Internal.getNodesFromValue(z,'BCdimNS')
                        if node != []:
                            for n in node:
                                infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                        count += 1
                    elif (z[0] in zones_LBM) and (name_other in zones_NS):
                        #Dans ce cas on a une CL periodique avec une frontiere NSLBM
                        node = Internal.getNodeFromName(m, 'PointRange')
                        per_range = Internal.getValue(node)
                        per_range = sum(per_range.tolist(),[])
                        C._addBC2Zone(z,'match_PerX'+z[0],'BCReconsLBM',per_range)
                        node = Internal.getNodesFromValue(z,'BCReconsLBM')
                        if node != []:
                            for n in node:
                                infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                        count += 1
        print('Found {} X-periodic NSLBM boundary conditions corresponding to {} NSLBM interfaces.'.format(count,int(count/2)))
        #Internal._rmNodesByType(t,'ZoneBC_t')
        C._rmBCOfType(t,'BCdimNS')
        C._rmBCOfType(t,'BCReconsLBM')

    count = 0
    if TY != 0.:
        isPerio=1
        t = X.connectMatchPeriodic(t,translation=[0.,TY,0.])
        for z in Internal.getZones(t):
            match = Internal.getNodesFromType(z,'GridConnectivity1to1_t')
            for m in match:
                #On ne s'interesse qu'aux nouvelles connections qui sont forcement
                #issues de la periodicite selon Y.
                if m[0] not in present_match:
                    present_match.append(m[0])
                    name_other = Internal.getValue(m)
                    if (z[0] in zones_NS) and (name_other in zones_LBM):
                        #Dans ce cas on a une CL periodique avec une frontiere NSLBM
                        node = Internal.getNodeFromName(m, 'PointRange')
                        per_range = Internal.getValue(node)
                        per_range = sum(per_range.tolist(),[])
                        C._addBC2Zone(z,'match_PerY'+z[0],'BCdimNS',per_range)
                        node = Internal.getNodesFromValue(z,'BCdimNS')
                        if node != []:
                            for n in node:
                                infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                        count += 1
                    elif (z[0] in zones_LBM) and (name_other in zones_NS):
                        #Dans ce cas on a une CL periodique avec une frontiere NSLBM
                        node = Internal.getNodeFromName(m, 'PointRange')
                        per_range = Internal.getValue(node)
                        per_range = sum(per_range.tolist(),[])
                        C._addBC2Zone(z,'match_PerY'+z[0],'BCReconsLBM',per_range)
                        node = Internal.getNodesFromValue(z,'BCReconsLBM')
                        if node != []:
                            for n in node:
                                infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                        count += 1
        print('Found {} Y-periodic NSLBM boundary conditions corresponding to {} NSLBM interfaces.'.format(count,int(count/2)))
        #Internal._rmNodesByType(t,'ZoneBC_t')
        C._rmBCOfType(t,'BCdimNS')
        C._rmBCOfType(t,'BCReconsLBM')

    # 2.3 Creation des zones periodiques en X et Y
    #     Pour assurer les bons transferts, on va dupliquer les zones par
    #     translation en X et Y.
    #===========================================================================
    # Avant la duplication, on garde une trace des zones d'origine
    zonesNamesORIG=[]
    for z in Internal.getZones(t): zonesNamesORIG.append(z[0])

    # Creation des zones periodiques : periodicite par translation en x et y
    dictOfTempPerNodes={}
    if isPerio:
        C._addPeriodicZones__(t[2][1])

        # Rajout des coins car pour la LBM il faut avoir les bonnes VDF dans
        # les coins sinon valeur aberrantes.
        zonesDup=[]; zonesDupNames=[]
        for z in Internal.getZones(t[2][1]):
            isperiod = Internal.getNodeFromName1(z, 'TempPeriodicZone')
            if isperiod is not None:
                C._rmBCOfType(z,'BCMatch')
                zonesDup.append(z)
                zonesDupNames.append(z[0])
        zcorners=[]
        for z in Internal.getZones(t):
            isperiod = Internal.getNodeFromName1(z, 'TempPeriodicZone')
            if isperiod is None:
                for j in [-1,1]:
                    for i in [-1,1]:
                        zdup = T.translate(z,(i*TX,j*TY,0.))
                        C._rmBCOfType(zdup,"BCMatch")
                        zdupname = C.getZoneName(Internal.getName(zdup)+'_dup')
                        zdup[0] = zdupname
                        zonesDup.append(zdup)
                        #
                        zonesDup = X.connectMatch(zonesDup)
                        zdup=zonesDup[-1]
                        gc = Internal.getNodeFromType(zdup,'GridConnectivity1to1_t')
                        if gc is not None:
                            zdonorname = Internal.getValue(gc)
                            Internal.createChild(zdup,'TempPeriodicZone','UserDefinedData_t',value=zdonorname,children=None)
                            zcorners.append(zdup)
                        del zonesDup[-1]

        del zonesDup
        t[2][1][2]+= zcorners
        del zcorners

        dictOfTempPerNodes={}
        for z in Internal.getZones(t):
            pernode = Internal.getNodeFromName1(z,'TempPeriodicZone')
            if pernode is not None:
                dictOfTempPerNodes[z[0]]=pernode
        Internal._rmNodesFromName2(t[2][1],'TempPeriodicZone')

    # duplication par periodicite en Z
    if isPerio: C._rmBCOfType(t,"BCMatch")
    #en fait, il faudrait enlever la notion de periodicite, noeud GridConnectivityProperty
    # pour ne pas redupliquer des zones

    # 2.2 Gestion de la periodicite du domaine complet selon Z
    #     Si on a une pêriodicie entre un domaine NS et LBM alors on stocke les
    #     informations pour creer la BC plus tard.
    #===========================================================================
    count = 0
    if TZ !=0.:
        isPerio = 1
        t = X.connectMatchPeriodic(t,translation=[0.,0.,TZ])
        for z in Internal.getZones(t):
            if z[0] in zonesNamesORIG:
                match = Internal.getNodesFromType(z,'GridConnectivity1to1_t')
                for m in match:
                    val = Internal.getNodeFromType(m,'Periodic_t')
                    if val is not None:
                        name_other = Internal.getValue(m)
                        if (z[0] in zones_NS) and (name_other in zones_LBM):
                            node = Internal.getNodeFromName(m,'PointRange')
                            per_range = Internal.getValue(node)
                            per_range = sum(per_range.tolist(),[])
                            C._addBC2Zone(z,'match_PerZ'+z[0],'BCdimNS',per_range)
                            node = Internal.getNodesFromValue(z,'BCdimNS')
                            if node !=[]:
                                for n in node:
                                    infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                            count += 1
                        elif (z[0] in zones_LBM) and (name_other in zones_NS):
                            node = Internal.getNodeFromName(m,'PointRange')
                            per_range = Internal.getValue(node)
                            per_range = sum(per_range.tolist(),[])
                            C._addBC2Zone(z,'math_PerZ'+z[0],'BCReconsLBM',per_range)
                            node = Internal.getNodesFromValue(z,'BCReconsLBM')
                            if node != []:
                                for n in node:
                                    infos_couplage_zones[z[0]].append(Internal.copyNode(n))
                            count += 1
        print('Found {} Z-periodic NSLBM boundary conditions correspondig to {} NSLBM interfaces.'.format(count,int(count/2.)))
        #Internal._rmNodesByType(t,'ZoneBC_t')
        C._rmBCOfType(t,'BCdimNS')
        C._rmBCOfType(t,'BCReconsLBM')

        C._addPeriodicZones__(t[2][1])
        for zn in dictOfTempPerNodes:
            TempNode = dictOfTempPerNodes[zn]
            z = Internal.getNodeFromName2(t,zn)
            z[2].append(TempNode)
        C._rmBCOfType(t,"BCMatch")

    t = X.connectMatch(t)
    #Internal._addGhostCells(t, t, 2)
    Internal._addGhostCells(t, t, 2, adaptBCs=1, fillCorner=1)
    C._rmBCOfType(t,"BCMatch")
    C._fillEmptyBCWith(t,'overlap','BCOverlap',dim=3)
    X._applyBCOverlaps(t, depth=2, loc='centers')

    tc = C.node2Center(t)
    #
    if isPerio: C._removeDuplicatedPeriodicZones__(t)
    #
    tc = X.setInterpData(t,tc,nature=1, loc='centers', storage='inverse',sameName=1,itype='chimera')
    #tc = X.setInterpData(t,tc,nature=1, loc='centers', storage='inverse',sameName=1,dim=3)

    # on remet chez la zone initiale les donnees d interp
    if isPerio:
        for z in Internal.getZones(tc):
            if z[0] not in zonesNamesORIG:
                parentz,noz = Internal.getParentOfNode(tc,z)
                zdnames = z[0].split('_')
                znameorig = zdnames[0]
                zorig = Internal.getNodeFromNameAndType(tc,znameorig,'Zone_t')
                IDS = Internal.getNodesFromType(z,'ZoneSubRegion_t')
                for ID in IDS:
                    ID_orig = Internal.getNodeFromNameAndType(zorig,ID[0],'ZoneSubRegion_t')
                    if ID_orig is None:
                        zorig[2].append(ID)
                    else:
                        PL_NEW  = Internal.getNodeFromName(ID,'PointList')
                        PL_ORIG = Internal.getNodeFromName(ID_orig,'PointList')
                        PL_ORIG[1] = numpy.concatenate((PL_ORIG[1],PL_NEW[1]))

                        PLD_NEW = Internal.getNodeFromName(ID,'PointListDonor')
                        PLD_ORIG = Internal.getNodeFromName(ID_orig,'PointListDonor')
                        PLD_ORIG[1] = numpy.concatenate((PLD_ORIG[1],PLD_NEW[1]))

                        COEF_NEW = Internal.getNodeFromName(ID,'InterpolantsDonor')
                        COEF_ORIG = Internal.getNodeFromName(ID_orig,'InterpolantsDonor')
                        COEF_ORIG[1] = numpy.concatenate((COEF_ORIG[1],COEF_NEW[1]))

                        TYPE_NEW = Internal.getNodeFromName(ID,'InterpolantsType')
                        TYPE_ORIG = Internal.getNodeFromName(ID_orig,'InterpolantsType')
                        TYPE_ORIG[1] = numpy.concatenate((TYPE_ORIG[1],TYPE_NEW[1]))

                del parentz[2][noz]

    C._rmBCOfType(t,"BCOverlap")

    count_NSLBM = []
    for i,z in enumerate(Internal.getZones(tc)):
        count_NSLBM.append(0)
        IDS = Internal.getNodesFromType(z,'ZoneSubRegion_t')
        for ID in IDS:
            NSLBM  = Internal.getNodeFromName(ID,'NSLBM')
            if NSLBM is not None:
                count_NSLBM[i] = count_NSLBM[i] + 1
    print(count_NSLBM)

    for i,z in enumerate(Internal.getZones(t)):
        dim = Internal.getZoneDim(z)
        dim = dim[1:4]
        if infos_couplage_zones[z[0]] != []:
            for bc in infos_couplage_zones[z[0]]:
                bc_name = Internal.getName(bc)
                print('Copying BC {} to zone {}.'.format(bc_name,z[0]))
                bc_type = Internal.getValue(bc)
                node = Internal.getNodeFromName(bc, 'PointRange')
                bc_range = Internal.getValue(node)
                bc_range = sum(bc_range.tolist(),[])
                print(nb_vert[i])
                print(bc_range)
                for k in range(3):
                    if(bc_range[2*k] == nb_vert[i][k]) : bc_range[2*k] = dim[k]
                    if(bc_range[2*k+1] == nb_vert[i][k]) : bc_range[2*k+1] = dim[k]
                print(bc_range)
                C._addBC2Zone(z,bc_name,bc_type,bc_range)
                count_NSLBM[i] = count_NSLBM[i] - 1
    #
    # EN COURS DE DVPT
    #
    bc_coins = {}
    dim_zones = {}
    for i,z in enumerate(Internal.getZones(tc)):
        bc_coins[z[0]] = []
        dim = Internal.getZoneDim(z); dim = dim[1:4]
        dim_zones[z[0]] = dim

    for i,z in enumerate(Internal.getZones(tc)):
        if z[0] in zones_NS and count_NSLBM[i]!=0:
            #Dans ce cas il reste des raccords non traites : coins
            IDS = Internal.getNodesFromType(z,'ZoneSubRegion_t')
            for ID in IDS:
                node_NSLBM  = Internal.getNodeFromName(ID,'NSLBM')
                if node_NSLBM is not None:
                    zone_dest = Internal.getValue(ID)
                    print(z[0], zone_dest)
                    PL_ORIG = Internal.getNodeFromName(ID,'PointListDonor')
                    list_pl = Internal.getValue(PL_ORIG)
                    #test des coins a revoir
                    coin = 2*2*dim_zones[zone_dest][-1]
                    if len(list_pl) == coin:
                        k = list_pl//((dim_zones[zone_dest][0])*(dim_zones[zone_dest][1]))
                        j = (list_pl-k*(dim_zones[zone_dest][0])*(dim_zones[zone_dest][1]))//(dim_zones[zone_dest][0])
                        i = list_pl - j*(dim_zones[zone_dest][0]) - k*dim_zones[zone_dest][0]*dim_zones[zone_dest][1]
                        bc_range = [i.min()+1,i.max()+1,j.min()+1,j.max()+1,k.min()+1,k.max()+1]
                        print(bc_range)
                        bc_coins[zone_dest].append(bc_range)

    for i,z in enumerate(Internal.getZones(t)):
        if bc_coins[z[0]]!=[]:
            print('BC COINS', z[0])
            bc_name = 'adim_coin'
            bc_type = 'BCadimcoins'
            dim = Internal.getZoneDim(z)
            dim = dim[1:4]
            for k in range(len(bc_coins[z[0]])):
                if bc_coins[z[0]][k][0:4]==[1,2,1,2]:
                    #Coin en bas a droite : dans ce cas imin
                    bc_range = [1,1,1,dim[1],1,dim[-1]]
                elif bc_coins[z[0]][k][0:4]==[dim[0]-2,dim[0]-1,1,2]:
                    #Coin en bas a gauche : dans ce cas jmin
                    bc_range = [1,dim[0],1,1,1,dim[-1]]
                elif bc_coins[z[0]][k][0:4]==[1,2,dim[1]-2,dim[1]-1]:
                    #Coin en haut a droite : dans ce cas jmax
                    bc_range = [1,dim[0],dim[1],dim[1],1,dim[-1]]
                elif bc_coins[z[0]][k][0:4]==[dim[0]-2,dim[0]-1,dim[1]-2,dim[1]-1]:
                    #Coin en haut a gauche : dans ce cas imax
                    bc_range = [dim[0],dim[0],1,dim[1],1,dim[-1]]
                print(bc_range)
                C._addBC2Zone(z,bc_name,bc_type,bc_range)

    Internal._rmNodesFromType(tc,"GridCoordinates_t")

    if isinstance(tc_out, str):
        C.convertPyTree2File(tc, tc_out)
    if isinstance(t_out, str):
        C.convertPyTree2File(t, t_out)

    return t, tc

#===============================================================================
class CLBMNS(Common):
    """ Preparation pour les simulations couplees FastLBM - FastS """
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["alexandre.suss@onera.fr"]
        self.cartesian = True

    #Prepare
    def prepare(self, t_case, t_out=None, tc_out=None, NP=0, solvertype=None,translation=[0.,0.,0.], NG=2):
        #if NP is None: NP = Cmpi.size
        #if NP == 0: print('Preparing for a sequential computation.')
        #else: print('Preparing for a computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, solvertype=solvertype, translation=translation, NP=NP, format=self.data['format'], NG=NG)
        return ret
