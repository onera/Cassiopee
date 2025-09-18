# - gestion du compactage et des transferts compacts -
from . import connector
import numpy

import Converter.Internal as Internal

from . import OversetData as XOD

ibm_lbm_variables_1 ='Q_'
ibm_lbm_variables_2 ='Qstar_'
ibm_lbm_variables_3 ='Qneq_'
NEQ_LBM =  89

#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la base
# fonctionne avec ___setInterpTransfer
#==============================================================================
def miseAPlatDonorTree__(t, tc, graph=None, list_graph=None, nbpts_linelets=0):
    if isinstance(graph, list):
        ###########################IMPORTANT ######################################
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        graphliste=True
    else:
        graphliste=False

    import Converter.Mpi as Cmpi
    rank = Cmpi.rank

    zones = Internal.getZones(t)

    if graph is not None and graphliste==False:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
        if 'graphID_Unsteady' in graph:
            graphID_U = graph['graphID_Unsteady']
            graphID_S = graph['graphID_Steady']
        else:
            graphID_U = None; graphID_S = None
    elif graph is not None and graphliste==True:
        procDict  = graph[0]['procDict']
        graphID   = graph[0]['graphID']
        graphIBCD = graph[0]['graphIBCD']
        graphID_U = None; graphID_S = None
    else:
        procDict=None; graphID=None; graphIBCD=None; graphID_U = None; graphID_S = None

    size_int  = 0
    size_real = 0
    listproc  = []
    rac       = []
    rac_inst  = []
    sizeI     = []
    sizeR     = []
    sizeNb    = []
    sizeNbD   = []
    sizeNbFlu = []
    sizeType  = []
    nrac      = 0

    ordered_subRegions=[]
    neq_subRegions=[]
    No_zoneD=[]
    MeshTypeD=[]
    infoZoneList={}
    inst={}
    numero_max =-100000000
    numero_min = 100000000

    ntab_int    = 18
    sizefluxCons=7

    bases = Internal.getNodesFromType1(tc, 'CGNSBase_t')  # noeud
    c     = 0
    for base in bases:

        model    = 'NSLaminar'
        a        = Internal.getNodeFromName2(base, 'GoverningEquations')
        if a is not None: model = Internal.getValue(a)
        if model=="NSLaminar" or model=="Euler": neq_base=5
        elif model=="NSTurbulent": neq_base=6
        elif model=='LBMLaminar':
            neq_base = Internal.getNodeFromName2(zones[0] , 'Parameter_int')[1][NEQ_LBM]

        zones_tc = Internal.getZones(base)
        for z in zones_tc:

            model_z = None
            a = Internal.getNodeFromName2(z, 'GoverningEquations')
            if a is not None:
                model_z = Internal.getValue(a)
                if model_z=="NSLaminar" or model_z=="Euler": neq_trans=5
                elif model_z=="NSTurbulent": neq_trans=6
                elif model_z=="LBMLaminar":
                    neq_trans = Internal.getNodeFromName2(zones[c], 'Parameter_int')[1][NEQ_LBM]

            else: neq_trans = neq_base

            subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
            meshtype   = 1
            zonetype   = Internal.getNodeFromType1(z, 'ZoneType_t')
            tmp        = Internal.getValue(zonetype)
            if tmp != "Structured": meshtype = 2
            infoZoneList[zones[c][0]] = [ c , None ] # init par defaut
            for s in subRegions:
                zRname = Internal.getValue(s)
                proc = 0
                if procDict is not None: proc = procDict[zRname]

                #tri des pas de temps instationnaire
                #  1) les stationnaires
                #  2) les instationnaires regroupes par pas de temps
                if '#' not in s[0]:
                    ordered_subRegions.append(s)
                    neq_subRegions.append(neq_trans)
                    No_zoneD.append(c)
                    MeshTypeD.append(meshtype)
                    #print('RANk=',rank, 'NODonneuse=', c, s[0], z[0])
                else:
                    numero_iter = int( s[0].split('#')[1].split('_')[0] )
                    if numero_iter < numero_min : numero_min = numero_iter
                    if numero_iter > numero_max : numero_max = numero_iter

                    if numero_iter in inst:
                        sub = inst[ numero_iter ][0]
                        sub = sub + [s]
                        Noz = inst[ numero_iter ][1]
                        Noz = Noz + [c]
                        mesh= inst[ numero_iter ][2]
                        mesh= mesh+ [meshtype]
                        dest= inst[ numero_iter ][3]
                        dest= dest+ [proc]
                        neqtrans= inst[ numero_iter ][4]
                        neqtrans= neqtrans+ [neq_trans]
                        inst[ numero_iter ]= [ sub , Noz , mesh, dest, neqtrans ]
                    else:
                        inst[ numero_iter ]= [ [s],[c],[meshtype], [proc], [neq_trans] ]

                TimeLevelNumber = len(inst)

                if TimeLevelNumber != 1+numero_max-numero_min and len(inst) != 0:
                    raise ValueError("miseAPlatDonorTree__: missing timestep in tc : %d %d %d")%(numero_max,numero_min, TimeLevelNumber)


                count_ID  = 0
                count_IBC = 0
                # alloc memoire
                pointlist     =  Internal.getNodeFromName1(s, 'PointList')
                pointlistD    =  Internal.getNodeFromName1(s, 'PointListDonor')
                InterpD       =  Internal.getNodeFromName1(s, 'InterpolantsDonor')
                Interptype    =  Internal.getNodeFromName1(s, 'InterpolantsType')
                RotationAngle =  Internal.getNodeFromName1(s, 'RotationAngle')
                RotationCenter=  Internal.getNodeFromName1(s, 'RotationCenter')
                prange        =  Internal.getNodeFromName1(s, 'PointRange')      # Besoin des point range pour l'explicite local
                pranged       =  Internal.getNodeFromName1(s, 'PointRangeDonor') # Besoin des point range pour l'explicite local
                direction     =  Internal.getNodeFromName1(s, 'DirReceveur')     # Besoin des directions pour l'explicite local
                directiond    =  Internal.getNodeFromName1(s, 'DirDonneur')      # Besoin des point directions pour l'explicite local
                transfo       =  Internal.getNodeFromName1(s, 'Transform')       # Besoin du transform pour l'explicite local
                pt_pivot      =  Internal.getNodeFromName1(s, 'PointPivot')      # Besoin du point pivot pour l'explicite local (conservativite)
                profondeur    =  Internal.getNodeFromName1(s, 'Profondeur')      # Besoin de la profondeur pour l'explicite local (nearmatch)
                ratio         =  Internal.getNodeFromName1(s, 'NMratio')         # Besoin des ratios entre les pas d espace des zones donneuse et receveuse (exp local)
                levelrcv      =  Internal.getNodeFromName1(s, 'LevelZRcv')       # Niveau en temps zone receveuse (exp local)
                leveldnr      =  Internal.getNodeFromName1(s, 'LevelZDnr')       # Niveau en temps zone donneuse (exp local)

                utau          =  Internal.getNodeFromName1(s, 'utau')
                temp_local    =  Internal.getNodeFromName1(s, 'Temperature')
                wmodel_local  =  Internal.getNodeFromName1(s, 'Density_WM')
                gradxP        =  Internal.getNodeFromName1(s, 'gradxPressure')
                gradxU        =  Internal.getNodeFromName1(s, 'gradxVelocityX')
                xcInit        =  Internal.getNodeFromName1(s, 'CoordinateX_PC#Init')
                motion_type   =  Internal.getNodeFromName1(s, 'MotionType')
                kcurv         =  Internal.getNodeFromName1(s, XOD.__KCURV__)
                sd1           =  Internal.getNodeFromName1(s, 'StagnationEnthalpy')
                yline         =  Internal.getNodeFromName1(s, 'CoordinateN_ODE')

                Nbpts         =  numpy.shape(pointlist[ 1])[0]
                Nbpts_D       =  numpy.shape(pointlistD[1])[0]
                Nbpts_InterpD =  numpy.shape(InterpD[ 1  ])[0]

                sname         = s[0][0:2]

                # recherche raccord conservatif
                zd = zones[c]
                nbfluCons = 0
                bcs = Internal.getNodesFromType2(zd, 'BC_t')
                bc_info = []
                for c4, bc in enumerate(bcs):
                    bcname = bc[0].split('_')
                    btype = Internal.getValue(bc)
                    if 'BCFluxOctreeF' == btype and bcname[1]==zRname: nbfluCons+=1 #flux donneur(Fine)
                    if 'BCFluxOctreeC' == btype:                                    #flux receveur (Coarse)
                        idir_tg=-10
                        if   bcname[2][0:4]=='imin': idir_tg=2
                        elif bcname[2][0:4]=='imax': idir_tg=1
                        elif bcname[2][0:4]=='jmin': idir_tg=4
                        elif bcname[2][0:4]=='jmax': idir_tg=3
                        elif bcname[2][0:4]=='kmin': idir_tg=6
                        elif bcname[2][0:4]=='kmax': idir_tg=5
                        bc_info.append( [ bcname[1], idir_tg, c4 ] )
                infoZoneList[zd[0]] = [ c , bc_info ] #[ No zoneD, bcs(conservatif)]

                if nbfluCons !=0: print("flux cons: Z_D/R",  z[0], zRname, 'nflux=', nbfluCons)

                sname = s[0][0:2]
                # cas ou les vitesses n'ont pas ete ajoutees lors du prep (ancien tc)
                if sname == 'IB':
                    vx = Internal.getNodeFromName1(s, 'VelocityX')
                    if vx is None:
                        density = Internal.getNodeFromName1(s, 'Density')
                        nIBC    = density[1].shape[0]
                        vxNP = numpy.zeros((nIBC),numpy.float64)
                        vyNP = numpy.zeros((nIBC),numpy.float64)
                        vzNP = numpy.zeros((nIBC),numpy.float64)
                        s[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
                        s[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
                        s[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])
                    if model == "LBMLaminar":
                        density = Internal.getNodeFromName1(s, 'Density')
                        nIBC    = density[1].shape[0]

                        qloc_1  = Internal.getNodeFromName1(s, ibm_lbm_variables_1 + str(1))
                        if qloc_1 is None:
                            for f_i in range (1,neq_trans+1):
                                s[2].append([ibm_lbm_variables_1 + str(f_i) , numpy.zeros((nIBC),numpy.float64) , [], 'DataArray_t'])
                        qloc_1  = Internal.getNodeFromName1(s, ibm_lbm_variables_1 + str(1))

                        qloc_2 = Internal.getNodeFromName1(s, ibm_lbm_variables_2 + str(1))
                        if qloc_2 is None:
                            for f_i in range (1,neq_trans+1):
                                s[2].append([ibm_lbm_variables_2 + str(f_i) , numpy.zeros((nIBC),numpy.float64) , [], 'DataArray_t'])
                        qloc_2 = Internal.getNodeFromName1(s, ibm_lbm_variables_2 + str(1))

                        qloc_3 = Internal.getNodeFromName1(s, ibm_lbm_variables_3 + str(1))
                        if qloc_3 is None:
                            for f_i in range (1,neq_trans+1):
                                s[2].append([ibm_lbm_variables_3 + str(f_i) , numpy.zeros((nIBC),numpy.float64) , [], 'DataArray_t'])
                        qloc_3 = Internal.getNodeFromName1(s, ibm_lbm_variables_3 + str(1))
                # DBX: supp utau en euler
                #if utau is None:
                #   utauNP  = numpy.zeros((Nbpts_D),numpy.float64)
                #   yplusNP = numpy.zeros((Nbpts_D),numpy.float64)
                #   Internal.createUniqueChild(s, 'utau'  , 'DataArray_t', utauNP )
                #   Internal.createUniqueChild(s, 'yplus' , 'DataArray_t', yplusNP )
                #   utau =  Internal.getNodeFromName1(s, 'utau')

                # on recupere le nombre de type different
                #typecell = Interptype[1][0]
                #Nbtype = [ typecell ]
                #for i in range(Nbpts_D):
                #  if Interptype[1][i] not in Nbtype: Nbtype.append(Interptype[1][i])
                #print('nb type',  len(Nbtype), s[0],z[0], Nbtype)
                nbType = numpy.unique(Interptype[1])
                nbTypeSize = nbType.size

                size_IBC =  0
                ntab_IBC = 11+3 #On ajoute dorenavant les vitesses dans l'arbre tc pour faciliter le post
                if utau is not None: ntab_IBC += 2
                if temp_local is not None: ntab_IBC += 2
                if wmodel_local is not None: ntab_IBC += 6
                if gradxP is not None: ntab_IBC += 3
                if gradxU is not None: ntab_IBC += 9
                if xcInit is not None: ntab_IBC += 9 # 3 for each type IBM point - 3 wall points, 3 target points, & 3 image points
                if motion_type is not None: ntab_IBC += 11 #MotionType,transl_speed(x3),axis_pnt(x3),axis_vct(x3),omega
                if kcurv is not None: ntab_IBC += 1
                if sd1 is not None: ntab_IBC += 5
                if yline is not None: ntab_IBC += (7*nbpts_linelets+2)
                if sname == 'IB' and model == "LBMLaminar":
                    if qloc_1 is not None: ntab_IBC += neq_trans
                    if qloc_2 is not None: ntab_IBC += neq_trans
                    if qloc_3 is not None: ntab_IBC += neq_trans
                if sname == 'IB':
                    size_IBC   = Nbpts_D*ntab_IBC
                    count_IBC += 1
                else:
                    count_ID  += 1

                # periodicite azimutale
                rotation = 0
                if RotationAngle is not None: rotation +=3
                if RotationCenter is not None: rotation +=3

                nrac  =  nrac + 1
                if proc not in listproc:
                    listproc.append(proc)
                    rac.append(1)
                    if '#' in s[0]: rac_inst.append(1)
                    else          : rac_inst.append(0)
                    sizeI.append(    Nbpts_D*2     + Nbpts   + nbTypeSize+1 + nbfluCons*sizefluxCons)
                    sizeR.append(    Nbpts_InterpD + size_IBC + rotation                            )
                    sizeNbD.append(  Nbpts_D                                                        )
                    sizeNb.append(   Nbpts                                                          )
                    sizeType.append( Nbpts_D                 + nbTypeSize+1                         )
                    sizeNbFlu.append( nbfluCons*sizefluxCons)
                else:
                    pos           = listproc.index(proc)
                    rac[pos]      = rac[pos] + 1
                    if '#' in s[0]: rac_inst[pos]= rac_inst[pos] + 1
                    sizeI[pos]    = sizeI[pos]     + Nbpts_D*2     + Nbpts   + nbTypeSize+1  + nbfluCons*sizefluxCons
                    sizeR[pos]    = sizeR[pos]     + Nbpts_InterpD + size_IBC + rotation
                    sizeNbD[pos]  = sizeNbD[pos]   + Nbpts_D
                    sizeNb[pos]   = sizeNb[pos]    + Nbpts
                    sizeType[pos] = sizeType[pos]  + Nbpts_D                 + nbTypeSize+1
                    sizeNbFlu[pos]= sizeNbFlu[pos] + nbfluCons*sizefluxCons
            c += 1

    base     = Internal.getNodeFromType1(tc, 'CGNSBase_t')  # noeud
    model    = 'NSLaminar'
    a        = Internal.getNodeFromName2(base, 'GoverningEquations')
    if a is not None: model = Internal.getValue(a)

    NbP2P     = len(listproc)  #nombre Comm MPI point a Point pour envoi
    sizeproc  = []
    ntab_tot  = ntab_int

    if nrac != 0 and prange is not None: ntab_tot += 27

    for i in range(NbP2P): sizeproc.append(4 + TimeLevelNumber*2 + ntab_tot*rac[i] + sizeI[i])

    size_int  =  2 + NbP2P + sum(sizeproc)
    size_real =  sum(sizeR)


    if not graphliste: # Si le graph n est pas une liste, on n'est pas en explicite local
        #on determine la liste des processus pour lequel rank  est Receveur
        graphIBCrcv=[];graphIDrcv=[]
        if graphIBCD is not None:
            #on recupere les infos Steady
            graphIBCrcv_=[]; pos_IBC=[]; S_IBC= 1; graphloc=[]
            S_IBC = _procSource(rank, S_IBC,  pos_IBC, graphIBCD, graphloc, graphIBCrcv_)

            graphIBCrcv  = pos_IBC + graphIBCrcv_

        if graphID_U is not None:
            #on recupere les infos Steady
            graphIDrcv_=[];graphrcv_S=[]; pos_ID=[]; S_ID=TimeLevelNumber + 1
            S_ID = _procSource(rank, S_ID, pos_ID, graphID_S, graphrcv_S, graphIDrcv_)
            #on ajoute les infos UNsteady
            for nstep in range(numero_min,numero_max+1):
                graphloc=[]
                S_ID = _procSource(rank, S_ID, pos_ID, graphID_U[nstep], graphloc, graphIDrcv_, filterGraph=graphrcv_S)

            graphIDrcv   = pos_ID  + graphIDrcv_

        else:
            #on recupere les infos ID Steady
            graphIDrcv_=[];graphloc=[]; pos_ID=[]; S_ID=1
            if graphID is not None:
                S_ID = _procSource(rank, S_ID, pos_ID, graphID, graphloc, graphIDrcv_)

                graphIDrcv = pos_ID + graphIDrcv_

    else:  # le graph est une liste, on est en explicite local, 1 graphe par ss ite

        graphIBCrcv_=[]; graphIDrcv_=[]; pos_ID=[]; pos_IBC=[]; S_IBC=len(graph); S_ID=len(graph)
        for nstep in range(0,len(graph)):

            graphIBCD_= graph[nstep]['graphIBCD']
            stokproc  = []
            S_IBC     = _procSource(rank, S_IBC, pos_IBC, graphIBCD_, stokproc, graphIBCrcv_)

            graphID_  = graph[nstep]['graphID']
            stokproc  = []
            S_ID      = _procSource(rank, S_ID, pos_ID, graphID_, stokproc, graphIDrcv_)

        graphIBCrcv  = pos_IBC + graphIBCrcv_
        graphIDrcv   = pos_ID  + graphIDrcv_


    ## Recuperation info zone
    infoZoneList = Cmpi.allgatherDict(infoZoneList)

    #print("len graphIBCrcv is",len(graphIBCrcv))
    #print("len graphIDrcv  is",len(graphIDrcv))
    #print("pos_IBC is",pos_IBC)
    #print("pos_ID  is",pos_ID)

    npass_transfert=1  # valeur par defaut. Pourra valoir 2 si raccord ordre 4 ou 5

    param_int  = numpy.empty(size_int + len(graphIDrcv) + len(graphIBCrcv) + 2, dtype=Internal.E_NpyInt)
    param_real = numpy.empty(size_real, dtype=numpy.float64)
    Internal.createUniqueChild(tc, 'Parameter_int' , 'DataArray_t', param_int)
    if size_real != 0:
        Internal.createUniqueChild(tc, 'Parameter_real', 'DataArray_t', param_real)

    if len(graphIBCrcv) == 0:
        _graphIBC = numpy.zeros(1, dtype=Internal.E_NpyInt)
    else:
        _graphIBC = numpy.asarray([len(graphIBCrcv)]+graphIBCrcv, dtype=Internal.E_NpyInt)

    _graphID = numpy.asarray([len(graphIDrcv)] +graphIDrcv, dtype=Internal.E_NpyInt)

    param_int[2                 :3+len(graphIBCrcv)                ] = _graphIBC
    param_int[3+len(graphIBCrcv):4+len(graphIBCrcv)+len(graphIDrcv)] = _graphID

    # print("param_int is ",param_int[0:2+len(graphIBCrcv)+len(graphIDrcv)+1])

    #
    #initialisation numpy
    #
    param_int[0] = 0  #flag pour init transfert couche C (X. Juvigny)
    param_int[1] = NbP2P
    size_ptflux = []
    size_ptlist = []
    size_ptlistD= []
    size_ptType = []
    nb_rac      = []
    size_coef   = []
    adr_coef    = []   # pour cibler debut de echange dans param_real

    shift_graph = len(graphIDrcv) + len(graphIBCrcv) + 4
    # print("shift_graph is ",shift_graph)
    shift_coef  =0
    shift       = shift_graph # le shift prend en compte la postion des graphs (ID+IBC) entre la address contenant NbP2P et
    for i in range(NbP2P):
        adr_coef.append(shift_coef)                    #adresse echange dans param_real
        shift_coef = shift_coef + sizeR[i]

        param_int[i+shift_graph] = NbP2P + shift              #adresse echange
        shift =  shift  + sizeproc[i]
        size_ptflux.append(0)
        size_ptlist.append(0)
        size_ptlistD.append(0)
        size_ptType.append(0)
        size_coef.append(0)
        nb_rac.append(0)

    for iter in range(numero_min, numero_max+1):
        ordered_subRegions =  ordered_subRegions + inst[ iter ][0]
        No_zoneD           =  No_zoneD           + inst[ iter ][1]
        MeshTypeD          =  MeshTypeD          + inst[ iter ][2]
        neq_subRegions     =  neq_subRegions     + inst[ iter ][4]

    # loop sur les raccords tries
    c        = 0
    Nbtot    = 0
    S = 0
    for s in ordered_subRegions:

        zRname   = Internal.getValue(s)

        NozoneD  = No_zoneD[c]
        meshtype = MeshTypeD[c]
        neq_loc  = neq_subRegions[c]

        #recherche raccord conservatif
        zd = zones[NozoneD]
        nbfluCons=0
        bcs = Internal.getNodesFromType2(zd, 'BC_t')
        no_bc=[]
        for c1, bc in enumerate(bcs):
            bcname = bc[0].split('_')
            btype = Internal.getValue(bc)
            if 'BCFluxOctreeF' == btype and bcname[1] == zRname:
                nbfluCons += 1
                no_bc.append(c1)
                param_int_zD = Internal.getNodeFromName2( zd, 'Parameter_int' )[1]
                #print("bc",  s[0], bc[0])

        #if nbfluCons !=0: print("Verif flux cons: Z_D/R",  zd[0], zRname, 'nflux=', nbfluCons, 'NoZoneD:', NozoneD)

        proc = 0
        if procDict is not None: proc = procDict[zRname]
        pos  = listproc.index(proc)
        pt_ech = param_int[ pos + shift_graph]      # adresse debut raccord pour l'echange pos
        pt_coef= adr_coef[pos] + size_coef[pos]     # adresse debut coef

        pointlist     =  Internal.getNodeFromName1(s, 'PointList')
        pointlistD    =  Internal.getNodeFromName1(s, 'PointListDonor')
        Interptype    =  Internal.getNodeFromName1(s, 'InterpolantsType')
        InterpD       =  Internal.getNodeFromName1(s, 'InterpolantsDonor')
        RotationAngle =  Internal.getNodeFromName1(s, 'RotationAngle')
        RotationCenter=  Internal.getNodeFromName1(s, 'RotationCenter')
        prange        =  Internal.getNodeFromName1(s, 'PointRange')       # Besoin des point range pour l'explicite local
        pranged       =  Internal.getNodeFromName1(s, 'PointRangeDonor')  # Besoin des point range pour l'explicite local
        direction     =  Internal.getNodeFromName1(s, 'DirReceveur')      # Besoin des directions pour l'explicite local
        directiond    =  Internal.getNodeFromName1(s, 'DirDonneur')  # Besoin des point directions pour l'explicite local
        transfo       =  Internal.getNodeFromName1(s, 'Transform')  # Besoin du transform pour l'explicite local(conservativite)
        pt_pivot      =  Internal.getNodeFromName1(s, 'PointPivot')  # Besoin du point pivot pour l'explicite local (conservativite)
        profondeur    =  Internal.getNodeFromName1(s, 'Profondeur')  # Besoin de la profondeur pour l'explicite local (nearmatch)
        ratio         =  Internal.getNodeFromName1(s, 'NMratio') # Besoin des ratios entre les pas d espace des zones donneuse et receveuse (exp local)
        levelrcv      =  Internal.getNodeFromName1(s, 'LevelZRcv') # Niveau en temps zone receveuse (exp local)
        leveldnr      =  Internal.getNodeFromName1(s, 'LevelZDnr') # Niveau en temps zone donneuse (exp local)

        #print(zRname, nb_rac[pos])

        Nbpts         =  numpy.shape(pointlist[ 1])[0]
        Nbpts_D       =  numpy.shape(pointlistD[1])[0]
        Nbpts_InterpD =  numpy.shape(InterpD[ 1  ])[0]

        param_int[ pt_ech    ] = proc                    #proc destination
        param_int[ pt_ech +1 ] = rac[pos]                #nbr de raccord a envoyer a proc destination
        param_int[ pt_ech +2 ] = rac_inst[pos]           #nbr de raccord instat a envoyer a proc destination
        nrac_steady            = rac[pos] - rac_inst[pos]

        param_int[ pt_ech +3 ] = TimeLevelNumber         #nbr pas de temps instationnaire dans raccord
        nrac_inst_deb  =  nrac_steady
        for i in range(TimeLevelNumber):

            # len(inst[i][3])  = list destination du rac pour le temps i
            NracInsta=0
            for procSearch in inst[i+numero_min][3]:
                if procSearch==proc: NracInsta+=1

            nrac_inst_fin  = nrac_inst_deb + NracInsta

            #print('NracInsta=',NracInsta,'TimeLevel=',i, 'dest=',proc)

            param_int[ pt_ech +4 + i                  ] = nrac_inst_deb   #unsteadyJoins(No rac debut)
            param_int[ pt_ech +4 + i + TimeLevelNumber] = nrac_inst_fin   #unsteadyJoins(No rac fin)

            nrac_inst_deb  = nrac_inst_fin

        iadr2    = pt_ech + 4 + TimeLevelNumber*2
        iadr     = iadr2 + nb_rac[pos]             # ptr echange + dest + nrac ...
        iadr_ptr = iadr2 + ntab_tot*rac[pos]

        param_int[ iadr            ] = Nbpts
        param_int[ iadr + rac[pos] ] = Nbpts_InterpD

        #on recupere le nombre de type different
        #NbType = numpy.unique(Interptype[1])
        #NbTypeSize = NbType.size
        typecell = Interptype[1][0]
        Nbtype= [ typecell ]
        for i in range(Nbpts_D):
            if Interptype[1][i] not in Nbtype: Nbtype.append(Interptype[1][i])

        #Si le type zero existe, on le place a la fin: sinon adressage openmp=boom dans donorPts
        if 0 in Nbtype: Nbtype += [Nbtype.pop( Nbtype.index( 0 ) )]

        param_int[ iadr+rac[pos]*2 ] = len(Nbtype)

        param_int[ iadr+rac[pos]*3 ] = -1
        size_IBC = 0

        zsrname = s[0]
        sname = zsrname[0:2]
        xc=None;yc=None;zc=None; xi=None;yi=None;zi=None; xw=None;yw=None;zw=None;density=None;pressure=None
        vx=None; vy=None; vz=None
        utau=None;yplus=None;kcurv=None
        ptxc=0;ptyc=0;ptzc=0;ptxi=0;ptyi=0;ptzi=0;ptxw=0;ptyw=0;ptzw=0;ptdensity=0;ptpressure=0
        ptvx=0;ptvy=0;ptvz=0
        ptutau=0;ptyplus=0;ptkcurv=None
        yline=None; uline=None; nutildeline=None; psiline=None; matmline=None; matline=None; matpline=None
        alphasbetaline=None; indexline=None
        ptyline=0; ptuline=0; ptnutildeline=0; ptpsiline=0; ptmatmline=0; ptmatline=0; ptmatpline=0
        ptalphasbetaline=0; ptindexline=0

        temp_local=None;pttemp_local=0;
        temp_extra_local =None;pttemp_extra_local =0;
        temp_extra2_local=None;pttemp_extra2_local=0;

        wmodel_local_dens=None;wmodel_local_velx=None;wmodel_local_vely=None;
        wmodel_local_velz=None;wmodel_local_temp=None;wmodel_local_sanu=None;
        pt_dens_wm_local=0;pt_velx_wm_local=0;pt_vely_wm_local=0;
        pt_velz_wm_local=0;pt_temp_wm_local=0;pt_sanu_wm_local=0;

        gradxP=None;gradyP=None;gradzP=None;
        ptgradxP=0;ptgradyP=0;ptgradzP=0;

        gradxU=None;gradyU=None;gradzU=None;
        ptgradxU=0;ptgradyU=0;ptgradzU=0;

        gradxV=None;gradyV=None;gradzV=None;
        ptgradxV=0;ptgradyV=0;ptgradzV=0;

        gradxW=None;gradyW=None;gradzW=None;
        ptgradxW=0;ptgradyW=0;ptgradzW=0;

        xcInit=None;ycInit=None;zcInit=None;
        xiInit=None;yiInit=None;ziInit=None;
        xwInit=None;ywInit=None;zwInit=None;
        ptxcInit=0;ptycInit=0;ptzcInit=0;
        ptxiInit=0;ptyiInit=0;ptziInit=0;
        ptxwInit=0;ptywInit=0;ptzwInit=0;

        motion_type=None;
        transl_speedX=None;transl_speedY=None;transl_speedZ=None;
        axis_pntX=None;axis_pntY=None;axis_pntZ=None;
        axis_vctX=None;axis_vctY=None;axis_vctZ=None;
        omega=None;
        ptmotion_type=0;
        pttransl_speedX=0;pttransl_speedY=0;pttransl_speedZ=0;
        ptaxis_pntX=0;ptaxis_pntY=0;ptaxis_pntZ=0;
        ptaxis_vctX=0;ptaxis_vctY=0;ptaxis_vctZ=0;
        ptomega=0;

        sd1=None;sd2=None;sd3=None;sd4=None;sd5=None
        ptd1=0;ptd2=0;ptd3=0;ptd4=0;ptd5=0

        qloc_1=[None]*(neq_loc)
        ptqloc_1=[0]*(neq_loc)

        qloc_2=[None]*(neq_loc)
        ptqloc_2=[0]*(neq_loc)

        qloc_3=[None]*(neq_loc)
        ptqloc_3=[0]*(neq_loc)

        if sname == 'IB':
            zsrname = zsrname.split('_')
            if len(zsrname) < 3:
                #print('Warning: miseAPlatDonorTree: non consistent with the version of IBM preprocessing.')
                if model=='Euler':
                    print('Assuming IBC type is wallslip.')
                    param_int[iadr+rac[pos]*3]  = 0
                elif model=='LBMLaminar':
                    print('Assuming IBC type is no-slip.')
                    param_int[iadr+rac[pos]*3]  = 1
                else:
                    print('Assuming IBC type is Musker wall model.')
                    param_int[iadr+rac[pos]*3]  = 3
            else:
                if "Mobile" in  zsrname[2]:
                    param_int[iadr+rac[pos]*3]  = 7  # musker paroi en rotation
                else:
                    param_int[iadr+rac[pos]*3]  = int(zsrname[1]) # 'IBCD_type_zonename'

            IBCType = param_int[iadr+rac[pos]*3]
            #print('len zsrname', len(zsrname),param_int[iadr+rac[pos]*3] )

#           print('IBCType = ', IBCType)
            xc        = Internal.getNodeFromName1(s , 'CoordinateX_PC')
            yc        = Internal.getNodeFromName1(s , 'CoordinateY_PC')
            zc        = Internal.getNodeFromName1(s , 'CoordinateZ_PC')
            xi        = Internal.getNodeFromName1(s , 'CoordinateX_PI')
            yi        = Internal.getNodeFromName1(s , 'CoordinateY_PI')
            zi        = Internal.getNodeFromName1(s , 'CoordinateZ_PI')
            xw        = Internal.getNodeFromName1(s , 'CoordinateX_PW')
            yw        = Internal.getNodeFromName1(s , 'CoordinateY_PW')
            zw        = Internal.getNodeFromName1(s , 'CoordinateZ_PW')

            density   = Internal.getNodeFromName1(s , 'Density')
            pressure  = Internal.getNodeFromName1(s , 'Pressure')

            ptxc      = pt_coef + Nbpts_InterpD
            ptyc      = pt_coef + Nbpts_InterpD + Nbpts_D
            ptzc      = pt_coef + Nbpts_InterpD + Nbpts_D*2
            ptxi      = pt_coef + Nbpts_InterpD + Nbpts_D*3
            ptyi      = pt_coef + Nbpts_InterpD + Nbpts_D*4
            ptzi      = pt_coef + Nbpts_InterpD + Nbpts_D*5
            ptxw      = pt_coef + Nbpts_InterpD + Nbpts_D*6
            ptyw      = pt_coef + Nbpts_InterpD + Nbpts_D*7
            ptzw      = pt_coef + Nbpts_InterpD + Nbpts_D*8
            ptdensity = pt_coef + Nbpts_InterpD + Nbpts_D*9
            ptpressure= pt_coef + Nbpts_InterpD + Nbpts_D*10

            size_IBC  = 11*Nbpts_D
            inc = 11

            vx = Internal.getNodeFromName1(s, 'VelocityX')
            vy = Internal.getNodeFromName1(s, 'VelocityY')
            vz = Internal.getNodeFromName1(s, 'VelocityZ')
            ptvx = pt_coef + Nbpts_InterpD + Nbpts_D*inc
            ptvy = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
            ptvz = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
            size_IBC +=3*Nbpts_D; inc += 3

            utau = Internal.getNodeFromName1(s, 'utau')
            yplus = Internal.getNodeFromName1(s, 'yplus')
            if utau is not None:
                ptutau    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptyplus   = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                size_IBC += 2*Nbpts_D; inc += 2

            kcurv = Internal.getNodeFromName1(s, XOD.__KCURV__)
            if kcurv is not None:
                ptkcurv = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC  += Nbpts_D; inc += 1

            temp_local      = Internal.getNodeFromName1(s, 'Temperature')
            if temp_local is not None:
                pttemp_local = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC   += Nbpts_D; inc += 1

            temp_extra_local = Internal.getNodeFromName1(s, 'TemperatureWall')
            if temp_extra_local is not None:
                pttemp_extra_local = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC   += Nbpts_D; inc += 1

            temp_extra2_local = Internal.getNodeFromName1(s, 'WallHeatFlux')
            if temp_extra2_local is not None:
                pttemp_extra2_local = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC   += Nbpts_D; inc += 1

            wmodel_local_dens = Internal.getNodeFromName1(s, 'Density_WM')
            wmodel_local_velx = Internal.getNodeFromName1(s, 'VelocityX_WM')
            wmodel_local_vely = Internal.getNodeFromName1(s, 'VelocityY_WM')
            wmodel_local_velz = Internal.getNodeFromName1(s, 'VelocityZ_WM')
            wmodel_local_temp = Internal.getNodeFromName1(s, 'Temperature_WM')
            wmodel_local_sanu = Internal.getNodeFromName1(s, 'TurbulentSANuTilde_WM')
            if wmodel_local_dens is not None:
                pt_dens_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                pt_velx_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                pt_vely_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                pt_velz_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+3)
                pt_temp_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+4)
                pt_sanu_wm_local    = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+5)
                size_IBC   += 6*Nbpts_D; inc += 6


            gradxP = Internal.getNodeFromName1(s, 'gradxPressure')
            gradyP = Internal.getNodeFromName1(s, 'gradyPressure')
            gradzP = Internal.getNodeFromName1(s, 'gradzPressure')
            if gradxP is not None:
                ptgradxP = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptgradyP = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                ptgradzP = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                size_IBC   += 3*Nbpts_D; inc += 3

            gradxU = Internal.getNodeFromName1(s, 'gradxVelocityX')
            gradyU = Internal.getNodeFromName1(s, 'gradyVelocityX')
            gradzU = Internal.getNodeFromName1(s, 'gradzVelocityX')
            if gradxU is not None:
                ptgradxU = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptgradyU = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                ptgradzU = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                size_IBC   += 3*Nbpts_D; inc += 3

            gradxV = Internal.getNodeFromName1(s, 'gradxVelocityY')
            gradyV = Internal.getNodeFromName1(s, 'gradyVelocityY')
            gradzV = Internal.getNodeFromName1(s, 'gradzVelocityY')
            if gradxV is not None:
                ptgradxV = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptgradyV = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                ptgradzV = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                size_IBC   += 3*Nbpts_D; inc += 3

            gradxW = Internal.getNodeFromName1(s, 'gradxVelocityZ')
            gradyW = Internal.getNodeFromName1(s, 'gradyVelocityZ')
            gradzW = Internal.getNodeFromName1(s, 'gradzVelocityZ')
            if gradxW is not None:
                ptgradxW = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptgradyW = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                ptgradzW = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                size_IBC   += 3*Nbpts_D; inc += 3


            xcInit    = Internal.getNodeFromName1(s , 'CoordinateX_PC#Init')
            ycInit    = Internal.getNodeFromName1(s , 'CoordinateY_PC#Init')
            zcInit    = Internal.getNodeFromName1(s , 'CoordinateZ_PC#Init')
            xiInit    = Internal.getNodeFromName1(s , 'CoordinateX_PI#Init')
            yiInit    = Internal.getNodeFromName1(s , 'CoordinateY_PI#Init')
            ziInit    = Internal.getNodeFromName1(s , 'CoordinateZ_PI#Init')
            xwInit    = Internal.getNodeFromName1(s , 'CoordinateX_PW#Init')
            ywInit    = Internal.getNodeFromName1(s , 'CoordinateY_PW#Init')
            zwInit    = Internal.getNodeFromName1(s , 'CoordinateZ_PW#Init')
            if xcInit is not None:
                ptxcInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+0)
                ptycInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                ptzcInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)

                ptxiInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+3)
                ptyiInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+4)
                ptziInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+5)

                ptxwInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+6)
                ptywInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+7)
                ptzwInit = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+8)
                size_IBC   += 9*Nbpts_D; inc += 9

            motion_type  = Internal.getNodeFromName1(s , 'MotionType')
            transl_speedX = Internal.getNodeFromName1(s , 'transl_speedX')
            transl_speedY = Internal.getNodeFromName1(s , 'transl_speedY')
            transl_speedZ = Internal.getNodeFromName1(s , 'transl_speedZ')
            axis_pntX     = Internal.getNodeFromName1(s , 'axis_pntX')
            axis_pntY     = Internal.getNodeFromName1(s , 'axis_pntY')
            axis_pntZ     = Internal.getNodeFromName1(s , 'axis_pntZ')
            axis_vctX     = Internal.getNodeFromName1(s , 'axis_vctX')
            axis_vctY     = Internal.getNodeFromName1(s , 'axis_vctY')
            axis_vctZ     = Internal.getNodeFromName1(s , 'axis_vctZ')
            omega        = Internal.getNodeFromName1(s , 'omega')
            if motion_type is not None:
                ptmotion_type   = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+0)
                pttransl_speedX = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                pttransl_speedY = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+2)
                pttransl_speedZ = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+3)
                ptaxis_pntX     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+4)
                ptaxis_pntY     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+5)
                ptaxis_pntZ     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+6)
                ptaxis_vctX     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+7)
                ptaxis_vctY     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+8)
                ptaxis_vctZ     = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+9)
                ptomega         = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+10)
                size_IBC   += 11*Nbpts_D; inc += 11

            sd1 = Internal.getNodeFromName1(s, 'StagnationEnthalpy')
            if sd1 is not None:
                ptd1    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC += Nbpts_D; inc += 1
            sd2 = Internal.getNodeFromName1(s, 'StagnationPressure')
            if sd2 is not None:
                ptd2    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC += Nbpts_D; inc += 1
            sd3 = Internal.getNodeFromName1(s, 'dirx')
            if sd3 is not None:
                ptd3    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC += Nbpts_D; inc += 1
            sd4 = Internal.getNodeFromName1(s, 'diry')
            if sd4 is not None:
                ptd4    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC += Nbpts_D; inc += 1
            sd5 = Internal.getNodeFromName1(s, 'dirz')
            if sd5 is not None:
                ptd5    = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                size_IBC += Nbpts_D; inc += 1

            if model=="LBMLaminar":
                qloc_1[0] = Internal.getNodeFromName1(s, ibm_lbm_variables_1 + str(1))
                if qloc_1[0] is not None:
                    ptqloc_1[0] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                    size_IBC += Nbpts_D; inc += 1
                    for f_i in range(1,neq_loc):
                        qloc_1[f_i]   = Internal.getNodeFromName1(s, ibm_lbm_variables_1 + str(f_i+1))
                        ptqloc_1[f_i] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                        size_IBC += Nbpts_D; inc += 1

                qloc_2[0] = Internal.getNodeFromName1(s, ibm_lbm_variables_2 + str(1))
                if qloc_2[0] is not None:
                    ptqloc_2[0] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                    size_IBC += Nbpts_D; inc += 1
                    for f_i in range(1,neq_loc):
                        qloc_2[f_i]   = Internal.getNodeFromName1(s, ibm_lbm_variables_2 + str(f_i+1))
                        ptqloc_2[f_i] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                        size_IBC += Nbpts_D; inc += 1

                qloc_3[0] = Internal.getNodeFromName1(s, ibm_lbm_variables_3 + str(1))
                if qloc_3[0] is not None:
                    ptqloc_3[0] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                    size_IBC += Nbpts_D; inc += 1
                    for f_i in range(1,neq_loc):
                        qloc_3[f_i]   = Internal.getNodeFromName1(s, ibm_lbm_variables_3 + str(f_i+1))
                        ptqloc_3[f_i] = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                        size_IBC += Nbpts_D; inc += 1

            yline          = Internal.getNodeFromName1(s, 'CoordinateN_ODE')
            uline          = Internal.getNodeFromName1(s, 'VelocityT_ODE')
            nutildeline    = Internal.getNodeFromName1(s, 'TurbulentSANuTilde_ODE')
            psiline        = Internal.getNodeFromName1(s, 'Psi_ODE')
            matmline       = Internal.getNodeFromName1(s, 'Matm_ODE')
            matline        = Internal.getNodeFromName1(s, 'Mat_ODE')
            matpline       = Internal.getNodeFromName1(s, 'Matp_ODE')
            alphasbetaline = Internal.getNodeFromName1(s, 'alphasbeta_ODE')
            indexline      = Internal.getNodeFromName1(s, 'index_ODE')
            if yline is not None:
                ptyline          = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptuline          = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptnutildeline    = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptpsiline        = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptmatmline       = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptmatline        = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                ptmatpline       = pt_coef + Nbpts_InterpD + Nbpts_D*inc; inc += nbpts_linelets
                size_IBC += 7*nbpts_linelets*Nbpts_D;
                ptalphasbetaline = pt_coef + Nbpts_InterpD + Nbpts_D*inc
                ptindexline      = pt_coef + Nbpts_InterpD + Nbpts_D*(inc+1)
                size_IBC += 2*Nbpts_D; inc += 2


        tmp = Internal.getNodeFromName1(s, 'ZoneRole')
        if tmp[1][0] == b'D': param_int[ iadr +rac[pos]*4 ] = 0   # role= Donor
        else                : param_int[ iadr +rac[pos]*4 ] = 1   # role= Receiver

        param_int[ iadr +rac[pos]*5 ] = NozoneD                    # No zone donneuse

        lstD  = iadr_ptr + size_ptlistD[pos]
        ptTy  = iadr_ptr + sizeNbD[pos] + size_ptType[pos]
        lst   = iadr_ptr + sizeNbD[pos] + sizeType[pos] + size_ptlist[pos]
        ptFlux= iadr_ptr + sizeNbD[pos] + sizeType[pos] + sizeNb[pos] + size_ptflux[pos]
        param_int[ iadr +rac[pos]*6 ] = lst                                                                      # PointlistAdr
        param_int[ iadr +rac[pos]*7 ] = ptTy                                                                     # TypAdr
        param_int[ iadr +rac[pos]*12] = lstD                                                                     # PointlistDAdr
        param_int[ iadr +rac[pos]*17] = ptFlux                                                                   # Ptr_flux

        #
        #construction info flux conservatif
        #
        zd = zones[NozoneD]
        nbfluCons=0
        bcs   = Internal.getNodesFromType2(zd, 'BC_t')
        nb_bc = len(bcs)
        no_bc = []
        for c1, bc in enumerate(bcs):
            bcname = bc[0].split('_')
            btype = Internal.getValue(bc)
            if 'BCFluxOctreeF' == btype and bcname[1] == zRname:
                nbfluCons+=1
                no_bc.append(c1)
                param_int_zD  = Internal.getNodeFromName2( zd, 'Parameter_int' )[1]
                param_real_zD = Internal.getNodeFromName2( zd, 'Parameter_real' )[1]

        for c2 , no in enumerate(no_bc):
            pt_bcs = param_int_zD[70]  #70=PT_BC
            pt_bc  = param_int_zD[pt_bcs + 1 + no]
            idir_bc= param_int_zD[pt_bc + 1]
            sz_bc  =(param_int_zD[pt_bc +3]-param_int_zD[pt_bc+2]+1)*(param_int_zD[pt_bc +5]-param_int_zD[pt_bc+4]+1)*(param_int_zD[pt_bc +7]-param_int_zD[pt_bc+6]+1)
            nbdata = param_int_zD[pt_bc + 8]
            iptsize_data = pt_bc + 9
            ipt_data = param_int_zD[pt_bcs + 1 + no + nb_bc]

            #print("nbdata: ", nbdata, "sz0", param_int_zD[iptsize_data], "sz1", param_int_zD[iptsize_data+1], "ratio",  param_real_zD[ipt_data], param_real_zD[ipt_data+2], param_real_zD[ipt_data+3])

            pt_flu = ptFlux + c2*sizefluxCons
            param_int[ pt_flu   ]= idir_bc
            param_int[ pt_flu +1]= sz_bc
            param_int[ pt_flu +2]= no

            #determination No Bc sur zone receuveuse a partir du nom du noeud bc
            for bc in infoZoneList[zRname][1]:
                if bc[0] == zd[0] and bc[1] == idir_bc: param_int[pt_flu+3] = bc[2]

            param_int[ pt_flu +4]= int( param_real_zD[ipt_data  ] )
            param_int[ pt_flu +5]= int( param_real_zD[ipt_data+1] )
            param_int[ pt_flu +6]= int( param_real_zD[ipt_data+2] )

            print("Flux compact: zr=", zRname, "zd=", zd[0], 'idir', idir_bc, 'sz_bc', sz_bc, 'ptbc', pt_bc, 'NobcR', param_int[pt_flu +3], 'ratioK:', param_int[ pt_flu +6] )
        #
        #Fin flux conservatif
        #

        Nbtot += Nbpts

        param_int[ ptTy ] = len(Nbtype)
        noi       = 0
        nocoef    = 0
        sizecoef  = 0
        shift_typ = 1 + len(Nbtype)
        ctyp      = 0
        l0        = 0

        #recopie dans tableau a plat + tri par type
        if len(Nbtype) == 1:
            npass = triMonoType(Nbpts_D, Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                                Interptype, pointlist, pointlistD, param_int,
                                ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw,
                                ptdensity,ptpressure, ptkcurv,
                                ptvx, ptvy, ptvz,
                                ptutau,ptyplus,
                                pttemp_local, pttemp_extra_local,pttemp_extra2_local,
                                pt_dens_wm_local,pt_velx_wm_local,pt_vely_wm_local,
                                pt_velz_wm_local,pt_temp_wm_local,pt_sanu_wm_local,
                                ptgradxP, ptgradyP, ptgradzP,
                                ptgradxU, ptgradyU, ptgradzU,
                                ptgradxV, ptgradyV, ptgradzV,
                                ptgradxW, ptgradyW, ptgradzW,
                                ptxcInit,ptycInit,ptzcInit,ptxiInit,ptyiInit,ptziInit,ptxwInit,ptywInit,ptzwInit,
                                ptmotion_type,
                                pttransl_speedX,pttransl_speedY,pttransl_speedZ,
                                ptaxis_pntX,ptaxis_pntY,ptaxis_pntZ,
                                ptaxis_vctX,ptaxis_vctY,ptaxis_vctZ,
                                ptomega,
                                ptd1,ptd2,ptd3,ptd4,ptd5,
                                ptyline,ptuline,ptnutildeline,ptpsiline,ptmatmline,ptmatline,ptmatpline,
                                ptalphasbetaline,ptindexline,
                                xc,yc,zc,xi,yi,zi,xw,yw,zw,
                                density,pressure, kcurv,
                                vx, vy, vz,
                                utau, yplus,
                                temp_local, temp_extra_local, temp_extra2_local,
                                wmodel_local_dens,wmodel_local_velx,wmodel_local_vely,
                                wmodel_local_velz,wmodel_local_temp,wmodel_local_sanu,
                                gradxP, gradyP, gradzP,
                                gradxU, gradyU, gradzU,
                                gradxV, gradyV, gradzV,
                                gradxW, gradyW, gradzW,
                                xcInit,ycInit,zcInit,xiInit,yiInit,ziInit,xwInit,ywInit,zwInit,
                                motion_type,
                                transl_speedX,transl_speedY,transl_speedZ,
                                axis_pntX,axis_pntY,axis_pntZ,
                                axis_vctX,axis_vctY,axis_vctZ,
                                omega,
                                sd1,sd2,sd3,sd4,sd5,
                                yline,uline,nutildeline,psiline,matmline,matline,matpline,
                                alphasbetaline,indexline,
                                InterpD,param_real,ptqloc_1,qloc_1,ptqloc_2,qloc_2,ptqloc_3,qloc_3,neq_loc,model,nbpts_linelets)

        else:
            npass = triMultiType(Nbpts_D,Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                                 Interptype, pointlist, pointlistD, param_int,
                                 ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw,
                                 ptdensity,ptpressure, ptkcurv,
                                 ptvx, ptvy, ptvz,
                                 ptutau,ptyplus,
                                 pttemp_local, pttemp_extra_local,pttemp_extra2_local,
                                 pt_dens_wm_local,pt_velx_wm_local,pt_vely_wm_local,
                                 pt_velz_wm_local,pt_temp_wm_local,pt_sanu_wm_local,
                                 ptgradxP, ptgradyP, ptgradzP,
                                 ptgradxU, ptgradyU, ptgradzU,
                                 ptgradxV, ptgradyV, ptgradzV,
                                 ptgradxW, ptgradyW, ptgradzW,
                                 ptxcInit,ptycInit,ptzcInit,ptxiInit,ptyiInit,ptziInit,ptxwInit,ptywInit,ptzwInit,
                                 ptmotion_type,
                                 pttransl_speedX,pttransl_speedY,pttransl_speedZ,
                                 ptaxis_pntX,ptaxis_pntY,ptaxis_pntZ,
                                 ptaxis_vctX,ptaxis_vctY,ptaxis_vctZ,
                                 ptomega,
                                 ptd1,ptd2,ptd3,ptd4,ptd5,
                                 ptyline,ptuline,ptnutildeline,ptpsiline,ptmatmline,ptmatline,ptmatpline,
                                 ptalphasbetaline,ptindexline,
                                 xc,yc,zc,xi,yi,zi,xw,yw,zw,
                                 density,pressure, kcurv,
                                 vx, vy, vz,
                                 utau,yplus,
                                 temp_local, temp_extra_local, temp_extra2_local,
                                 wmodel_local_dens,wmodel_local_velx,wmodel_local_vely,
                                 wmodel_local_velz,wmodel_local_temp,wmodel_local_sanu,
                                 gradxP, gradyP, gradzP,
                                 gradxU, gradyU, gradzU,
                                 gradxV, gradyV, gradzV,
                                 gradxW, gradyW, gradzW,
                                 xcInit,ycInit,zcInit,xiInit,yiInit,ziInit,xwInit,ywInit,zwInit,
                                 motion_type,
                                 transl_speedX,transl_speedY,transl_speedZ,
                                 axis_pntX,axis_pntY,axis_pntZ,
                                 axis_vctX,axis_vctY,axis_vctZ,
                                 omega,
                                 sd1,sd2,sd3,sd4,sd5,
                                 yline,uline,nutildeline,psiline,matmline,matline,matpline,
                                 alphasbetaline,indexline,
                                 InterpD,param_real,ptqloc_1,qloc_1,ptqloc_2,qloc_2,ptqloc_3,qloc_3,neq_loc,model)

        if npass==2: npass_transfert=2

        pointlist[ 1] = param_int[ lst             : lst              + Nbpts         ]    # supression numpy initial pointlist
        Interptype[1] = param_int[ ptTy + shift_typ: ptTy + shift_typ + Nbpts_D       ]    # supression numpy initial interpolantType
        pointlistD[1] = param_int[ lstD            : lstD             + Nbpts_D       ]    # supression numpy initial pointlistDonor
        InterpD[   1] = param_real[ pt_coef        : pt_coef          + Nbpts_InterpD ]    # supression numpy initial interpDonor

        #if s[0] == 'ID_cart3' and z[0]=='cart1': print('verif',  InterpD[   1][0], pt_coef,numpy.shape(InterpD[ 1  ]))

        if sname == 'IB':
            xc[1]       = param_real[ ptxc: ptxc+ Nbpts_D ]
            yc[1]       = param_real[ ptyc: ptyc+ Nbpts_D ]
            zc[1]       = param_real[ ptzc: ptzc+ Nbpts_D ]
            xi[1]       = param_real[ ptxi: ptxi+ Nbpts_D ]
            yi[1]       = param_real[ ptyi: ptyi+ Nbpts_D ]
            zi[1]       = param_real[ ptzi: ptzi+ Nbpts_D ]                                      # supression numpy initial IBC
            xw[1]       = param_real[ ptxw: ptxw+ Nbpts_D ]
            yw[1]       = param_real[ ptyw: ptyw+ Nbpts_D ]
            zw[1]       = param_real[ ptzw: ptzw+ Nbpts_D ]
            density[1]  = param_real[ ptdensity : ptdensity + Nbpts_D ]
            pressure[1] = param_real[ ptpressure: ptpressure+ Nbpts_D ]

            vx[1]       = param_real[ ptvx: ptvx+ Nbpts_D ]
            vy[1]       = param_real[ ptvy: ptvy+ Nbpts_D ]
            vz[1]       = param_real[ ptvz: ptvz+ Nbpts_D ]

            if utau is not None:
                utau[1]  = param_real[ ptutau : ptutau + Nbpts_D ]
                yplus[1] = param_real[ ptyplus: ptyplus + Nbpts_D ]

            if temp_local is not None:
                temp_local[1]       = param_real[ pttemp_local       : pttemp_local       + Nbpts_D ]
            if temp_extra_local is not None:
                temp_extra_local[1] = param_real[ pttemp_extra_local : pttemp_extra_local + Nbpts_D ]
            if temp_extra2_local is not None:
                temp_extra2_local[1] = param_real[ pttemp_extra2_local: pttemp_extra2_local + Nbpts_D ]

            if wmodel_local_dens is not None:
                wmodel_local_dens[1]       = param_real[ pt_dens_wm_local    : pt_dens_wm_local    + Nbpts_D ]
                wmodel_local_velx[1]       = param_real[ pt_velx_wm_local    : pt_velx_wm_local    + Nbpts_D ]
                wmodel_local_vely[1]       = param_real[ pt_vely_wm_local    : pt_vely_wm_local    + Nbpts_D ]
                wmodel_local_velz[1]       = param_real[ pt_velz_wm_local    : pt_velz_wm_local    + Nbpts_D ]
                wmodel_local_temp[1]       = param_real[ pt_temp_wm_local    : pt_temp_wm_local    + Nbpts_D ]
                wmodel_local_sanu[1]       = param_real[ pt_sanu_wm_local    : pt_sanu_wm_local    + Nbpts_D ]


            if gradxP is not None:
                gradxP[1] = param_real[ ptgradxP : ptgradxP + Nbpts_D ]
                gradyP[1] = param_real[ ptgradyP : ptgradyP + Nbpts_D ]
                gradzP[1] = param_real[ ptgradzP : ptgradzP + Nbpts_D ]

            if gradxU is not None:
                gradxU[1] = param_real[ ptgradxU : ptgradxU + Nbpts_D ]
                gradyU[1] = param_real[ ptgradyU : ptgradyU + Nbpts_D ]
                gradzU[1] = param_real[ ptgradzU : ptgradzU + Nbpts_D ]

            if gradxV is not None:
                gradxV[1] = param_real[ ptgradxV : ptgradxV + Nbpts_D ]
                gradyV[1] = param_real[ ptgradyV : ptgradyV + Nbpts_D ]
                gradzV[1] = param_real[ ptgradzV : ptgradzV + Nbpts_D ]

            if gradxW is not None:
                gradxW[1] = param_real[ ptgradxW : ptgradxW + Nbpts_D ]
                gradyW[1] = param_real[ ptgradyW : ptgradyW + Nbpts_D ]
                gradzW[1] = param_real[ ptgradzW : ptgradzW + Nbpts_D ]

            if xcInit is not None:
                xcInit = param_real[ ptxcInit : ptxcInit + Nbpts_D ]
                ycInit = param_real[ ptycInit : ptycInit + Nbpts_D ]
                zcInit = param_real[ ptzcInit : ptzcInit + Nbpts_D ]

                xiInit = param_real[ ptxiInit : ptxiInit + Nbpts_D ]
                yiInit = param_real[ ptyiInit : ptyiInit + Nbpts_D ]
                ziInit = param_real[ ptziInit : ptziInit + Nbpts_D ]

                xwInit = param_real[ ptxwInit : ptxwInit + Nbpts_D ]
                ywInit = param_real[ ptywInit : ptywInit + Nbpts_D ]
                zwInit = param_real[ ptzwInit : ptzwInit + Nbpts_D ]

            if motion_type is not None:
                motion_type  = param_real[ ptmotion_type : ptmotion_type + Nbpts_D ]

                transl_speedX = param_real[ pttransl_speedX : pttransl_speedX + Nbpts_D ]
                transl_speedY = param_real[ pttransl_speedY : pttransl_speedY + Nbpts_D ]
                transl_speedZ = param_real[ pttransl_speedZ : pttransl_speedZ + Nbpts_D ]

                axis_pntX     = param_real[ ptaxis_pntX : ptaxis_pntX + Nbpts_D ]
                axis_pntY     = param_real[ ptaxis_pntY : ptaxis_pntY + Nbpts_D ]
                axis_pntZ     = param_real[ ptaxis_pntZ : ptaxis_pntZ + Nbpts_D ]

                axis_vctX     = param_real[ ptaxis_vctX : ptaxis_vctX + Nbpts_D ]
                axis_vctY     = param_real[ ptaxis_vctY : ptaxis_vctY + Nbpts_D ]
                axis_vctZ     = param_real[ ptaxis_vctZ : ptaxis_vctZ + Nbpts_D ]

                omega        = param_real[ ptomega : ptomega + Nbpts_D ]

            if kcurv is not None:
                kcurv[1]  = param_real[ ptkcurv : ptkcurv + Nbpts_D ]

            if sd1 is not None: sd1[1]  = param_real[ ptd1 : ptd1 + Nbpts_D ]
            if sd2 is not None: sd2[1]  = param_real[ ptd2 : ptd2 + Nbpts_D ]
            if sd3 is not None: sd3[1]  = param_real[ ptd3 : ptd3 + Nbpts_D ]
            if sd4 is not None: sd4[1]  = param_real[ ptd4 : ptd4 + Nbpts_D ]
            if sd5 is not None: sd5[1]  = param_real[ ptd5 : ptd5 + Nbpts_D ]


            if yline is not None:
                yline[1]          = param_real[ ptyline : ptyline + Nbpts_D*nbpts_linelets ]
                uline[1]          = param_real[ ptuline : ptuline + Nbpts_D*nbpts_linelets ]
                nutildeline[1]    = param_real[ ptnutildeline : ptnutildeline + Nbpts_D*nbpts_linelets ]
                psiline[1]        = param_real[ ptpsiline : ptpsiline + Nbpts_D*nbpts_linelets ]
                matmline[1]       = param_real[ ptmatmline : ptmatmline + Nbpts_D*nbpts_linelets ]
                matline[1]        = param_real[ ptmatline : ptmatline + Nbpts_D*nbpts_linelets ]
                matpline[1]       = param_real[ ptmatpline : ptmatpline + Nbpts_D*nbpts_linelets ]
                alphasbetaline[1] = param_real[ ptalphasbetaline : ptalphasbetaline + Nbpts_D ]
                indexline[1]      = param_real[ ptindexline : ptindexline + Nbpts_D ]

            if model=="LBMLaminar":
                if qloc_1[0] is not None:
                    for f_i in range(0,neq_loc):
                        qloc_1[f_i][1]   = param_real[ ptqloc_1[f_i] : ptqloc_1[f_i] + Nbpts_D ]
                if qloc_2[0] is not None:
                    for f_i in range(0,neq_loc):
                        qloc_2[f_i][1]   = param_real[ ptqloc_2[f_i] : ptqloc_2[f_i] + Nbpts_D ]
                if qloc_3[0] is not None:
                    for f_i in range(0,neq_loc):
                        qloc_3[f_i][1]   = param_real[ ptqloc_3[f_i] : ptqloc_3[f_i] + Nbpts_D ]

        param_int[ iadr +rac[pos]*8 ] = adr_coef[pos] + size_coef[pos]          # PtcoefAdr

        tmp = Internal.getNodeFromName1(s , 'GridLocation')
        if tmp[1][4] == b'C': param_int[ iadr +rac[pos]*9 ] = 1   # location= CellCenter
        else                : param_int[ iadr +rac[pos]*9 ] = 0   # location= Vertex

        param_int[ iadr +rac[pos]*10 ] = Nbpts_D

        NoR = infoZoneList[zRname][0]

        param_int[ iadr +rac[pos]*11  ]= NoR  # No zone receveuse

        #print('rac',s[0],'zoneR=',zRname,'NoR=',NoR ,'adr=',iadr +rac[pos]*11, 'NoD=', param_int[ iadr-1 +rac[pos]*5 ], 'adr=',iadr-1 +rac[pos]*5,'rank=',rank,'dest=', proc)

        #print 'model=',model,'zoneR',zones_tc[param_int[ iadr +rac[pos]*11  ]][0], 'NoR=', NoR, 'NoD=', c
        tmp =  Internal.getNodeFromName1(s , 'RANSLES')
        if tmp is not None: param_int[ iadr +rac[pos]*13  ] = min (5, neq_loc)   # RANSLES
        else:               param_int[ iadr +rac[pos]*13  ] = neq_loc

        tmp =  Internal.getNodeFromName1(s , 'NSLBM')
        if tmp is not None:
            if neq_loc==5: param_int[ iadr +rac[pos]*13  ] = 11  # NS vers LBM
            else:          param_int[ iadr +rac[pos]*13  ] = -5  # LBM vers NS

        #print('raccord',s[0], 'neq=',neq_loc, param_int[ iadr +rac[pos]*13  ] , 'NoDR=',NozoneD, param_int[ iadr +rac[pos]*11  ])

        # raccord periodique avec rotation
        if RotationAngle is not None:
            param_int[ iadr +rac[pos]*14  ] = 1
            shiftRotation                   = 6
            ptdeb =   pt_coef + Nbpts_InterpD
            param_real[ ptdeb   : ptdeb+3 ] = RotationAngle[1][0:3]
            param_real[ ptdeb+3 : ptdeb+6 ] = RotationCenter[1][0:3]
            RotationAngle[1]  =    param_real[ ptdeb   : ptdeb+3]
            RotationCenter[1] =    param_real[ ptdeb+3 : ptdeb+6]

        else:
            param_int[ iadr +rac[pos]*14  ] = 0
            shiftRotation                   = 0

        # raccord instationnaire
        #Ivan: rac*15: inutile?? On affecte levelReceveur pour preparer nettoyage dtloc
        param_int[ iadr + rac[pos]*15 ] = 1  # valeur par defaut niveau temps zone receveuse
        if levelrcv is not None:
            param_int[ iadr +rac[pos]*15 ] = int(levelrcv[1])

        # nombre de flux conservatif a corriger
        param_int[ iadr + rac[pos]*16 ] = nbfluCons

        '''
        param_int[ iadr +rac[pos]*15  ] =-1
        if '#' in s[0]:
            numero_iter = int( s[0].split('#')[1].split('_')[0] )
            param_int[ iadr +rac[pos]*15  ] = numero_iter
        '''

        if nrac != 0:

            if pranged is not None and prange is not None : #Si on est en explicite local, on rajoute des choses dans param_int

                adr = iadr2 + rac[pos]*ntab_int + 27*nb_rac[pos]
                param_int[ adr     : adr+ 6 ] = prange[1]
                param_int[ adr+ 7  : adr+ 13] = pranged[1]
                param_int[ adr+ 14 : adr+ 17] = transfo[1]
                param_int[ adr+ 17 : adr+ 20] = pt_pivot[1]
                param_int[ adr+ 20 : adr+ 21] = profondeur[1]
                param_int[ adr+ 21 : adr+ 24] = ratio[1]
                param_int[ adr+ 6]  = direction[1]
                param_int[ adr+ 13] = directiond[1]
                param_int[ adr+ 24] = int(levelrcv[1])
                param_int[ adr+ 25] = int(leveldnr[1])
                param_int[ adr+ 26] = S

                if int(levelrcv[1]) != int(leveldnr[1]):
                    S=S + 3*5*(pranged[1][1]-pranged[1][0]+1)*(pranged[1][3]-pranged[1][2]+1)*(pranged[1][5]-pranged[1][4]+1)*(profondeur[1][0]+1)

        ### Verifier position et choix entre Nbpts et NbptsD
        size_ptlistD[pos] = size_ptlistD[pos] + Nbpts_D
        size_ptlist[pos]  = size_ptlist[pos]  + Nbpts
        size_ptType[pos]  = size_ptType[pos]  + Nbpts_D + len(Nbtype)+1
        size_ptflux[pos]  = size_ptflux[pos]  + nbfluCons*sizefluxCons

        size_coef[pos] = size_coef[pos] + Nbpts_InterpD + size_IBC + shiftRotation
        nb_rac[pos]    = nb_rac[pos] + 1

        c += 1

    tmp       = Internal.getNodeFromName1(t, '.Solver#ownData')
    dtloc     = Internal.getNodeFromName1(tmp, '.Solver#dtloc')
    if dtloc is not None: dtloc[1][12] = npass_transfert

    return None

#==============================================================================
# tri multitype
#==============================================================================
def triMultiType(Nbpts_D, Nbpts, Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                 Interptype, pointlist, pointlistD, param_int,
                 ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw,
                 ptdensity,ptpressure, ptkcurv,
                 ptvx, ptvy, ptvz,
                 ptutau,ptyplus,
                 pttemp_local, pttemp_extra_local, pttemp_extra2_local,
                 pt_dens_wm_local,pt_velx_wm_local,pt_vely_wm_local,
                 pt_velz_wm_local,pt_temp_wm_local,pt_sanu_wm_local,
                 ptgradxP, ptgradyP, ptgradzP,
                 ptgradxU, ptgradyU, ptgradzU,
                 ptgradxV, ptgradyV, ptgradzV,
                 ptgradxW, ptgradyW, ptgradzW,
                 ptxcInit,ptycInit,ptzcInit,ptxiInit,ptyiInit,ptziInit,ptxwInit,ptywInit,ptzwInit,
                 ptmotion_type,
                 pttransl_speedX,pttransl_speedY,pttransl_speedZ,
                 ptaxis_pntX,ptaxis_pntY,ptaxis_pntZ,
                 ptaxis_vctX,ptaxis_vctY,ptaxis_vctZ,
                 ptomega,
                 ptd1,ptd2,ptd3,ptd4,ptd5,
                 ptyline,ptuline,ptnutildeline,ptpsiline,ptmatmline,ptmatline,ptmatpline,
                 ptalphasbetaline,ptindexline,
                 xc,yc,zc,xi,yi,zi,xw,yw,zw,
                 density,pressure,kcurv,
                 vx, vy, vz,
                 utau,yplus,
                 temp_local, temp_extra_local, temp_extra2_local,
                 wmodel_local_dens,wmodel_local_velx,wmodel_local_vely,
                 wmodel_local_velz,wmodel_local_temp,wmodel_local_sanu,
                 gradxP, gradyP, gradzP,
                 gradxU, gradyU, gradzU,
                 gradxV, gradyV, gradzV,
                 gradxW, gradyW, gradzW,
                 xcInit,ycInit,zcInit,xiInit,yiInit,ziInit,xwInit,ywInit,zwInit,
                 motion_type,
                 transl_speedX,transl_speedY,transl_speedZ,
                 axis_pntX,axis_pntY,axis_pntZ,
                 axis_vctX,axis_vctY,axis_vctZ,
                 omega,
                 sd1,sd2,sd3,sd4,sd5,
                 yline,uline,nutildeline,psiline,matmline,matline,matpline,
                 alphasbetaline,indexline,
                 InterpD,param_real,ptqloc_1,qloc_1,ptqloc_2,qloc_2,ptqloc_3,qloc_3,neq_loc,model):

    npass_transfert = 1

    for ntype in Nbtype:
        noi_old   = 0
        nocoef_old= 0
        l         = 0
        for i in range(Nbpts_D):
            ltype = Interptype[1][i]
            if meshtype == 1:
                if ltype == 1: sizecoef=1
                elif ltype == 2: sizecoef=8
                elif ltype ==44:
                    sizecoef        =12
                    npass_transfert = 2
                elif ltype == 3: sizecoef=9
                elif ltype == 5:
                    sizecoef        =15
                    npass_transfert = 2
                elif ltype == 22: sizecoef=4
                elif ltype == 4: sizecoef=8

            else:
                if ltype == 1: sizecoef=1
                elif ltype == 4: sizecoef=4

            if ltype == ntype:
                # recopie interpolantType
                param_int[ ptTy + shift_typ + l + l0 ] = ltype

                # recopie pointlist
                if ntype != 0:
                    param_int[  lst + noi] = pointlist[ 1][ noi_old ]
                    noi     = noi     + 1
                    noi_old = noi_old + 1
                else:
                    ncfLoc   = pointlist[ 1][ noi_old ]
                    sizecoef = ncfLoc
                    param_int[  lst + noi] = ncfLoc
                    param_int[ lst+noi+1: lst+noi+1+ncfLoc] = pointlist[1][ noi_old+1: noi_old+1+ncfLoc]
                    noi     = noi     + 1 + ncfLoc
                    noi_old = noi_old + 1 + ncfLoc

                # recopie pointListDonor
                param_int[ lstD +  l + l0] = pointlistD[1][i]
                # recopie Maillage IBC
                if sname == 'IB':
                    param_real[ ptxc      + l + l0 ]= xc[1][i]
                    param_real[ ptyc      + l + l0 ]= yc[1][i]
                    param_real[ ptzc      + l + l0 ]= zc[1][i]
                    param_real[ ptxi      + l + l0 ]= xi[1][i]
                    param_real[ ptyi      + l + l0 ]= yi[1][i]
                    param_real[ ptzi      + l + l0 ]= zi[1][i]
                    param_real[ ptxw      + l + l0 ]= xw[1][i]
                    param_real[ ptyw      + l + l0 ]= yw[1][i]
                    param_real[ ptzw      + l + l0 ]= zw[1][i]
                    param_real[ ptdensity + l + l0 ]= density[1][i]
                    param_real[ ptpressure+ l + l0 ]= pressure[1][i]

                    param_real[ ptvx      + l + l0 ]= vx[1][i]
                    param_real[ ptvy      + l + l0 ]= vy[1][i]
                    param_real[ ptvz      + l + l0 ]= vz[1][i]

                    if utau is not None:
                        param_real[ ptutau    + l + l0 ]= utau[1][i]
                        param_real[ ptyplus   + l + l0 ]= yplus[1][i]

                    if temp_local is not None:
                        param_real[ pttemp_local       + l + l0 ]= temp_local[1][i]
                    if temp_extra_local is not None:
                        param_real[ pttemp_extra_local + l + l0 ]= temp_extra_local[1][i]
                    if temp_extra2_local is not None:
                        param_real[ pttemp_extra2_local + l + l0 ]= temp_extra2_local[1][i]

                    if wmodel_local_dens is not None:
                        param_real[ pt_dens_wm_local       + l + l0 ]= wmodel_local_dens[1][i]
                        param_real[ pt_velx_wm_local       + l + l0 ]= wmodel_local_velx[1][i]
                        param_real[ pt_vely_wm_local       + l + l0 ]= wmodel_local_vely[1][i]
                        param_real[ pt_velz_wm_local       + l + l0 ]= wmodel_local_velz[1][i]
                        param_real[ pt_temp_wm_local       + l + l0 ]= wmodel_local_temp[1][i]
                        param_real[ pt_sanu_wm_local       + l + l0 ]= wmodel_local_sanu[1][i]

                    if gradxP is not None:
                        param_real[ ptgradxP + l + l0 ]= gradxP[1][i]
                        param_real[ ptgradyP + l + l0 ]= gradyP[1][i]
                        param_real[ ptgradzP + l + l0 ]= gradzP[1][i]

                    if gradxU is not None:
                        param_real[ ptgradxU + l + l0 ]= gradxU[1][i]
                        param_real[ ptgradyU + l + l0 ]= gradyU[1][i]
                        param_real[ ptgradzU + l + l0 ]= gradzU[1][i]

                    if gradxV is not None:
                        param_real[ ptgradxV + l + l0 ]= gradxV[1][i]
                        param_real[ ptgradyV + l + l0 ]= gradyV[1][i]
                        param_real[ ptgradzV + l + l0 ]= gradzV[1][i]

                    if gradxW is not None:
                        param_real[ ptgradxW + l + l0 ]= gradxW[1][i]
                        param_real[ ptgradyW + l + l0 ]= gradyW[1][i]
                        param_real[ ptgradzW + l + l0 ]= gradzW[1][i]

                    if xcInit is not None:
                        param_real[ ptxcInit + l + l0 ]= xcInit[1][i]
                        param_real[ ptycInit + l + l0 ]= ycInit[1][i]
                        param_real[ ptzcInit + l + l0 ]= zcInit[1][i]

                        param_real[ ptxiInit + l + l0 ]= xiInit[1][i]
                        param_real[ ptyiInit + l + l0 ]= yiInit[1][i]
                        param_real[ ptziInit + l + l0 ]= ziInit[1][i]

                        param_real[ ptxwInit + l + l0 ]= xwInit[1][i]
                        param_real[ ptywInit + l + l0 ]= ywInit[1][i]
                        param_real[ ptzwInit + l + l0 ]= zwInit[1][i]

                    if motion_type is not None:
                        param_real[ ptmotion_type  + l + l0 ] = motion_type[1][i]

                        param_real[ pttransl_speedX + l + l0 ] = transl_speedX[1][i]
                        param_real[ pttransl_speedY + l + l0 ] = transl_speedY[1][i]
                        param_real[ pttransl_speedZ + l + l0 ] = transl_speedZ[1][i]

                        param_real[ ptaxis_pntX     + l + l0 ] = axis_pntX[1][i]
                        param_real[ ptaxis_pntY     + l + l0 ] = axis_pntY[1][i]
                        param_real[ ptaxis_pntZ     + l + l0 ] = axis_pntZ[1][i]

                        param_real[ ptaxis_vctX     + l + l0 ] = axis_vctX[1][i]
                        param_real[ ptaxis_vctY     + l + l0 ] = axis_vctY[1][i]
                        param_real[ ptaxis_vctZ     + l + l0 ] = axis_vctZ[1][i]

                        param_real[ ptomega        + l + l0 ] = omega[1][i]


                    if kcurv is not None:
                        param_real[ ptkcurv + l + l0 ]= kcurv[1][i]

                    if sd1 is not None:
                        param_real[ ptd1   + l + l0 ]= sd1[1][i]
                        param_real[ ptd2   + l + l0 ]= sd2[1][i]
                        param_real[ ptd3   + l + l0 ]= sd3[1][i]
                        param_real[ ptd4   + l + l0 ]= sd4[1][i]
                        param_real[ ptd5   + l + l0 ]= sd5[1][i]

                    if yline is not None:
                        param_real[ ptyline          + l + l0 ]= yline[1][i]
                        param_real[ ptuline          + l + l0 ]= uline[1][i]
                        param_real[ ptnutildeline    + l + l0 ]= nutildeline[1][i]
                        param_real[ ptpsiline        + l + l0 ]= psiline[1][i]
                        param_real[ ptmatmline       + l + l0 ]= matmline[1][i]
                        param_real[ ptmatline        + l + l0 ]= matline[1][i]
                        param_real[ ptmatpline       + l + l0 ]= matpline[1][i]
                        param_real[ ptalphasbetaline + l + l0 ]= alphasbetaline[1][i]
                        param_real[ ptindexline      + l + l0 ]= indexline[1][i]

                    if model == 'LBMLaminar':
                        if qloc_1[0] is not None:
                            for f_i in range (0, neq_loc):
                                param_real[ ptqloc_1[f_i]   + l + l0 ]= qloc_1[f_i][1][i]
                        if qloc_2[0] is not None:
                            for f_i in range (0, neq_loc):
                                param_real[ ptqloc_2[f_i]   + l + l0 ]= qloc_2[f_i][1][i]
                        if qloc_3[0] is not None:
                            for f_i in range (0, neq_loc):
                                param_real[ ptqloc_3[f_i]   + l + l0 ]= qloc_3[f_i][1][i]
                #recopie  InterpD
                param_real[ pt_coef + nocoef: pt_coef + nocoef+sizecoef] = InterpD[1][ nocoef_old: nocoef_old+sizecoef]
                nocoef     = nocoef     + sizecoef
                nocoef_old = nocoef_old + sizecoef
                l += 1

            else:
                if ntype != 0:
                    noi_old = noi_old + 1
                else:
                    ncfLoc  = pointlist[1][ noi_old ]
                    noi_old = noi_old + 1 + ncfLoc

                nocoef_old += sizecoef

        l0 = l0 + l

        param_int[ ptTy + ctyp +1 ] = l
        ctyp                        = ctyp +1

    return npass_transfert

#==============================================================================
# tri monotype
#==============================================================================
def triMonoType(Nbpts_D, Nbpts, Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                Interptype, pointlist, pointlistD, param_int,
                ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw,
                ptdensity,ptpressure,ptkcurv,
                ptvx, ptvy, ptvz,
                ptutau,ptyplus,
                pttemp_local, pttemp_extra_local, pttemp_extra2_local,
                pt_dens_wm_local,pt_velx_wm_local,pt_vely_wm_local,
                pt_velz_wm_local,pt_temp_wm_local,pt_sanu_wm_local,
                ptgradxP, ptgradyP, ptgradzP,
                ptgradxU, ptgradyU, ptgradzU,
                ptgradxV, ptgradyV, ptgradzV,
                ptgradxW, ptgradyW, ptgradzW,
                ptxcInit,ptycInit,ptzcInit,ptxiInit,ptyiInit,ptziInit,ptxwInit,ptywInit,ptzwInit,
                ptmotion_type,
                pttransl_speedX,pttransl_speedY,pttransl_speedZ,
                ptaxis_pntX,ptaxis_pntY,ptaxis_pntZ,
                ptaxis_vctX,ptaxis_vctY,ptaxis_vctZ,
                ptomega,
                ptd1,ptd2,ptd3,ptd4,ptd5,
                ptyline,ptuline,ptnutildeline,ptpsiline,ptmatmline,ptmatline,ptmatpline,
                ptalphasbetaline,ptindexline,
                xc,yc,zc,xi,yi,zi,xw,yw,zw,
                density,pressure,kcurv,
                vx, vy, vz,
                utau,yplus,
                temp_local, temp_extra_local, temp_extra2_local,
                wmodel_local_dens,wmodel_local_velx,wmodel_local_vely,
                wmodel_local_velz,wmodel_local_temp,wmodel_local_sanu,
                gradxP, gradyP, gradzP,
                gradxU, gradyU, gradzU,
                gradxV, gradyV, gradzV,
                gradxW, gradyW, gradzW,
                xcInit,ycInit,zcInit,xiInit,yiInit,ziInit,xwInit,ywInit,zwInit,
                motion_type,
                transl_speedX,transl_speedY,transl_speedZ,
                axis_pntX,axis_pntY,axis_pntZ,
                axis_vctX,axis_vctY,axis_vctZ,
                omega,
                sd1,sd2,sd3,sd4,sd5,
                yline,uline,nutildeline,psiline,matmline,matline,matpline,
                alphasbetaline,indexline,
                InterpD, param_real,ptqloc_1,qloc_1,ptqloc_2,qloc_2,ptqloc_3,qloc_3,neq_loc,model,nbpts_linelets):

    ntype     = Nbtype[0]
    noi_old   = 0
    nocoef_old= 0
    l         = 0
    ltype     = Interptype[1][0]

    npass_transfert = 1
    if ltype ==44 or ltype == 5: npass_transfert = 2

    #recopieinterpolantType
    ideb =  ptTy + shift_typ
    val = float(ltype)
    connector.initNuma(None, param_int, ideb, Nbpts_D, 1, val)

    # recopie pointlist
    ideb = lst
    connector.initNuma(pointlist[1], param_int, ideb, Nbpts, 1, val)
    #recopie pointListDonor
    ideb = lstD
    connector.initNuma(pointlistD[1], param_int, ideb, Nbpts_D, 1, val)
    #recopie Maillage IBC
    if sname == 'IB':
        connector.initNuma(xc[1], param_real, ptxc, Nbpts_D , 0, val)
        connector.initNuma(yc[1], param_real, ptyc, Nbpts_D , 0, val)
        connector.initNuma(zc[1], param_real, ptzc, Nbpts_D , 0, val)
        connector.initNuma(xi[1], param_real, ptxi, Nbpts_D , 0, val)
        connector.initNuma(yi[1], param_real, ptyi, Nbpts_D , 0, val)
        connector.initNuma(zi[1], param_real, ptzi, Nbpts_D , 0, val)
        connector.initNuma(xw[1], param_real, ptxw, Nbpts_D , 0, val)
        connector.initNuma(yw[1], param_real, ptyw, Nbpts_D , 0, val)
        connector.initNuma(zw[1], param_real, ptzw, Nbpts_D , 0, val)

        connector.initNuma(density[1], param_real, ptdensity , Nbpts_D , 0, val)
        connector.initNuma(pressure[1], param_real, ptpressure, Nbpts_D , 0, val)

        connector.initNuma(vx[1], param_real, ptvx, Nbpts_D , 0, val)
        connector.initNuma(vy[1], param_real, ptvy, Nbpts_D , 0, val)
        connector.initNuma(vz[1], param_real, ptvz, Nbpts_D , 0, val)

        if utau is not None:
            connector.initNuma(utau[1], param_real, ptutau, Nbpts_D, 0, val)
            connector.initNuma(yplus[1], param_real, ptyplus, Nbpts_D , 0, val)

        if temp_local is not None:
            connector.initNuma(temp_local[1]      , param_real, pttemp_local       , Nbpts_D , 0, val)
        if temp_extra_local is not None:
            connector.initNuma(temp_extra_local[1], param_real, pttemp_extra_local , Nbpts_D , 0, val)
        if temp_extra2_local is not None:
            connector.initNuma(temp_extra2_local[1], param_real, pttemp_extra2_local , Nbpts_D , 0, val)

        if wmodel_local_dens is not None:
            connector.initNuma(wmodel_local_dens[1]      , param_real, pt_dens_wm_local       , Nbpts_D , 0, val)
            connector.initNuma(wmodel_local_velx[1]      , param_real, pt_velx_wm_local       , Nbpts_D , 0, val)
            connector.initNuma(wmodel_local_vely[1]      , param_real, pt_vely_wm_local       , Nbpts_D , 0, val)
            connector.initNuma(wmodel_local_velz[1]      , param_real, pt_velz_wm_local       , Nbpts_D , 0, val)
            connector.initNuma(wmodel_local_temp[1]      , param_real, pt_temp_wm_local       , Nbpts_D , 0, val)
            connector.initNuma(wmodel_local_sanu[1]      , param_real, pt_sanu_wm_local       , Nbpts_D , 0, val)


        if gradxP is not None:
            connector.initNuma(gradxP[1] , param_real, ptgradxP , Nbpts_D , 0, val)
            connector.initNuma(gradyP[1] , param_real, ptgradyP , Nbpts_D , 0, val)
            connector.initNuma(gradzP[1] , param_real, ptgradzP , Nbpts_D , 0, val)

        if gradxU is not None:
            connector.initNuma(gradxU[1] , param_real, ptgradxU , Nbpts_D , 0, val)
            connector.initNuma(gradyU[1] , param_real, ptgradyU , Nbpts_D , 0, val)
            connector.initNuma(gradzU[1] , param_real, ptgradzU , Nbpts_D , 0, val)

        if gradxV is not None:
            connector.initNuma(gradxV[1] , param_real, ptgradxV , Nbpts_D , 0, val)
            connector.initNuma(gradyV[1] , param_real, ptgradyV , Nbpts_D , 0, val)
            connector.initNuma(gradzV[1] , param_real, ptgradzV , Nbpts_D , 0, val)

        if gradxW is not None:
            connector.initNuma(gradxW[1] , param_real, ptgradxW , Nbpts_D , 0, val)
            connector.initNuma(gradyW[1] , param_real, ptgradyW , Nbpts_D , 0, val)
            connector.initNuma(gradzW[1] , param_real, ptgradzW , Nbpts_D , 0, val)

        if xcInit is not None:
            connector.initNuma(xcInit[1] , param_real, ptxcInit , Nbpts_D , 0, val)
            connector.initNuma(ycInit[1] , param_real, ptycInit , Nbpts_D , 0, val)
            connector.initNuma(zcInit[1] , param_real, ptzcInit , Nbpts_D , 0, val)

            connector.initNuma(xiInit[1] , param_real, ptxiInit , Nbpts_D , 0, val)
            connector.initNuma(yiInit[1] , param_real, ptyiInit , Nbpts_D , 0, val)
            connector.initNuma(ziInit[1] , param_real, ptziInit , Nbpts_D , 0, val)

            connector.initNuma(xwInit[1] , param_real, ptxwInit , Nbpts_D , 0, val)
            connector.initNuma(ywInit[1] , param_real, ptywInit , Nbpts_D , 0, val)
            connector.initNuma(zwInit[1] , param_real, ptzwInit , Nbpts_D , 0, val)

        if motion_type is not None:
            connector.initNuma(motion_type[1]  , param_real, ptmotion_type  , Nbpts_D , 0, val)

            connector.initNuma(transl_speedX[1] , param_real, pttransl_speedX , Nbpts_D , 0, val)
            connector.initNuma(transl_speedY[1] , param_real, pttransl_speedY , Nbpts_D , 0, val)
            connector.initNuma(transl_speedZ[1] , param_real, pttransl_speedZ , Nbpts_D , 0, val)

            connector.initNuma(axis_pntX[1]     , param_real, ptaxis_pntX     , Nbpts_D , 0, val)
            connector.initNuma(axis_pntY[1]     , param_real, ptaxis_pntY     , Nbpts_D , 0, val)
            connector.initNuma(axis_pntZ[1]     , param_real, ptaxis_pntZ     , Nbpts_D , 0, val)

            connector.initNuma(axis_vctX[1]     , param_real, ptaxis_vctX     , Nbpts_D , 0, val)
            connector.initNuma(axis_vctY[1]     , param_real, ptaxis_vctY     , Nbpts_D , 0, val)
            connector.initNuma(axis_vctZ[1]     , param_real, ptaxis_vctZ     , Nbpts_D , 0, val)

            connector.initNuma(omega[1]        , param_real, ptomega        , Nbpts_D , 0, val)

        if kcurv is not None:
            connector.initNuma(kcurv[1] , param_real, ptkcurv , Nbpts_D , 0, val)

        if sd1 is not None:
            connector.initNuma(sd1[1], param_real, ptd1 , Nbpts_D , 0, val)
        if sd2 is not None:
            connector.initNuma(sd2[1], param_real, ptd2 , Nbpts_D , 0, val)
        if sd3 is not None:
            connector.initNuma(sd3[1] , param_real, ptd3 , Nbpts_D , 0, val)
        if sd4 is not None:
            connector.initNuma(sd4[1] , param_real, ptd4 , Nbpts_D , 0, val)
        if sd5 is not None:
            connector.initNuma(sd5[1] , param_real, ptd5 , Nbpts_D , 0, val)

        if yline is not None:
            connector.initNuma(yline[1]          , param_real, ptyline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(uline[1]          , param_real, ptuline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(nutildeline[1]    , param_real, ptnutildeline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(psiline[1]        , param_real, ptpsiline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(matmline[1]       , param_real, ptmatmline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(matline[1]        , param_real, ptmatline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(matpline[1]       , param_real, ptmatpline, Nbpts_D*nbpts_linelets, 0, val)
            connector.initNuma(alphasbetaline[1] , param_real, ptalphasbetaline, Nbpts_D, 0, val)
            connector.initNuma(indexline[1]      , param_real, ptindexline, Nbpts_D, 0, val)

        if model=='LBMLaminar':
            if qloc_1[0] is not None:
                for f_i in range (0,neq_loc):
                    connector.initNuma(qloc_1[f_i][1] , param_real, ptqloc_1[f_i] , Nbpts_D , 0, val)
            if qloc_2[0] is not None:
                for f_i in range (0,neq_loc):
                    connector.initNuma(qloc_2[f_i][1] , param_real, ptqloc_2[f_i] , Nbpts_D , 0, val)
            if qloc_3[0] is not None:
                for f_i in range (0,neq_loc):
                    connector.initNuma(qloc_3[f_i][1] , param_real, ptqloc_3[f_i] , Nbpts_D , 0, val)
    # recopie  InterpD
    connector.initNuma(InterpD[1] , param_real, pt_coef , Nbpts_InterpD , 0, val)

    param_int[ ptTy + ctyp +1 ] = Nbpts_D

    return npass_transfert

#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la zone donneuse
# fonctionne avec ___setInterpTransfer
#==============================================================================
def miseAPlatDonorZone__(zones, tc, procDict):
    zones_tc = Internal.getZones(tc)
    #[AJ]
    base     = Internal.getNodeFromType1(tc, 'CGNSBase_t')  # noeud
    model    = 'NSLaminar'
    a        = Internal.getNodeFromName2(base, 'GoverningEquations')
    if a is not None: model = Internal.getValue(a)

    neq_loc = 5
    if model=='NSTurbulent':
        neq_loc = 6
    elif model=='LBMLaminar':
        neq_loc = Internal.getNodeFromName2(zones[0] , 'Parameter_int')[1][NEQ_LBM]

    for z in zones_tc:
        racs      =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        size_int  = 0
        size_real = 0
        count_ID  = 0
        count_IBC = 0
        # alloc memoire
        for rac in racs:
            pointlist    =  Internal.getNodeFromName1(rac, 'PointList')
            pointlistD   =  Internal.getNodeFromName1(rac, 'PointListDonor')
            InterpD      =  Internal.getNodeFromName1(rac, 'InterpolantsDonor')
            utau         =  Internal.getNodeFromName1(rac, 'utau')
            temp_local   =  Internal.getNodeFromName1(rac, 'Temperature')
            wmodel_local =  Internal.getNodeFromName1(rac, 'Density_WM')
            gradxP       =  Internal.getNodeFromName1(rac, 'gradxPressure')
            gradxU       =  Internal.getNodeFromName1(rac, 'gradxVelocityX')
            xcInit       =  Internal.getNodeFromName1(rac, 'CoordinateX_PC#Init')
            motion_type  =  Internal.getNodeFromName1(rac, 'MotionType')
            kcurv        =  Internal.getNodeFromName1(rac, 'KCurv')
            sd1          =  Internal.getNodeFromName1(rac, 'StagnationEnthalpy')
            yline        =  Internal.getNodeFromName1(rac, 'CoordinateN_ODE')
            qloc_1       =  Internal.getNodeFromName1(rac, ibm_lbm_variables_1 + str(1))

            ntab_IBC   = 11+3 #On ajoute dorenavant les vitesses dans l'arbre tc pour le post
            if utau is not None: ntab_IBC += 2
            if temp_local is not None: ntab_IBC += 2
            if wmodel_local is not None: ntab_IBC += 6
            if gradxP is not None: ntab_IBC += 3
            if gradxU is not None: ntab_IBC += 9
            if xcInit is not None: ntab_IBC += 9
            if motion_type is not None: ntab_IBC += 11 #MotionType,transl_speed(x3),axis_pnt(x3),axis_vct(x3),omega
            if kcurv is not None: ntab_IBC += 1
            if sd1 is not None: ntab_IBC += 5
            if yline is not None: ntab_IBC += (7*nbpts_linelets+2)
            if qloc_1 is not None: ntab_IBC += neq_loc

            Nbpts        =  numpy.shape(pointlist[ 1])[0]
            Nbpts_D      =  numpy.shape(pointlistD[1])[0]
            Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

            size_int   =  size_int + 7 + Nbpts_D*2 + Nbpts
            size_real  =  size_real+   Nbpts_InterpD
            sname = rac[0][0:2]
            if sname == 'IB':
                size_real = size_real +Nbpts_D*ntab_IBC
                count_IBC = count_IBC +1
            else:
                count_ID  = count_ID  +1
            #print('nbpt, nbpt_donor', sname,Nbpts,Nbpts_InterpD)

        size_int = size_int + 2 + (count_IBC + count_ID)*2  # 2: nbr rac ID et IBC, stockage adresse debut raccord et coef
        param_int  = numpy.empty(size_int , dtype=Internal.E_NpyInt)
        param_real = numpy.empty(size_real, dtype=numpy.float64)
        Internal.createUniqueChild(z, 'Parameter_int' , 'DataArray_t', param_int )
        if size_real !=0 :
            Internal.createUniqueChild(z, 'Parameter_real', 'DataArray_t', param_real)# recopie pointlis

        #print('size int et real', size_int, size_real)
        #print('ID et IBC', count_ID, count_IBC)

        param_int[0] = count_ID
        param_int[1] = count_IBC
        #initialisation
        c        = 0
        size_rac = 0
        size_coef= 0
        for rac in racs:
            pt_rac = 2 + (count_ID + count_IBC)*2 + size_rac # adresse debut raccord
            pt_coef= size_coef                               # adresse debut coef
            #print 'indice', pt_rac , Nbpts,count_ID,count_IBC, size_rac

            pointlist    =  Internal.getNodeFromName1(rac, 'PointList')
            pointlistD   =  Internal.getNodeFromName1(rac, 'PointListDonor')
            Interptype   =  Internal.getNodeFromName1(rac, 'InterpolantsType')
            InterpD      =  Internal.getNodeFromName1(rac, 'InterpolantsDonor')
            Nbpts        =  numpy.shape(pointlist[ 1])[0]
            Nbpts_D      =  numpy.shape(pointlistD[1])[0]
            Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

            param_int[ 2+c                    ] = pt_rac
            param_int[ 2+c+count_ID+count_IBC ] = pt_coef
            param_int[ pt_rac                 ] = Nbpts
            param_int[ pt_rac +1              ] = Nbpts_D
            param_int[ pt_rac +2              ] = Nbpts_InterpD

            typecell = Interptype[1][0]
            #on cree un tableau qui contient les elements non egaux au premier element
            b = Interptype[1] [ Interptype[1] != typecell ]
            if len(b) == 0:   param_int[ pt_rac +3 ] = 0    # type homogene
            else:             param_int[ pt_rac +3 ] = 1    # type melange

            tmp =  Internal.getNodeFromName1(rac, 'ZoneRole')
            if tmp[1][0] == b'D': param_int[ pt_rac +4 ] = 0           # role= Donor
            else                : param_int[ pt_rac +4 ] = 1           # role= Receiver
            tmp =  Internal.getNodeFromName1(rac, 'GridLocation')
            if tmp[1][4] == b'C': param_int[ pt_rac +5 ] = 1           # location= CellCenter
            else                : param_int[ pt_rac +5 ] = 0           # location= Vertex

            zrcvname = Internal.getValue(rac)
            no_zone = 0
            for z0 in zones:
                if z0[0] == zrcvname: param_int[ pt_rac +6 ]= no_zone  # No zone raccord
                no_zone += 1

            ideb =  pt_rac +7
            param_int[ ideb:ideb + Nbpts   ] = pointlist[ 1][0:Nbpts  ]           # recopie pointlist
            pointlist[ 1]                    = param_int[ ideb : ideb + Nbpts ]   # supression numpy initial

            ideb =  ideb  + Nbpts
            param_int[ ideb:ideb+ Nbpts_D  ] = pointlistD[1][0:Nbpts_D]           # recopie pointlistdonor
            pointlistD[ 1]                   = param_int[ ideb : ideb + Nbpts_D ] # supression numpy initial

            ideb =  ideb  + Nbpts_D
            param_int[ ideb:ideb + Nbpts_D ] = Interptype[1][0:Nbpts_D]           #recopieinterpolantType
            Interptype[ 1]                   = param_int[ ideb : ideb + Nbpts_D ] # supression numpy initial

            size_rac   =  size_rac + 7 + Nbpts_D*2 + Nbpts

            param_real[ pt_coef:pt_coef + Nbpts_InterpD   ] = InterpD[1][0:Nbpts_InterpD]
            ### supression numpy initial
            InterpD[1] = param_real[ pt_coef:pt_coef + Nbpts_InterpD ]

            size_coef = size_coef + Nbpts_InterpD

            sname = rac[0][0:2]
            if sname == 'IB':
                if utau is not None:
                    var_ibc=['CoordinateX_PC','CoordinateY_PC','CoordinateZ_PC','CoordinateX_PI','CoordinateY_PI','CoordinateZ_PI','CoordinateX_PW','CoordinateY_PW','CoordinateZ_PW', 'Density','Pressure','VelocityX','VelocityY','VelocityZ','utau','yplus']
                    if gradxP is not None:
                        var_ibc.append('gradxPressure')
                        var_ibc.append('gradyPressure')
                        var_ibc.append('gradzPressure')
                    if gradxU is not None:
                        var_ibc.append('gradxVelocityX')
                        var_ibc.append('gradyVelocityX')
                        var_ibc.append('gradzVelocityX')
                    if gradxV is not None:
                        var_ibc.append('gradxVelocityY')
                        var_ibc.append('gradyVelocityY')
                        var_ibc.append('gradzVelocityY')
                    if gradxW is not None:
                        var_ibc.append('gradxVelocityZ')
                        var_ibc.append('gradyVelocityZ')
                        var_ibc.append('gradzVelocityZ')
                    if gradxW is not None:
                        var_ibc.append('gradxVelocityZ')
                        var_ibc.append('gradyVelocityZ')
                        var_ibc.append('gradzVelocityZ')
                    if xcInit is not None:
                        var_ibc.append('CoordinateX_PC#Init')
                        var_ibc.append('CoordinateY_PC#Init')
                        var_ibc.append('CoordinateZ_PC#Init')

                        var_ibc.append('CoordinateX_PI#Init')
                        var_ibc.append('CoordinateY_PI#Init')
                        var_ibc.append('CoordinateZ_PI#Init')

                        var_ibc.append('CoordinateX_PW#Init')
                        var_ibc.append('CoordinateY_PW#Init')
                        var_ibc.append('CoordinateZ_PW#Init')
                    if motion_type is not None:
                        var_ibc.append('MotionType')

                        var_ibc.append('transl_speedX')
                        var_ibc.append('transl_speedY')
                        var_ibc.append('transl_speedZ')

                        var_ibc.append('axis_pntX')
                        var_ibc.append('axis_pntY')
                        var_ibc.append('axis_pntZ')

                        var_ibc.append('axis_vctX')
                        var_ibc.append('axis_vctY')
                        var_ibc.append('axis_vctZ')

                        var_ibc.append('omega')
                    if kcurv is not None:
                        var_ibc.append('KCurv')
                    if yline is not None:
                        var_ibc.append('CoordinateN_ODE')
                        var_ibc.append('VelocityT_ODE')
                        var_ibc.append('TurbulentSANuTilde_ODE')
                        var_ibc.append('Psi_ODE')
                        var_ibc.append('Matm_ODE')
                        var_ibc.append('Mat_ODE')
                        var_ibc.append('Matp_ODE')
                        var_ibc.append('alphasbeta_ODE')
                        var_ibc.append('index_ODE')
                else:
                    var_ibc=['CoordinateX_PC','CoordinateY_PC','CoordinateZ_PC','CoordinateX_PI','CoordinateY_PI','CoordinateZ_PI','CoordinateX_PW','CoordinateY_PW','CoordinateZ_PW', 'Density','Pressure','VelocityX','VelocityY','VelocityZ']

                    if xcInit is not None:
                        var_ibc.append('CoordinateX_PC#Init')
                        var_ibc.append('CoordinateY_PC#Init')
                        var_ibc.append('CoordinateZ_PC#Init')

                        var_ibc.append('CoordinateX_PI#Init')
                        var_ibc.append('CoordinateY_PI#Init')
                        var_ibc.append('CoordinateZ_PI#Init')

                        var_ibc.append('CoordinateX_PW#Init')
                        var_ibc.append('CoordinateY_PW#Init')
                        var_ibc.append('CoordinateZ_PW#Init')

                    if motion_type is not None:
                        var_ibc.append('MotionType')

                        var_ibc.append('transl_speedX')
                        var_ibc.append('transl_speedY')
                        var_ibc.append('transl_speedZ')

                        var_ibc.append('axis_pntX')
                        var_ibc.append('axis_pntY')
                        var_ibc.append('axis_pntZ')

                        var_ibc.append('axis_vctX')
                        var_ibc.append('axis_vctY')
                        var_ibc.append('axis_vctZ')

                        var_ibc.append('omega')

                count_ibc = 0
                ideb      = pt_coef + Nbpts_InterpD
                for v_ibc in var_ibc:
                    tmp                            = Internal.getNodeFromName1(rac, v_ibc)
                    param_real[ ideb:ideb+ Nbpts_D]= tmp[1][0:Nbpts_D]
                    tmp[1]                         = param_real[ ideb:ideb+ Nbpts_D ]
                    ideb                           = ideb + Nbpts_D
                    count_ibc += 1

                size_coef = size_coef + count_ibc*Nbpts_D

            c = c+1

#===============================================================================
# General transfers: Chimera + IBC - inplace version optimiser par arbre tc compacte par zone donneuse
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to be
# consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                = None: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# IN: bcType (IBC only) 0: glissement
#                       1: adherence
#                       2: loi de paroi log
#                       3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE)+ronultideSA
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t)+nultideSA
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p)+nultideSA
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def ___setInterpTransfers(aR, topTreeD,
                          variables=[],
                          variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                          bcType=0, varType=1, storage=-1, compact=0,
                          Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                          Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage != 1 or compact == 0:
        raise ValueError("___setInterpTransfers: Mode receveur a coder. Mode compact obligatoire: set compact=1 in Fast.warmup.")

    # Recup des donnees a partir des zones donneuses
    zones     = Internal.getZones(aR)
    zonesD    = Internal.getZones(topTreeD)
    param_int = Internal.getNodeFromName2(topTreeD, 'Parameter_int' )[1]
    param_real= Internal.getNodeFromName2(topTreeD, 'Parameter_real')[1]
    connector.___setInterpTransfers(zones, zonesD, variables, param_int, param_real, varType, bcType, Gamma, Cv, MuS, Cs, Ts )
    return None

#===============================================================================
# calcul les processus source pour rank ( graphrcv_) et la position dans param_int (pos_list)
# filterGraph:: pour les cas unsteady, permet de ne pas ajouter les sources fournies par raccord steady
#
def _procSource(rank, S_pos, pos_list, graph, graphloc, graphrcv_, filterGraph=None):

    pos_list.append(S_pos)
    k_pos= 0
    for proc in graph.keys():
        for n in graph[proc].keys():
            if n == rank:
                if filterGraph is None:
                    graphloc.append(proc)
                    k_pos  +=1
                elif proc not in filterGraph:
                    graphloc.append(proc)
                    k_pos  +=1

    graphrcv_.append(k_pos)
    S_pos += k_pos +1
    for proc in graphloc: graphrcv_.append(proc)
    #print("Spos=", S_pos,graphrcv_, filterGraph )

    return S_pos
