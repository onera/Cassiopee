# - gestion du compactage et des transferts compacts - 
from . import connector
import numpy

try: import Converter.Internal as Internal
except: raise ImportError("Connector.compactTransfers requires Converter module.")

try: range = xrange
except: pass

#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la base
# fonctionne avec ___setInterpTransfer
#==============================================================================
def miseAPlatDonorTree__(zones, tc, graph=None, list_graph=None):
    try: import Converter.GhostCells as GhostCells
    except: raise ImportError("setInterpDataForGhostCells__ requires Converter.GhostCells module.")

    if graph is not None :
        # Suppose mpirun
        procDict  = graph['procDict']
        procList  = graph['procList']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
        import Converter.Mpi as Cmpi
        rank = Cmpi.rank
    else: 
        procDict=None; procList=None; graphID=None; graphIBCD=None; rank=0

    # if Cmpi is not None and rank == 0: 
    #    print "GRAPH IBC IS : ",graph['graphIBCD']
    #    print "graphID IS :",graph['graphID']
    
    # print "procDict IS : ",procDict    
    size_int  = 0
    size_real = 0
    listproc  = []
    rac       = []
    rac_inst  = []
    sizeI     = []
    sizeR     = []
    sizeNbD   = []
    sizeType  = []
    nrac      = 0

    ordered_subRegions =[]
    No_zoneD =[]
    MeshTypeD=[]
    inst = {}
    numero_max =-100000000
    numero_min = 100000000

    zones_tc = Internal.getZones(tc)
    c        = 0
    for z in zones_tc:
      #print z[0]
      subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
      meshtype   = 1
      zonetype   = Internal.getNodeFromType1(z, 'ZoneType_t')
      tmp        = Internal.getValue(zonetype)
      dimPb      = Internal.getZoneDim(z)[4]
      if tmp != "Structured": meshtype = 2
      for s in subRegions:
         zRname = Internal.getValue(s)
         proc = 0
         if procDict is not None: proc = procDict[zRname]

         #tri des pas de temps instationnaire
         #  1) les stationnaires
         #  2) les instationnaires regroupes par pas de temps
         if '#' not in s[0]:
              ordered_subRegions.append(s)
              No_zoneD.append(c)
              MeshTypeD.append(meshtype)
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
                inst[ numero_iter ]=  [ sub , Noz , mesh, dest ]
            else:
                inst[ numero_iter ]= [ [s],[c],[meshtype], [proc] ]

         TimeLevelNumber = len(inst)

         if TimeLevelNumber != 1+numero_max-numero_min and len(inst)!= 0: 
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
         prange        =  Internal.getNodeFromName1(s, 'PointRange')       # Besoin des point range pour l'explicite local
         pranged       =  Internal.getNodeFromName1(s, 'PointRangeDonor')  # Besoin des point range pour l'explicite local 
         direction     =  Internal.getNodeFromName1(s, 'DirReceveur')       # Besoin des directions pour l'explicite local
         directiond    =  Internal.getNodeFromName1(s, 'DirDonneur')  # Besoin des point directions pour l'explicite local 
         transfo       =  Internal.getNodeFromName1(s, 'Transform')  # Besoin du transform pour l'explicite local
         pt_pivot      =  Internal.getNodeFromName1(s, 'PointPivot')  # Besoin du point pivot pour l'explicite local (conservativite)
         profondeur    =  Internal.getNodeFromName1(s, 'Profondeur')  # Besoin de la profondeur pour l'explicite local (nearmatch) 
         ratio         =  Internal.getNodeFromName1(s, 'NMratio') # Besoin des ratios entre les pas d espace des zones donneuse et receveuse (exp local)
         levelrcv      =  Internal.getNodeFromName1(s, 'LevelZRcv') # Niveau en temps zone receveuse (exp local)
         leveldnr      =  Internal.getNodeFromName1(s, 'LevelZDnr') # Niveau en temps zone donneuse (exp local)
         #print 'prange= ',prange

         Nbpts        =  numpy.shape(pointlist[ 1])[0]
         Nbpts_D      =  numpy.shape(pointlistD[1])[0]
         Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

         sname = s[0][0:2]
         utau = Internal.getNodeFromName1(s, 'utau')
         sd1 = Internal.getNodeFromName1(s, 'StagnationEnthalpy')
         
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
         #print 'nb type',  len(Nbtype), s[0],z[0], Nbtype
         nbType = numpy.unique(Interptype[1])
         nbTypeSize = nbType.size
         #print len(Nbtype), Nbtype2.size

         size_IBC =  0
         ntab_IBC = 11+3 #On ajoute dorenavant les vitesses dans l'arbre tc pour faciliter le post
         if utau is not None: ntab_IBC += 2
         if sd1 is not None: ntab_IBC += 5
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
                sizeI.append(    Nbpts_D*2     + Nbpts   + nbTypeSize+1 )
                sizeR.append(    Nbpts_InterpD + size_IBC + rotation    )
                sizeNbD.append(  Nbpts_D                                )
                sizeType.append( Nbpts_D                 + nbTypeSize+1 )
         else:
                pos           = listproc.index(proc)
                rac[pos]      = rac[pos] + 1
                if '#' in s[0]: rac_inst[pos]= rac_inst[pos] + 1
                sizeI[pos]    = sizeI[pos]     + Nbpts_D*2     + Nbpts   + nbTypeSize+1
                sizeR[pos]    = sizeR[pos]     + Nbpts_InterpD + size_IBC + rotation
                sizeNbD[pos]  = sizeNbD[pos]   + Nbpts_D
                sizeType[pos] = sizeType[pos]  + Nbpts_D                 + nbTypeSize+1

      c += 1

    base     = Internal.getNodeFromType1(tc, 'CGNSBase_t')  # noeud
    model    = 'NSLaminar'
    a        = Internal.getNodeFromName2(base, 'GoverningEquations')
    if a is not None: model = Internal.getValue(a)

    NbP2P     = len(listproc)
    sizeproc  = []
    ntab_int  = 16

    if nrac != 0:
        if prange is not None:
             ntab_int = ntab_int + 27
            

    for i in range(NbP2P): sizeproc.append(5 + TimeLevelNumber*2 + ntab_int*rac[i] + sizeI[i])
           
    size_int  =  1 + NbP2P + sum(sizeproc)
    size_real =  sum(sizeR)

    ### list_graph est necessaire pour l'explicite local ####
    
    if list_graph is None: # Si list_graph n existe pas, on est pas en explicite local
        graphIBCrcv=[]
        if graphIBCD is not None:
            for proc in graphIBCD.keys():
                for n in graphIBCD[proc].keys():
                    if n == rank: graphIBCrcv.append(proc)

        graphIDrcv=[]
        if graphID is not None:
            for proc in graphID.keys():
                for n in graphID[proc].keys():
                    if n == rank: graphIDrcv.append(proc)
    
    else:  # Si list_graph existe, on est en explicite local, 1 graphe par ss ite      
        graphIBCrcv_=[]
        graphIDrcv_=[]
        pos_ID=[]
        pos_IBC=[]
        S_IBC=0
        S_ID=0

        for nstep in range(0,len(list_graph)):

            graphID_   = list_graph[nstep]['graphID']
            graphIBCD_ = list_graph[nstep]['graphIBCD']

            pos_IBC.append(len(list_graph)+S_IBC)
            k_IBC=0
            stokproc=[]
            for proc in graphIBCD_.keys():
               for n in graphIBCD_[proc].keys():
                     if n == rank:
                         k_IBC+=1
                         stokproc.append(proc)
            graphIBCrcv_.append(k_IBC)
            S_IBC += k_IBC+1
            for proc in stokproc:
                graphIBCrcv_.append(proc)              


            pos_ID.append(len(list_graph)+S_ID)            
            k_ID=0
            stokproc=[]
            for proc in graphID_.keys():
                for n in graphID_[proc].keys():
                     if n == rank:
                         k_ID+=1
                         stokproc.append(proc)
            graphIDrcv_.append(k_ID)
            S_ID += k_ID+1
            for proc in stokproc:
                graphIDrcv_.append(proc)

        graphIBCrcv  = pos_IBC + graphIBCrcv_    
        graphIDrcv   = pos_ID  + graphIDrcv_    

            #if sum(graphIBCrcv) < 1 : graphIBCrcv=[]
        
    
    #graphIDrcv=[0,3]
    #print "len graphIBCrcv is",len(graphIBCrcv)
    #print "len graphIDrcv  is",len(graphIDrcv)
    #print "pos_IBC is",pos_IBC
    #print "pos_ID  is",pos_ID

    param_int  = numpy.empty(size_int + len(graphIDrcv) + len(graphIBCrcv) + 2, dtype=numpy.int32  )
    param_real = numpy.empty(size_real, dtype=numpy.float64)
    Internal.createUniqueChild(tc, 'Parameter_int' , 'DataArray_t', param_int)
    if size_real !=0 : 
        Internal.createUniqueChild(tc, 'Parameter_real', 'DataArray_t', param_real)

    if len(graphIBCrcv) == 0:
         _graphIBC = numpy.zeros(1,dtype=numpy.int32)
    else:
        _graphIBC  = numpy.asarray([len(graphIBCrcv)]+graphIBCrcv,dtype=numpy.int32)  

    _graphID   = numpy.asarray([len(graphIDrcv)] +graphIDrcv ,dtype=numpy.int32)

    param_int[1                 :1+len(graphIBCrcv)+1] = _graphIBC
    param_int[2+len(graphIBCrcv):2+len(graphIBCrcv)+len(graphIDrcv)+1] = _graphID    

    # print "param_int is ",param_int[0:2+len(graphIBCrcv)+len(graphIDrcv)+1]
    # Dictionnaire pour optimisation
    znd = []
    for z in zones: znd.append(z[0])

    #
    #initialisation numpy
    #
    param_int[0] = NbP2P
    size_ptlist = []
    size_ptlistD= []
    size_ptType = []
    nb_rac      = []
    size_coef   = []
    adr_coef    = []   # pour cibler debut de echange dans param_real

    shift_graph = len(graphIDrcv) + len(graphIBCrcv) + 2 
    # print "shift_graph is ",shift_graph
    shift_coef  =0
    shift       = shift_graph # le shift prend en compte la postion des graphs (ID+IBC) entre la address contenant NbP2P et 
    for i in range(NbP2P):
       adr_coef.append(shift_coef)                    #adresse echange dans param_real
       shift_coef = shift_coef + sizeR[i]

       #param_int[i+1] = 1 + NbP2P + shift
       param_int[i+shift_graph+1] = 1 + NbP2P + shift              #adresse echange
       shift          =  shift  + sizeproc[i]
       size_ptlist.append(0)
       size_ptlistD.append(0)
       size_ptType.append(0)
       size_coef.append(0)
       nb_rac.append(0)

    for iter in range(numero_min,numero_max+1): 
        ordered_subRegions =  ordered_subRegions + inst[ iter ][0]
        No_zoneD           =  No_zoneD           + inst[ iter ][1]
        MeshTypeD          =  MeshTypeD          + inst[ iter ][2]

    # loop sur les raccords tries
    c        = 0
    Nbtot    = 0
    S = 0
    for s in ordered_subRegions:

       NozoneD  = No_zoneD[c]
       meshtype = MeshTypeD[c]

       zRname = Internal.getValue(s)
       proc = 0
       if procDict is not None: proc = procDict[zRname]
       pos  = listproc.index(proc) 
       #print 'pos= ',pos
       pt_ech = param_int[ pos + 1 + shift_graph]                 # adresse debut raccord pour l'echange pos
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

       #print zRname, nb_rac[pos]

       Nbpts         =  numpy.shape(pointlist[ 1])[0]
       Nbpts_D       =  numpy.shape(pointlistD[1])[0]
       Nbpts_InterpD =  numpy.shape(InterpD[ 1  ])[0]

       param_int[ pt_ech    ] = proc
       param_int[ pt_ech +1 ] = rac[pos]
       param_int[ pt_ech +2 ] = rac_inst[pos]
       nrac_steady            = rac[pos] - rac_inst[pos]

       param_int[ pt_ech +3 ] = TimeLevelNumber
       nrac_inst_deb  =  nrac_steady
       for i in range(TimeLevelNumber):
            # len(inst[i][0])  = nb de raccord instationnaire pour le temps i
            nrac_inst_fin  = nrac_inst_deb + len(inst[i+numero_min][0])

            param_int[ pt_ech +4 + i                  ] = nrac_inst_deb
            param_int[ pt_ech +4 + i + TimeLevelNumber] = nrac_inst_fin

            nrac_inst_deb  = nrac_inst_fin
 
       iadr = pt_ech + 4 + TimeLevelNumber*2 + nb_rac[pos]   # ptr echange + dest + nrac + norac
       iadr2= pt_ech + 4 + TimeLevelNumber*2 + 1

       param_int[ iadr            ] = Nbpts
       param_int[ iadr + rac[pos] ] = Nbpts_InterpD

       #on recupere le nombre de type different
       typecell = Interptype[1][0]
       Nbtype= [ typecell ]
       for i in range(Nbpts_D):
          if Interptype[1][i] not in Nbtype: Nbtype.append(Interptype[1][i])

       #Si le type zero existe, on le place a la fin: sinon adressage openmp=boom dans donorPts
       if 0 in Nbtype: Nbtype += [Nbtype.pop( Nbtype.index( 0 ) )]

       param_int[ iadr +rac[pos]*2 ] = len(Nbtype)

       param_int[ iadr +rac[pos]*3 ] = -1
       size_IBC = 0

       zsrname = s[0]
       sname = zsrname[0:2]
       xc=None;yc=None;zc=None; xi=None;yi=None;zi=None; xw=None;yw=None;zw=None;density=None;pressure=None;
       vx=None; vy=None; vz=None;
       utau=None;yplus=None;
       ptxc=0;ptyc=0;ptzc=0;ptxi=0;ptyi=0;ptzi=0;ptxw=0;ptyw=0;ptzw=0;ptdensity=0;ptpressure=0; 
       ptvx=0;ptvy=0;ptvz=0;
       ptutau=0;ptyplus=0;
       sd1=None;sd2=None;sd3=None;sd4=None;sd5=None
       ptd1=0;ptd2=0;ptd3=0;ptd4=0;ptd5=0
       if sname == 'IB': 
           zsrname = zsrname.split('_')
           if len(zsrname) < 3: 
                #print('Warning: miseAPlatDonorTree: non consistent with the version of IBM preprocessing.')
                if model=='Euler': 
                    print('Assuming IBC type is wallslip.')
                    param_int[iadr+rac[pos]*3]  = 0
                else: 
                    print('Assuming IBC type is Musker wall model.')
                    param_int[iadr+rac[pos]*3]  = 3
           else: param_int[iadr+rac[pos]*3]  = int(zsrname[1]) # 'IBCD_type_zonename'
          
           # print 'IBCType = ', IBCTypes[param_int[ iadr +rac[pos]*3 ]]
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

       tmp = Internal.getNodeFromName1(s, 'ZoneRole')
       if tmp[1][0] == b'D': param_int[ iadr +rac[pos]*4 ] = 0   # role= Donor
       else                : param_int[ iadr +rac[pos]*4 ] = 1   # role= Receiver
            
       param_int[ iadr +rac[pos]*5 ] = NozoneD                    # No zone donneuse

       lst                              = iadr + 1 -nb_rac[pos] +ntab_int*rac[pos] + sizeNbD[pos] + sizeType[pos] + size_ptlist[pos]  
       param_int[ iadr +rac[pos]*6    ] = lst                                                                      # PointlistAdr
       ptTy                             = iadr + 1 -nb_rac[pos] +ntab_int*rac[pos] + sizeNbD[pos] + size_ptType[pos]
       param_int[ iadr +rac[pos]*7    ] = ptTy                                                                     # TypAdr
       lstD                             = iadr + 1 -nb_rac[pos] +ntab_int*rac[pos] + size_ptlistD[pos]
       param_int[ iadr +rac[pos]*12 +1] = lstD                                                                     # PointlistDAdr

       Nbtot += Nbpts

       param_int[ ptTy  ] = len(Nbtype)
       noi       = 0
       nocoef    = 0
       sizecoef  = 0
       shift_typ = 1 + len(Nbtype)
       ctyp      = 0
       l0        = 0

       #recopie dans tableau a plat + tri par type
       if len(Nbtype) == 1:
           triMonoType(Nbpts_D, Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                       Interptype, pointlist, pointlistD, param_int,
                       ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, 
                       ptdensity,ptpressure,
                       ptvx, ptvy, ptvz,
                       ptutau,ptyplus,
                       ptd1,ptd2,ptd3,ptd4,ptd5,
                       xc,yc,zc,xi,yi,zi,xw,yw,zw, 
                       density,pressure,
                       vx, vy, vz,
                       utau,yplus,
                       sd1,sd2,sd3,sd4,sd5,
                       InterpD,param_real)
       else:
           triMultiType(Nbpts_D,Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                        Interptype, pointlist, pointlistD, param_int,
                        ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, 
                        ptdensity,ptpressure,
                        ptvx, ptvy, ptvz,
                        ptutau,ptyplus,
                        ptd1,ptd2,ptd3,ptd4,ptd5, 
                        xc,yc,zc,xi,yi,zi,xw,yw,zw, 
                        density,pressure,
                        vx, vy, vz,
                        utau,yplus,
                        sd1,sd2,sd3,sd4,sd5,
                        InterpD,param_real)

       pointlist[ 1] = param_int[ lst             : lst              + Nbpts         ]    # supression numpy initial pointlist
       Interptype[1] = param_int[ ptTy + shift_typ: ptTy + shift_typ + Nbpts_D       ]    # supression numpy initial interpolantType
       pointlistD[1] = param_int[ lstD            : lstD             + Nbpts_D       ]    # supression numpy initial pointlistDonor
       InterpD[   1] = param_real[ pt_coef        : pt_coef          + Nbpts_InterpD ]    # supression numpy initial interpDonor

       #if (s[0] == 'ID_cart3' and z[0]=='cart1'): print 'verif',  InterpD[   1][0], pt_coef,numpy.shape(InterpD[ 1  ])

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

           if sd1 is not None:
               sd1[1]  = param_real[ ptd1 : ptd1 + Nbpts_D ]
           if sd2 is not None:
               sd2[1]  = param_real[ ptd2 : ptd2 + Nbpts_D ]
           if sd3 is not None:
               sd3[1]  = param_real[ ptd3 : ptd3 + Nbpts_D ]
           if sd4 is not None:
               sd4[1]  = param_real[ ptd4 : ptd4 + Nbpts_D ]
           if sd5 is not None:
               sd5[1]  = param_real[ ptd5 : ptd5 + Nbpts_D ]
               
       param_int[ iadr +rac[pos]*8 ] = adr_coef[pos] + size_coef[pos]          # PtcoefAdr
          
       iadr = iadr +1
       param_int[ iadr +rac[pos]*8 ] = rac[pos]                  # nrac pour mpi

 
       tmp = Internal.getNodeFromName1(s , 'GridLocation')
       if tmp[1][4] == b'C': param_int[ iadr +rac[pos]*9 ] = 1   # location= CellCenter
       else                : param_int[ iadr +rac[pos]*9 ] = 0   # location= Vertex

       param_int[ iadr +rac[pos]*10 ] = Nbpts_D

       #chercher No zone receveuse grace a tc ou dico (si mpi)
       if procDict is None:
         param_int[ iadr +rac[pos]*11 ] = znd.index( zRname )         # No zone receveuse
       else:
         param_int[ iadr +rac[pos]*11  ]= procList[proc].index( zRname )  # No zone raccord

       #print 'rac', s[0], 'zoneR=', zRname, param_int[ iadr +rac[pos]*11 ], 'NozoneD=', zones_tc[No_zoneD[c]][0]

       # ATTENTION STEPHANIE A TOUT COMMENTE
       #a = Internal.getNodeFromName2(zones_tc[ param_int[ iadr +rac[pos]*11  ]  ], 'GoverningEquations')
       #if a is not None: model = Internal.getValue(a)

       #print 'model=',model,'zoneR',zones_tc[param_int[ iadr +rac[pos]*11  ]][0], 'NoR=', param_int[ iadr +rac[pos]*11  ], 'NoD=', c
       if model=='NSTurbulent': neq_loc = 6
       else                   : neq_loc = 5
            
       tmp =  Internal.getNodeFromName1(s , 'RANSLES')
       if tmp is not None: param_int[ iadr +rac[pos]*13  ] = min (5, neq_loc)   # RANSLES
       else:               param_int[ iadr +rac[pos]*13  ] = neq_loc
       
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

       #print 'rac*14= ',iadr +rac[pos]*14 

       # raccord instationnaire
       param_int[ iadr +rac[pos]*15  ] =-1
       if '#' in s[0]:
           numero_iter = int( s[0].split('#')[1].split('_')[0] )
           param_int[ iadr +rac[pos]*15  ] = numero_iter

       if nrac != 0:    

         if pranged is not None and prange is not None : #Si on est en explicite local, on rajoute des choses dans param_int
         
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] : iadr2 + rac[pos]*16 + 27*nb_rac[pos]+6 ] = prange[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos]  + 6] = direction[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 7 : iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 13 ] = pranged[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 13] = directiond[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 14 : iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 17 ] = transfo[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 17 : iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 20 ] = pt_pivot[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 20 : iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 21 ] = profondeur[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 21 : iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 24 ] = ratio[1]
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 24] = int(levelrcv[1])
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 25] = int(leveldnr[1])
           param_int[ iadr2 + rac[pos]*16 + 27*nb_rac[pos] + 26] = S

           if (int(levelrcv[1]) != int(leveldnr[1])):
               S=S + 3*5*(pranged[1][1]-pranged[1][0]+1)*(pranged[1][3]-pranged[1][2]+1)*(pranged[1][5]-pranged[1][4]+1)*(profondeur[1][0]+1)
               #print S
           #print [int(levelrcv[1])], [int(leveldnr[1])],'--'
           #print int([1]
       #print 'model=', model, 'tmp=', tmp, 'neq_loc=', neq_loc

       ### Verifier position et chois entre Nbpts et NbptsD 
       size_ptlistD[pos] = size_ptlistD[pos] + Nbpts_D
       size_ptlist[pos]  = size_ptlist[pos]  + Nbpts
       size_ptType[pos]  = size_ptType[pos]  + Nbpts_D + len(Nbtype)+1

       size_coef[pos] = size_coef[pos] + Nbpts_InterpD + size_IBC + shiftRotation
       nb_rac[pos]    = nb_rac[pos] + 1
       
       c += 1
    
    return None

#==============================================================================
# tri multitype
#==============================================================================
def triMultiType(Nbpts_D, Nbpts, Nbpts_InterpD,meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                 Interptype, pointlist, pointlistD, param_int,
                 ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, 
                 ptdensity,ptpressure,
                 ptvx, ptvy, ptvz,
                 ptutau,ptyplus,
                 ptd1,ptd2,ptd3,ptd4,ptd5,
                 xc,yc,zc,xi,yi,zi,xw,yw,zw, 
                 density,pressure,
                 vx, vy, vz,
                 utau,yplus,
                 sd1,sd2,sd3,sd4,sd5,
                 InterpD,param_real):
  for ntype in Nbtype:
    noi_old   = 0
    nocoef_old= 0
    l         = 0
    for i in range(Nbpts_D):
       ltype = Interptype[1][i]
       if meshtype == 1:
         if ltype == 1: sizecoef=1
         elif ltype == 2: sizecoef=8
         elif ltype == 3: sizecoef=9
         elif ltype == 4: sizecoef=8
         elif ltype == 5: sizecoef=15
         elif ltype == 22: sizecoef=4
       else:
         if ltype == 1: sizecoef=1
         elif ltype == 4: sizecoef=4

       if ltype == ntype:
            # recopieinterpolantType
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

               if sd1 is not None:
                   param_real[ ptd1   + l + l0 ]= sd1[1][i]
                   param_real[ ptd2   + l + l0 ]= sd2[1][i]
                   param_real[ ptd3   + l + l0 ]= sd3[1][i]
                   param_real[ ptd4   + l + l0 ]= sd4[1][i]
                   param_real[ ptd5   + l + l0 ]= sd5[1][i]
                   
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
    

  return None
#==============================================================================
# tri monotype
#==============================================================================
def triMonoType(Nbpts_D, Nbpts, Nbpts_InterpD,meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                Interptype, pointlist, pointlistD, param_int,
                ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, 
                ptdensity,ptpressure,
                ptvx, ptvy, ptvz,
                ptutau,ptyplus,
                ptd1,ptd2,ptd3,ptd4,ptd5,
                xc,yc,zc,xi,yi,zi,xw,yw,zw, 
                density,pressure,
                vx, vy, vz,
                utau,yplus,
                sd1,sd2,sd3,sd4,sd5,
                InterpD, param_real):

  ntype     = Nbtype[0]
  noi_old   = 0
  nocoef_old= 0
  l         = 0
  ltype     = Interptype[1][0]

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
           
  # recopie  InterpD
  connector.initNuma(InterpD[1] , param_real, pt_coef , Nbpts_InterpD , 0, val)

  param_int[ ptTy + ctyp +1 ] = Nbpts_D

  return None

#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la zone donneuse
# fonctionne avec ___setInterpTransfer
#==============================================================================
def miseAPlatDonorZone__(zones, tc, procDict):
    zones_tc = Internal.getZones(tc)
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
            sd1          =  Internal.getNodeFromName1(rac, 'StagnationEnthalpy')

            ntab_IBC   = 11+3 #On ajoute dorenavant les vitesses dans l'arbre tc pour le post
            if utau is not None: ntab_IBC += 2
            if sd1 is not None: ntab_IBC += 5

            Nbpts        =  numpy.shape(pointlist[ 1])[0]
            Nbpts_D      =  numpy.shape(pointlistD[1])[0]
            Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

            size_int   =  size_int + 7 + Nbpts_D*2 + Nbpts
            size_real  =  size_real+   Nbpts_InterpD
            sname = rac[0][0:2]
            if sname == 'IB': 
               size_real = size_real +Nbpts_D*ntab_IBC
               count_IBC = count_IBC +1;
            else: 
               count_ID  = count_ID  +1;
            #print 'nbpt, nbpt_donor', sname,Nbpts,Nbpts_InterpD

        size_int = size_int + 2 + (count_IBC + count_ID)*2  # 2: nbr rac ID et IBC, stockage adresse debut raccord et coef
        param_int  = numpy.empty(size_int , dtype=numpy.int32  )
        param_real = numpy.empty(size_real, dtype=numpy.float64)
        Internal.createUniqueChild(z, 'Parameter_int' , 'DataArray_t', param_int )
        if size_real !=0 :
           Internal.createUniqueChild(z, 'Parameter_real', 'DataArray_t', param_real)# recopie pointlis

        #print 'size int et real', size_int, size_real
        #print 'ID et IBC', count_ID, count_IBC

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
            no_zone =0
            for z in zones:
               if z[0] == zrcvname: param_int[ pt_rac +6 ]= no_zone  # No zone raccord                    
               no_zone = no_zone +1

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
                else:
                   var_ibc=['CoordinateX_PC','CoordinateY_PC','CoordinateZ_PC','CoordinateX_PI','CoordinateY_PI','CoordinateZ_PI','CoordinateX_PW','CoordinateY_PW','CoordinateZ_PW', 'Density','Pressure','VelocityX','VelocityY','VelocityZ']

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
