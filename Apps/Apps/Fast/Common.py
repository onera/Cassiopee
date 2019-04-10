# Class common to all types of FAST simulations
import Fast.PyTree as Fast
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
from Apps.App import App

try: range = xrange
except: pass

#================================================================================
# Redistribue un fichier in place sans com pour l'instant
# Change les noeuds procs seulement
#================================================================================
def distribute(t_in, NP):
    import Filter
    t  = Filter.convertFile2SkeletonTree(t_in, maxDepth=2,maxSize=6)
    t, stats = D2.distribute(t, NP, algorithm='graph', useCom=None)
    nodes = Internal.getNodesFromName(t, 'Proc')
    for n in nodes:
        p = Internal.getPath(n)
        Filter.writeNodesFromPaths(t_in, p, n)
    return t

# distribution + choix du nombre de coeurs optimum
def distributeOpt(t_in, tc_in, t_out, tc_out):
    if instance(t_in, str): t = C.convertFile2PyTree(t_in)
    else: t = t_in
    nptMaxPerCore = 4000000.

    nbpts=0.; maxipts=0.
    for zone in Internal.getZones(t):
        ncells = C.getNCells(zone)
        nbpts += ncells
        maxipts = max(maxipts, ncells)

    MaxNbProcs=int(nbpts/maxipts)+1

    MinNbProcs=MaxNbProcs
    for nbproc in range(2,MaxNbProcs+1):
        if nbpts/nbproc < nptMaxPerCore: MinNbProcs=min(nbproc,MinNbProcs)

    print('La distribution sera testee entre %d procs et %d procs.'%(MinNbProcs,MaxNbProcs))

    #MinNbProcs = 140
    #MaxNbProcs = 140
    listequ=[]
    varmax = 99.
    for nbproc in range(MinNbProcs,MaxNbProcs+1):
        if nbproc%28==0:
         print('Distribution sur %s procs.')
         stats=D2._distribute(t,nbproc,algorithm='fast')
         listequ.append([nbproc,stats['varMax']])
         if stats['varMax']<varmax: varmax = stats['varMax'];NP=nbproc

    import pprint; print(pprint.pformat(listequ))

    #import sys;sys.exit()
    #NP = 532
    stats = D2._distribute(t, NP)
    print(stats)
    print(NP)
    #D2._printProcStats(t,stats,NP)

    C.convertPyTree2File(t, t_out)
    if instance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    
    D2._copyDistribution(tc, t)
    C.convertPyTree2File(tc, tc_out)
    return t, tc, NP

#================================================================================
# en gros, warmup
#================================================================================
def setup(t_in, tc_in, numb, numz, NP=0, format='single'):
    if NP > 0:
        import Converter.Mpi as Cmpi
        import FastS.Mpi as FastS
        rank = Cmpi.rank; size = Cmpi.size
    else:
        import FastS.PyTree as FastS
        rank = 0; size = 1

    if NP != 0 and NP != size:
        raise ValueError, 'setup: you are running not on the prepared number of processors %d != %d'%(NP, size)

    t,tc,ts,graph = Fast.load(t_in, tc_in, split=format, NP=NP)

    # Numerics
    Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)
    (t, tc, metrics) = FastS.warmup(t, tc, graph)
    return t, tc, ts, metrics, graph

#============================================================================
# Ecrit le resultat
# t: arbre
# t_out: fichier de sortie
# it0: iteration correspondant a la fin du cacul
# time0: temps correspondant a la fin du calcul
# ===========================================================================
def finalize(t, t_out, it0=None, time0=None, NP=0, format='single'):
    if it0 is not None:
        Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0)
    if time0 is not None:
        Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
    Fast.save(t, t_out, split=format, NP=NP)

#=====================================================================================
# NP is the currently running number of processors
# IN: file names
#======================================================================================
def compute(t_in, tc_in, 
            t_out, tc_out,
            numb, numz,
            NIT, 
            NP=0, format='single'):
    if NP > 0:
        import Converter.Mpi as Cmpi
        import FastS.Mpi as FastS
        rank = Cmpi.rank; size = Cmpi.size
    else:
        import FastS.PyTree as FastS
        rank = 0; size = 1

    if NP != 0 and NP != size:
        raise ValueError('compute: you are running not on the prepared number of processors %d != %d'%(NP, size))

    t,tc,ts,graph = Fast.load(t_in, tc_in, split=format, NP=NP)

    # Numerics
    Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

    (t, tc, metrics) = FastS.warmup(t, tc, graph)

    it0 = 0; time0 = 0.
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'Time')
    if first is not None: time0 = Internal.getValue(first)
    time_step = Internal.getNodeFromName(t, 'time_step')
    time_step = Internal.getValue(time_step)

    for it in range(NIT):
        FastS._compute(t, metrics, it, tc, graph)
        if it%100 == 0:
            if rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        #FastS.display_temporal_criteria(t, metrics, it, format='double')
        time0 += time_step

    # time stamp
    Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
    Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
    Fast.save(t, t_out, split=format, NP=NP)
    if tc_out is not None: Fast.save(tc, tc_out, split=format, NP=NP)
    if NP > 0: Cmpi.barrier()
    return t, tc

#===============================================================================
class Common(App):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, NP=None, format=None, numb=None, numz=None):
        App.__init__(self)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.requires(['NP', 'format', 'numb', 'numz'])
        # default values
        if NP is not None: self.set(NP=NP)
        else: self.set(NP=0)
        if format is not None: self.set(format=format)
        else: self.set(format='single')
        if numb is not None: self.set(numb=numb)
        if numz is not None: self.set(numz=numz)

    # Compute nit iterations
    # peut etre lance en sequentiel ou en parallele
    def compute(self, t_in, tc_in, t_out, nit, tc_out=None):
        numb = self.data['numb']
        numz = self.data['numz']
        return compute(t_in, tc_in, t_out, tc_out,
                       numb, numz,
                       nit, 
                       NP=self.data['NP'], 
                       format=self.data['format'])

    # warm up et all
    def setup(self, t_in, tc_in):
        numb = self.data['numb']
        numz = self.data['numz']
        return setup(t_in, tc_in, numb, numz, self.data['NP'], self.data['format'])

    # Ecrit le fichier de sortie
    def finalize(self, t_out, it0=None, time0=None):
        finalize(t_out, it0, time0, self.data['NP'], self.data['format'])
