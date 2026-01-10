# Class common to all types of FAST simulations
import FastC.PyTree as FastC
import Converter.PyTree as C
import Converter.Internal as Internal
from Apps.App import App
import Converter.Filter as Filter
import Converter.Mpi as Cmpi

#================================================================================
# Redistribue les fichiers in place sans com pour l'instant
# Change les noeuds procs seulement
#=================================================================================
def _distribute(t_in, tc_in, NP=Cmpi.size):
    if isinstance(t_in, str):
        if Cmpi.rank == 0: _distributeFile(t_in, tc_in, NP)
    else: _distributeMem(t_in, tc_in, NP)
    return None

def _distributeFile(t_in, tc_in, NP):
    import Distributor2.PyTree as D2
    t = Filter.convertFile2SkeletonTree(t_in, maxDepth=3, maxFloatSize=6)
    tc = Filter.convertFile2SkeletonTree(tc_in, maxDepth=3, maxFloatSize=6)
    stats = D2._distribute(tc, NP, algorithm='graph', useCom='ID')
    D2._copyDistribution(t, tc)
    nodes = Internal.getNodesFromName(t, 'proc')
    for n in nodes:
        p = Internal.getPath(t, n)
        Filter.writeNodesFromPaths(t_in, p, n)
    nodes = Internal.getNodesFromName(tc, 'proc')
    for n in nodes:
        p = Internal.getPath(tc, n)
        Filter.writeNodesFromPaths(tc_in, p, n)
    return None

def _distributeMem(t, tc, NP):
    import Distributor2.PyTree as D2
    tbbc = Cmpi.createBBoxTree(tc)
    stats = D2._distribute(tbbc, NP, algorithm='graph', useCom='ID')
    D2._copyDistribution(tc, tbbc)
    D2._copyDistribution(t, tbbc)
    return None

#==================================================
# distribution + choix du nombre de coeurs optimum
#==================================================
def _distributeOpt(t_in, tc_in, corePerNode=28, nptMaxPerCore=4.e6):
    if isinstance(t_in, str):
        if Cmpi.rank == 0: _distributeOptFile(t_in, tc_in, corePerNode, nptMaxPerCore)
    else:
        #_distributeOptMem(t_in, tc_in, corePerNode, nptMaxPerCore)
        raise ValueError('Not implemented.')
    return None

def _distributeOptFile(t_in, tc_in, corePerNode=28, nptMaxPerCore=4.e6):
    import Distributor2.PyTree as D2
    t = Filter.convertFile2SkeletonTree(t_in, maxDepth=3, maxFloatSize=6)
    tc = Filter.convertFile2SkeletonTree(tc_in, maxDepth=3, maxFloatSize=6)

    nbpts=0.; maxipts=0.
    for zone in Internal.getZones(t):
        ncells = C.getNCells(zone)
        nbpts += ncells
        maxipts = max(maxipts, ncells)

    MaxNbProcs=int(nbpts/maxipts)+1

    MinNbProcs=MaxNbProcs
    for nbproc in range(2,MaxNbProcs+1):
        if nbpts*1./nbproc < nptMaxPerCore: MinNbProcs=min(nbproc,MinNbProcs)

    print('La distribution sera testee entre %d procs et %d procs.'%(MinNbProcs,MaxNbProcs))

    #MinNbProcs = 140
    #MaxNbProcs = 140
    listequ = []
    varmax = 99.
    NP = MinNbProcs
    for nbproc in range(MinNbProcs,MaxNbProcs+1):
        if nbproc%corePerNode==0:
            print('Distribution sur %s procs.')
            stats = D2._distribute(tc, nbproc, algorithm='graph', useCom='ID')
            listequ.append([nbproc,stats['varMax']])
            if stats['varMax']<varmax: varmax = stats['varMax']; NP=nbproc

    stats = D2._distribute(tc, NP, algorithm='graph', useCom='ID')
    print('Best distribution found on %d procs.'%NP)
    print(stats)
    #D2._printProcStats(t,stats,NP)

    D2._copyDistribution(t, tc)
    nodes = Internal.getNodesFromName(t, 'proc')
    for n in nodes:
        p = Internal.getPath(t, n)
        Filter.writeNodesFromPaths(t_in, p, n)
    nodes = Internal.getNodesFromName(tc, 'proc')
    for n in nodes:
        p = Internal.getPath(tc, n)
        Filter.writeNodesFromPaths(tc_in, p, n)

    return NP

#================================================================================
# en gros, warmup
#================================================================================
def setup(t_in, tc_in, numb, numz, format='single'):
    if Cmpi.size > 1:
        import FastS.Mpi as FastS
        FastC.HOOK = None
        rank = Cmpi.rank; size = Cmpi.size
    else:
        import FastS.PyTree as FastS
        FastC.HOOK = None
        rank = 0; size = 1

    t,tc,ts,graph = FastC.load(t_in, tc_in, split=format)

    # Numerics
    FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)
    (t, tc, metrics) = FastS.warmup(t, tc, graph)
    return t, tc, ts, metrics, graph

#============================================================================
# Ecrit le resultat
# t: arbre
# t_out: fichier de sortie
# it0: iteration correspondant a la fin du cacul
# time0: temps correspondant a la fin du calcul
# format: "single" ou "multiple"
# compress: si 1, compress le fichier pour le cartesien, 2, compressAll
# ===========================================================================
def finalize(t, t_out=None, it0=None, time0=None, format='single', compress=0):
    if it0 is not None:
        Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0)
    if time0 is not None:
        Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
    if t_out is not None and isinstance(t_out, str):
        FastC.save(t, t_out, split=format, compress=compress)

#=====================================================================================
# IN: t_in : nom du fichier t input ou arbre input
# IN: tc_in: nom du fichier tc input ou arbre tc
# IN: t_out, tc_out: noms de fichiers ou arbres de sortie
# IN: numb, numz: les data numeriques
# IN: NIT: nbre d'iterations
# format: single ou multiple
# compress: si 1, compress le fichier de sortie pour le cartesien, 2 compressAll
#======================================================================================
def compute(t_in, tc_in,
            t_out, tc_out,
            numb, numz,
            NIT,
            format='single', compress=0):
    if Cmpi.size > 1:
        import FastS.Mpi as FastS
        rank = Cmpi.rank; size = Cmpi.size
        tc, graph = FastC.loadTree(tc_in, graph=True, split=format)
    else:
        import FastS.PyTree as FastS
        rank = 0; size = 1 ; graph=None
        tc = FastC.loadTree(tc_in, split=format)

    t = FastC.loadTree(t_in, split=format)

    # Numerics
    FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

    (t, tc, metrics) = FastS.warmup(t, tc, graph)

    it0 = 0; time0 = 0.
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'Time')
    if first is not None: time0 = Internal.getValue(first)
    time_step = Internal.getNodeFromName(t, 'time_step')
    time_step = Internal.getValue(time_step)

    if 'modulo_verif' in numb: moduloVerif = numb['modulo_verif']
    else: moduloVerif = 200

    for it in range(NIT):
        FastS._compute(t, metrics, it, tc, graph)
        if it%moduloVerif == 0:
            if rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
            FastS.display_temporal_criteria(t, metrics, it, format='double')
        #if it%50 == 0:
        #    import CPlot.PyTree as CPlot
        #    CPlot.display(t, dim=2, mode='Scalar', scalarField='Density')
        time0 += time_step

    # time stamp
    Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
    Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
    if t_out is not None and isinstance(t_out,str):
        FastC.save(t, t_out, split=format, compress=compress)
    if tc_out is not None and isinstance(tc_out,str):
        FastC.save(tc, tc_out, split=format)
    if Cmpi.size > 1: Cmpi.barrier()
    return t, tc

#===============================================================================
class Common(App):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        App.__init__(self)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.requires(['format', 'numb', 'numz'])
        self.compress = 0
        # default values
        if format is not None: self.set(format=format)
        else: self.set(format='single')
        if numb is not None: self.set(numb=numb)
        if numz is not None: self.set(numz=numz)

    # Compute nit iterations
    # peut etre lance en sequentiel ou en parallele
    def compute(self, t_in, tc_in, t_out, nit, tc_out=None):
        numb = self.data['numb']
        numz = self.data['numz']
        compress = self.compress
        return compute(t_in, tc_in, t_out, tc_out,
                       numb, numz,
                       nit,
                       format=self.data['format'],
                       compress=compress)

    # warm up et all
    def setup(self, t_in, tc_in):
        numb = self.data['numb']
        numz = self.data['numz']
        return setup(t_in, tc_in, numb, numz, self.data['format'])

    # Ecrit le fichier de sortie
    def finalize(self, t_out=None, it0=None, time0=None):
        finalize(t_out, it0, time0, self.data['format'])

    # distribue fichiers ou en memoire
    def _distribute(self, t_in, tc_in, NP=Cmpi.size):
        return _distribute(t_in, tc_in, NP)

    # distribue fichiers ou en memoire, trouve le NP optimal
    def _distributeOpt(self, t_in, tc_in, corePerNode=28, nptMaxPerCore=4.e6):
        return _distributeOpt(t_in, tc_in, corePerNode, nptMaxPerCore)
