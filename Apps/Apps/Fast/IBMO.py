# Class for FastS "IBM"+"Overset" prepare and compute
import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Initiator.PyTree as I
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
from Apps.Fast.Common import Common

# return True si z est une zone chimera
def isZoneChimera(z):
    if z[3] != 'Zone_t': print('Warning: isZoneChimera: not a Zone node.')
    n = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
    if n is None: return False
    n = Internal.getNodeFromType1(n, 'GridConnectivity_t')
    if n is None: return False
    n = Internal.getNodeFromType1(n, 'GridConnectivityType_t')
    if n is None: return False
    if Internal.getValue(n) == 'Overset': return True
    return False

#================================================================================
# IBMO prepare
#
#================================================================================
def prepare(t_case, t_out, tc_out,   
            vmin=21, check=False, NP=0,
            frontType=1, expand=3):
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    # Extraction de la liste des dfars de tb
    zones = Internal.getZones(tb)
    dfarList = [10.]*len(zones)
    for c, z in enumerate(zones): 
        n = Internal.getNodeFromName2(z, 'dfar')
        if n is not None: dfarList[c] = Internal.getValue(n)*1.

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree defined in %s.'%FILE)
    model = Internal.getValue(model)

    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced
    
    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    # constuction de l'arbre des corps pour l'octree : tbo
    tbo = Internal.copyRef(tb)
    bases = Internal.getBases(tb)
    for b in bases:
        c = 0
        for z in b[2]:
            if z[3] == 'Zone_t':
                if isZoneChimera(z): 
                    ext = C.extractBCOfType(z, 'BCOverlap')
                    if len(ext)>0: b[2][c] = ext[0]
                    if len(ext)>1: b[2] += ext[1:]
            c += 1
            
    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank,
                         expand=expand)

    return tb

#====================================================================================
class IBMO(Common):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["stephanie@onera.fr", "ash@onera.fr"]
        self.cartesian = True
        
    # Prepare 
    def prepare(self, t_case, t_out, tc_out, 
                vmin=21, check=False, frontType=1, NP=None, expand=3):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for a computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, 
                      vmin=vmin, check=check, NP=NP,  
                      frontType=frontType, expand=expand)
        return ret
