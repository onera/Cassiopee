"""Extract doubly defined boundaries."""
# A doubly defined BC extractor

# IN: maillage volumique + BC walls + BC Overlap doubly defined
# Utilise setDoublyDefined -> cellN modifie
# Recopie le cellN proche aux centres dans le BCDataSet
# Extrait les BCs
# Fait un selectCells

import Converter.Internal as Internal
import Converter.PyTree as C
import Transform.PyTree as T
import Connector.PyTree as X
import Post.PyTree as P
from Apps.App import App

def _extrapOnBCDataSet(t, variables):
    # variables are extrapolated on BCDataSet **if not already present**
    for z in Internal.getZones(t):
        for bc in Internal.getNodesFromType2(z, 'BC_t'):
            allbcdset = Internal.getNodesFromType1(bc, 'BCDataSet_t')

            for bcdset in allbcdset:
                for bcd in Internal.getNodesFromType1(bcdset, 'BCData_t'):
                    fields = Internal.getNodesFromType1(bcd, 'DataArray_t')

                    for var in variables:
                        v = var.split(':',1)
                        if len(v) == 2:
                            if v[0] == 'centers' or v[0] == 'nodes': v = v[1]
                            else: v = var
                        else: v = var
                        found = 0
                        for field in fields:
                            if Internal.isName(field, v): found=1; break
                        if found == 0:
                            _setBCDataSet__(z, bc, [var])
            if allbcdset == []:
                for var in variables:
                    _setBCDataSet__(z, bc, [var])

def _setBCDataSet__(z, bc, variables):
    # IndexRange
    r = Internal.getNodeFromType1(bc, 'IndexRange_t')
    if r is not None and r[1].shape[0] > 1: rangew = r[1]
    w = Internal.range2Window(rangew)

    imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
    zw = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
    # zw = CP.subzoneWithReorder__(z, w) # Modif CW 04/11/16
    datas = Internal.getBCDataSet(z,bc)
    FSC = Internal.getNodeFromName1(zw, Internal.__FlowSolutionCenters__)
    FSN = Internal.getNodeFromName1(zw, Internal.__FlowSolutionNodes__)
    if datas == []:
        d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                  gridLocation='FaceCenter', parent=bc)
        d = Internal.newBCData('BCNeumann', parent=d)
        if variables == []:
            if FSC is not None:
                for fsc in FSC[2]:
                    if fsc[3] == 'DataArray_t':
                        varName = fsc[0]
                        Internal._createUniqueChild(d,varName,'DataArray_t',value=fsc[1])
            if FSN is not None:
                for fsn in FSN[2]:
                    if fsc[3] == 'DataArray_t':
                        varName = fsn[0]
                        Internal._createUniqueChild(d,varName,'DataArray_t',value=fsn[1])
        else:
            for var in variables:
                varl = var.split(':',1)
                if len(varl) > 1:
                    if varl[0] == 'centers':
                        varname = varl[1]
                        fs = Internal.getNodeFromName1(FSC,varname)
                    elif varl[0] == 'nodes':
                        varname = varl[1]
                        fs = Internal.getNodeFromName1(FSN,varname)
                    else:
                        varname = var
                        fs = Internal.getNodeFromName1(FSN,varname)
                else:
                    varname = var
                    fs = Internal.getNodeFromName1(FSN,varname)
                Internal._createUniqueChild(d,varname,'DataArray_t',value=fs[1])

    else:
        cont, c = Internal.getParentOfNode(z, datas[0])
        if variables == []:
            for fsc in FSC:
                varName = fsc[0]
                Internal._createUniqueChild(cont,varName,'DataArray_t',value=fsc[1])
            for fsn in FSN:
                varName = fsn[0]
                Internal._createUniqueChild(cont,varName,'DataArray_t',value=fsn[1])
        else:
            for var in variables:
                varl = var.split(':',1)
                if len(varl) > 1:
                    if varl[0] == 'centers':
                        varname = varl[1]
                        fs = Internal.getNodeFromName1(FSC,varname)
                    else:
                        varname = varl[1]
                        fs = Internal.getNodeFromName1(FSN,varname)
                else:
                    varname = varl[0]
                    fs = Internal.getNodeFromName1(FSN,varname)
                Internal._createUniqueChild(cont, varname, 'DataArray_t', value=fs[1])
    return None

def extract(t, bcType, compute):
    """Extract doubly defined boundaries."""
    if compute:
        # Calculs des doubly defined -> cellN
        C._initVars(t, 'centers:cellN', 1)
        t = X.applyBCOverlaps(t)
        t = X.setDoublyDefinedBC(t)
    _extrapOnBCDataSet(t, ['centers:cellN'])
    bcs = C.extractBCOfType(t, bcType)
    p = P.selectCells(bcs, '{centers:cellN}==1')
    return p

#====================================================================================
class ExtractDoublyDefined(App):
    """Extraction of boundaries taking doubly defined boundaries into account."""
    def __init__(self):
        App.__init__(self)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.requires(['dataIn', 'dataOut', 'bcType', 'compute'])
        # default values
        self.set(bcType='BCWall')
        self.set(compute=True)

    def run(self):
        self.check()
        t = self.readDataIn()
        p = extract(t, self.data['bcType'], self.data['compute'])
        self.writeDataOut(p)
