import Converter.Internal as Internal
import numpy
import Compressor
import Compressor.PyTree as Co

try: range = xrange
except: pass

# =============================================================================
# Initialization of dictionnary for unsteady preprocessing
# OUT : hook: list of dictionnaries for blanking, center and face interpolations
# =============================================================================
def initUnsteadyChimera(t):
    iterationB={}; iterationI={}; iterationIInt={}                          # first iteration of blanking, cell center interpolations, face center interpolations
    globalBlankingData={}; globalCellInterpData={}; globalFaceInterpData={} # data for blanking, cell center interpolations, face center interpolations
    refB={}; refI={}; refIInt={}                                            # reference for  blanking, cell center interpolations, face center interpolations
    zoneId={}                                                               # zoneId: elsA Id of zones (zones order in the CGNSTree)
    alreadyBlanked={}                                                       # alreadyBlanked: tells if key zone is already blanked or not

    # initialize dictionaries
    zones= Internal.getNodesFromType(t,'Zone_t')
    Id = 0
    for z in zones: zoneId[z[0]] = Id; Id = Id+1
    alreadyBlanked={}
    for z in zones:
        refB[z[0]] = None
        globalBlankingData[z[0]] = []
        alreadyBlanked[z[0]] = 0

    # build hook
    hook=[iterationB, iterationI, iterationIInt,globalBlankingData,globalCellInterpData,globalFaceInterpData,refB, refI, refIInt,zoneId, alreadyBlanked]
    return hook

#==============================================================================
# Compute unsteady interpolations
# IN: tp: pyTree
# IN: ite: current iteration
# IN: loc: localization of interpolation (cell or face)
# IN: nGhostCells: number of ghost cells by direction to take into account
# IN/OUT: hook
#==============================================================================
def computeUnsteadyInterp(tp, hook, ite,loc='cell', nGhostCells=2):
    if loc == 'cell':
        ListDonor = 'PointListDonor'
        ListExtC = 'PointListExtC'
        InterpolantsDonor = 'InterpolantsDonor'
        InterpolantsType = 'InterpolantsType'
        FaceDirection=None
        iteration = hook[1]
        listInterpData = hook[4]
        ref = hook[7]
    else:
        ListDonor = 'FaceListDonor'
        ListExtC = 'FaceListExtC'
        InterpolantsDonor = 'FaceInterpolantsDonor'
        InterpolantsType = 'FaceInterpolantsType'
        FaceDirection='FaceDirection'
        iteration = hook[2]
        listInterpData = hook[5]
        ref = hook[8]
    zones = Internal.getNodesFromType(tp,'Zone_t')
    for z in zones:
        donorName = z[0]
        dim = Internal.getZoneDim(z)
        subRegions = Internal.getNodesFromType(z, 'ZoneSubRegion_t')
        interpolations=[]
        interpData = {}
        for s in subRegions:
            if 'ID_' in s[0]:interpolations.append(s)
        # la zone a une ou plusieurs interpolations
        if interpolations != []:
            if donorName not in iteration:
                iteration[donorName] = ite+1
                listInterpData[donorName] = []
                ref[donorName] = {}
            for interp in interpolations:
                rcvName = '_'.join(interp[0].split('_')[1:])
                rcvId = hook[9][rcvName]
                # cells
                if Internal.getNodeFromName1(interp, ListExtC) is not None:
                    donorIndices = Internal.getNodeFromName1(interp, ListExtC)[1]; 
                    donorIndices = donorIndices.reshape((donorIndices.shape[0]))
                    if donorIndices.shape[0] != 0: # avoid interpolation regions with only orphan points
                        coefs = Internal.getNodeFromName1(interp, InterpolantsDonor)[1][:,0:7]
                        rcvIndices = Internal.getNodeFromName1(interp, ListDonor)[1]; rcvIndices = rcvIndices.reshape((rcvIndices.shape[0]))

                        periodicity =  Internal.getNodeFromName1(interp, InterpolantsType)[1]; periodicity= periodicity.reshape((periodicity.shape[0]))
                        if FaceDirection is not None:
                            faceDir = Internal.getNodeFromName1(interp, FaceDirection)[1]; faceDir= faceDir.reshape((faceDir.shape[0]))
                        # cell index => faceIndex
                        if FaceDirection is not None:
                            zRcv = Internal.getNodeFromName2(tp,rcvName) 
                            dimrcv = Internal.getZoneDim(zRcv)
                            imr = dimrcv[1];jmr = dimrcv[2];kmr = dimrcv[3]
                            imrg = imr-1 + 2*nGhostCells
                            jmrg = jmr-1 + 2*nGhostCells
                            kmrg = kmr-1 + 2*nGhostCells
                            nbintByDir = imrg*jmrg*kmrg
                            for i in range(len(rcvIndices)):
                                rk = rcvIndices[i]/((imr-1)*(jmr-1))
                                rj = (rcvIndices[i] - rk*(imr-1)*(jmr-1))/(imr-1)
                                ri = rcvIndices[i] - rk*(imr-1)*(jmr-1) - rj*(imr-1)+1
                                rj = rj+1; rk = rk+1
                                if faceDir[i] == 0: # imax (interface a droite)
                                    rcvIndices[i] = (rk-1 + nGhostCells)*imrg*jmrg+(rj-1 + nGhostCells)*imrg+ ri-1 + nGhostCells+1
                                elif faceDir[i] == 1: # imin (interface a gauche)
                                    rcvIndices[i] = (rk-1 + nGhostCells)*imrg*jmrg+(rj-1 + nGhostCells)*imrg+ ri-1 + nGhostCells
                                elif faceDir[i] == 2: # jmax
                                    rcvIndices[i] = (rk-1 + nGhostCells)*imrg*jmrg+(rj-1 + nGhostCells+1)*imrg+ ri-1 + nGhostCells + nbintByDir
                                elif faceDir[i] == 3: # jmin
                                    rcvIndices[i] = (rk-1 + nGhostCells)*imrg*jmrg+(rj-1 + nGhostCells)*imrg+ ri-1 + nGhostCells + nbintByDir
                                elif faceDir[i] == 4: # kmax
                                    rcvIndices[i] = (rk-1 + nGhostCells+1)*imrg*jmrg+(rj-1 + nGhostCells)*imrg+ ri-1 + nGhostCells + 2*nbintByDir
                                elif faceDir[i] == 5: # kmin
                                    rcvIndices[i] = (rk-1 + nGhostCells)*imrg*jmrg+(rj-1 + nGhostCells)*imrg+ ri-1 + nGhostCells + 2*nbintByDir
                        # First iteration of storage: full storage
                        if ite == 0 or listInterpData[donorName] == [] or rcvId not in listInterpData[donorName][-1].keys(): 
                            interpData[rcvId]={}
                            flag=0; i=0
                            if FaceDirection is None: # cell
                                for rcvIndex in rcvIndices:
                                    interpData[rcvId][(int)(rcvIndex)]=[flag,(int)(donorIndices[i]),(int)(periodicity[i])]+[(float)(c) for c in coefs[i]]; i = i+1
                            else: # face
                                for rcvIndex in rcvIndices:
                                    interpData[rcvId][(int)(rcvIndex)]=[flag,(int)(donorIndices[i]),(int)(periodicity[i])]+[(float)(c) for c in coefs[i]]+[(int)(faceDir[i])]; i = i+1                                
                        # delta storage
                        else: 
                            if FaceDirection is None: data=[rcvIndices,donorIndices,periodicity,coefs]
                            else: data=[rcvIndices,donorIndices,periodicity,coefs,faceDir]
                            delta = Co.deltaInterpolations(data, ref[donorName][rcvId],loc)
                            interpData[rcvId]=delta
                        # set reference of current iteration for next iteration
                        if FaceDirection is None: ref[donorName][rcvId] = [rcvIndices,donorIndices,periodicity, coefs]
                        else: ref[donorName][rcvId] = [rcvIndices,donorIndices,periodicity, coefs,faceDir]
            # Deals with rcv zones which have disappeared from interpolations
            for rcvId in ref[donorName]:
                if rcvId not in interpData:
                    if FaceDirection is None:
                        data=[numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= numpy.float64)]
                    else:
                        data=[numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= numpy.float64), numpy.array([], dtype= Internal.E_NpyInt)]
                    delta = Co.deltaInterpolations(data, ref[donorName][rcvId],loc)
                    interpData[rcvId]=delta
                    ref[donorName][rcvId] = data
            listInterpData[donorName].append(interpData)
        # la zone a aucune interpolation
        elif donorName in listInterpData:
            for rcvId in ref[donorName]:
                if FaceDirection is None:
                    data=[numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= numpy.float64)]
                else:
                    data=[numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= Internal.E_NpyInt), numpy.array([], dtype= numpy.float64), numpy.array([], dtype= Internal.E_NpyInt)]
                delta = Co.deltaInterpolations(data, ref[donorName][rcvId],loc)
                interpData[rcvId]=delta
                ref[donorName][rcvId] = data
            listInterpData[donorName].append(interpData)
    return

#==============================================================================
# Compute unsteady blanking
# IN: tp: pyTree
# IN: ite: current iteration
# IN/OUT: hook
#==============================================================================
def computeUnsteadyBlanking(tp, hook, ite):
    zones = Internal.getNodesFromType(tp,'Zone_t')
    for z in zones:
        dim = Internal.getZoneDim(z)
        hole = Internal.getNodesFromName(z, 'OversetHoles')
        if hole != []:
            if z[0] not in hook[0]: hook[0][z[0]] = ite+1
            list = Internal.getNodesFromName(hole[0], 'PointList')[0][1]
            index = globalIndex(dim, list,0) # Pour Cassiopee et elsA (pas de cellules fictives dans ces fichiers)
            if (hook[6][z[0]] is None): hook[6][z[0]] = index ; hook[3][z[0]].append(hook[6][z[0]])
            else: delta = Compressor.deltaIndex(index, hook[6][z[0]]) ; hook[6][z[0]] = index ; hook[3][z[0]].append(delta)
            hook[10][z[0]] = 1
        elif hook[10][z[0]] == 1:
            index=[0]
            delta = Compressor.deltaIndex(index, hook[6][z[0]]); hook[6][z[0]] = index ; hook[3][z[0]].append(delta)
    return

#==============================================================================
# Compute global indices, taking into account ghostcells if required
# ------------------------------------------------------------------
# IN: dim: zone dim
# IN: array1: index i,j,k
# OUT: array2: index ind = i + j*ni + k*ni*nj
#==============================================================================
def globalIndex(dim, array1,ghostcells):
    ni = dim[1]-1+2*ghostcells ; nj = dim[2]-1+2*ghostcells; nk = dim[3]-1+2*ghostcells ;
    nij = ni*nj
    s = array1.shape[1]
    a = numpy.empty( (s), dtype=Internal.E_NpyInt )
    for i in range(s):
        a[i] = (array1[0,i]-1+ghostcells) + (array1[1,i]-1+ghostcells)*ni + (array1[2,i]-1+ghostcells)*nij
    return a

#==============================================================================
# write files storing unsteady coefficients
# IN: hook
# IN: prefix: prefix of file name
# IN: total number of iterations
#==============================================================================
def writeUnsteadyCoefs(hook, prefix,nit,format="b"):
    # cells
    for name in hook[4]: # stocker le nom pour s y retrouver dans elsA
        if (hook[4][name] != []):
            if format == "b":
                filename = prefix+"_%d_%04d.bin"%(nit,hook[9][name])
            else:
                filename = prefix+"_%d_%04d.fmt"%(nit,hook[9][name])
            Compressor.writeUnsteadyCoefs(hook[1][name],hook[4][name], filename,"cell",format)
    # faces
    for name in hook[5]: # stocker le nom pour s y retrouver dans elsA
        if (hook[5][name] != []):
            if format == "b":
                filename = prefix+"_%d_%04d_Int.bin"%(nit,hook[9][name])
            else:
                filename = prefix+"_%d_%04d_Int.fmt"%(nit,hook[9][name])
            Compressor.writeUnsteadyCoefs(hook[2][name],hook[5][name], filename,"face",format)
#=========================================================================================
# Conversion of array of indices to files
# ---------------------------------------
# IN: iteration: entier designant la premiere iteration a laquelle le fichier doit etre lu
# IN: indices: une liste de numpy int32
# IN: fileName: le nom du fichier
# IN: format: bin_raw, fmt_raw
#=========================================================================================
def convertIndices2File(hook,fileName,iteration,format):
    for name in hook[3]:
        if hook[3][name] != []:
            convertIndices2File__(hook[0][name],hook[3][name],iteration, 'deltas_%s.fmt'%name, 'fmt_raw')

def convertIndices2File__(iteration,indices,nbTotalIterations,fileName,format):
    fileName=fileName+"_%d"%nbTotalIterations
    if format == 'bin_raw':
        tt = numpy.empty( (1,), dtype=Internal.E_NpyInt)
        f = open(fileName, "w") # open file for writing
        tt[0] = 1 # for endian test
        f.write(tt.tobytes())
        tt[0] = iteration # write first iteration of reading
        f.write(tt.tobytes())
        # write indices
        for i in indices:
            tt[0] = i.size
            f.write(tt.tobytes())
            f.write(i.tobytes())
        f.close()
    else:
        tt = numpy.empty( (1,), dtype=Internal.E_NpyInt)
        f = open(fileName, "w") # open file for writing

        tt[0] = iteration # write first iteration of reading
        tt.tofile(f, sep=" ", format='%d') ; f.write('\n')
        # write indices
        for i in indices:
            tt[0] = i.size
            tt.tofile(f, sep=" ", format='%d') ; f.write('\n')
            i.tofile(f, sep=" ", format='%d') ; f.write('\n')
        f.close()
