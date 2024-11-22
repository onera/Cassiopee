
import Converter.Mpi as Cmpi
from mpi4py import MPI
import Converter.Internal as Internal
import Converter.PyTree as C
import Converter.Distributed as CD

import Intersector.PyTree as XOR
from . import intersector

from collections import OrderedDict

import numpy
#import time as time

#==============================================================================
# getZonesRanks : Computes the dictionary zid <-> rank
# IN: zidDict : dictionary zname <-> zid
# IN: procDict : dictionary zname <-> rank
# OUT: Returns the dictionary zid <-> rank
#==============================================================================
def getZonesRanks(zidDict, procDict):
    zonerank = {}
    sorted_keys = sorted(zidDict, key=zidDict.get)
    #print(sorted_keys)
    for w in sorted_keys:
        zid = zidDict[w]
        rank = procDict[w]
        zonerank[zid]=rank
    #print(zonerank)
    return zonerank

#==============================================================================
# adaptCells : Adapts an unstructured mesh a with respect to a sensor
# IN: t : 3D NGON mesh
# IN: sensdata : sensor data (a bunch of vertices or a mesh for a geom sensor, a mesh for a xsensor, punctual values for a nodal or cell sensor)
# IN: sensor_type : geom_sensor (0) , xsensor (1), nodal_sensor (2), cell_sensor(3), xsensor2(4)
# IN smoothing_type : First-neighborhood (0) Shell-neighborhood(1)
# IN itermax : max number of level in the hierarchy
# IN: subdiv_type : isotropic currently
# IN: hmesh : hierarchical mesh hook
# IN: sensor : sensor hook
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, sensor_metric_policy = 0, subdiv_type=0, hmesh=None, sensor=None, com = MPI.COMM_WORLD, procDict=None, zidDict=None):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""
    tp = Internal.copyRef(t)
    _adaptCells(tp, sensdata, sensor_type, smoothing_type, itermax, sensor_metric_policy, subdiv_type, hmesh, sensor, com, procDict, zidDict)
    return tp

def _adaptCells(t, sensdata=None, sensor_type=0, smoothing_type=0, itermax=-1, sensor_metric_policy=0, subdiv_type=0, hmesh=None, sensor=None, com=MPI.COMM_WORLD, procDict=None, zidDict=None):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""

    if sensor_type == 1: sensor_type=4 # initial xsensor does not exist anymore

    if sensdata is None and sensor is None:
        print('INPUT ERROR : no source data to initialize a sensor')
        return

    if hmesh is None and sensor is not None:
        print ('INPUT ERROR : you must also give as an argument the hmesh targeted by the sensor')
        return

    if procDict is None or procDict == {}:
        print ('INPUT ERROR : you must also give as an argument the processors dictionary')
        return

    if zidDict is None or zidDict == {}:
        print ('INPUT ERROR : you must also give as an argument the zone id dictionary')
        return

    #tt0 = time.time()

    # FIXME : verifying both dicts have same keys ###############################################

    NBZ = len(Internal.getZones(t))
    #print('NBZ : '+str(NBZ))

    owesHmesh=0

    #tt0 = time.time()

    if hmesh is None:
        #print("create hm : ", NBZ)
        hmesh = XOR.createHMesh(t, subdiv_type)
        if hmesh == [None] : # no basic elts in t (assumed because hmesh creation as done [None]
            return
        owesHmesh=1

    owesSensor=0
    if sensor is None: 
        #print("create sensor")
        sensor = XOR.createSensor(hmesh, sensor_type, smoothing_type, itermax, sensor_metric_policy)
        owesSensor=1

    #if Cmpi.rank == 0 :
    #print('rank :' + str(Cmpi.rank)+' python : create HMesh and Sensor ' + str(time.time()-tt0))

    err=0
    if sensdata is not None:
        #print("assignData2Sensor")
        if sensor_type == 4:
            sensdata = C.convertArray2NGon(sensdata)
        err = XOR.assignData2Sensor(sensor, sensdata)
        if err == 1:
            print('INPUT ERROR : sensor data list must be sized as nb of sensors')
            return

    #
    zone_to_rid_to_list_owned = XOR.getJoinsPtLists(t, zidDict)

    #
    zonerank = getZonesRanks(zidDict, procDict)

    rid_to_zones = XOR.getRidToZones(t, zidDict)

    #print('adaptCells..')
    intersector.adaptCells_mpi(hmesh, sensor, zone_to_rid_to_list_owned, zonerank, rid_to_zones, com) #fxme agglo

    if owesHmesh == 1 : #and owesSensor == 1 :
        #print("_conformizeHMesh")
        _conformizeHMesh(t, hmesh, zidDict, procDict, rid_to_zones, zonerank, zone_to_rid_to_list_owned, com)

    if owesHmesh == 1:
    #   #print('delete owned hmesh')
        XOR.deleteHMesh(hmesh)
    if owesSensor == 1: 
        #print('delete owned sensor')
        XOR.deleteSensor(sensor)

#==============================================================================
# _conformizeHMesh : Converts the basic element leaves of a hierarchical mesh (hooks is a list of hooks to hiearchical zones) to a conformal polyhedral mesh.
#                   Each hiearchcial zone is referring to a zone in the original mesh t. So the mesh is replaced in the tree and the BCs/Joins/Fields are transferred.
# IN: t : PyTree before adaptation
# IN: hook : list of hooks to hiearchical zones (same size as nb of zones in t).
# OUT: Nothing 
#==============================================================================
def _conformizeHMesh(t, hooks, zidDict, procDict, rid_to_zones=None, zonerank=None, zone_to_rid_to_list_owned=None, com=MPI.COMM_WORLD):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: _conformizeHMesh(t, hooks)"""

    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    if zonerank is None:
        zonerank = getZonesRanks(zidDict, procDict)
        #print(zonerank)

    ### 1. UPDATE BC & JOIN POINTLISTS WITH CURRENT ENABLING STATUS AND MPI EXCHANGES
    if zone_to_rid_to_list_owned is None:
        zone_to_rid_to_list_owned = XOR.getJoinsPtLists(t, zidDict)

    zone_to_bcptlists = XOR.getBCsPtLists(t)

    if rid_to_zones is None:
        rid_to_zones = XOR.getRidToZones(t, zidDict)

    i = -1
    for z in zones:

        i +=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        if hooks[i] is None: continue

        zid = CD.getProperty(z, 'zid')

        # BC and Joins point list (owned side)
        bcptlists = zone_to_bcptlists[zid]
        rid_to_ptlist = zone_to_rid_to_list_owned[zid]

        fieldsC = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
        if fieldsC == []: fieldsC = None

        fieldsN = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fieldsN == []: fieldsN = None

        try:
            pfieldsF = Internal.getNodeFromName(z, 'CADData')
            fieldsF = Internal.getChildFromName(pfieldsF, 'fcadid')
            fieldsF = [fieldsF[1]]
            if fieldsF == []: fieldsF = None # Unnecessary ?

        except TypeError:
            fieldsF = None

        #tconf=time.time()
        #for now, always conformize the output
        conformize = 1
        res = intersector.conformizeHMesh(hooks[i], bcptlists, rid_to_ptlist, rid_to_zones, fieldsC, fieldsN, fieldsF, conformize)
        #dtconf += time.time() - tconf

        # res[0] : mesh
        # res[1] : List : updated bc pointlist
        # res[2] : Dict : updated 'jzone to PtList'
        # res[3] : List : center fields
        # res[4] : List : node fields
        # res[5] : List : face fields

        mesh = res[0]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        if len(res) < 2: continue

        # MAJ BCs
        bcptlists = res[1]
        #print(bcptlists)

        if bcptlists != []:
            XOR.updateBCPointLists2(z, bcptlists)
        else:
            C._deleteZoneBC__(z)

        # MAJ Joins
        rid_to_ptlist = res[2]

        zone_to_rid_to_list_owned[zid] = rid_to_ptlist #dico z_to_rid_to_ptlist_owned est mis à jour avant de le donner à exchangePointList

        if rid_to_ptlist != {}:
            XOR.updateJoinsPointLists3(z, zidDict, rid_to_ptlist, 'PointList')

        # MAJ FIELDS

        C._deleteFlowSolutions__(z)

        ## center fields
        fieldz = [res[3]]
        #print (fieldz)
        if fieldz != [None]:
            for f in fieldz:
                C.setFields([f], z, 'centers', False)

        # # ## node fields
        fieldz = [res[4]]
        #print (fieldz)
        if fieldz != [None]:
            for f in fieldz:
                C.setFields([f], z, 'nodes', False)

        ## face fields 
        fieldz = [res[5]]
        #print(fieldz)
        if fieldz != [None]:
            Internal.newDataArray('fcadid', value=fieldz[0], parent=Internal.getNodeFromName(z, 'CADData'))

    _exchangePointLists(t, hooks, zidDict, procDict, rid_to_zones, zonerank, zone_to_bcptlists, zone_to_rid_to_list_owned, com)



#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------
def _exchangePointLists(t, hooks, zidDict, procDict, rid_to_zones = None, zonerank=None, zone_to_bcptlists=None, zone_to_rid_to_list_owned=None, com = MPI.COMM_WORLD):

    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    if zonerank is None:
        zonerank = getZonesRanks(zidDict, procDict)
        #print(zonerank)

    if zone_to_rid_to_list_owned == None :
        zone_to_rid_to_list_owned = XOR.getJoinsPtLists(t, zidDict) #already up-to-date if conformizeHMesh is the caller or has been called

    if rid_to_zones is None:
        rid_to_zones = XOR.getRidToZones(t, zidDict)

    # MPI exchange
    # zone_to_rid_to_list_opp or zone_to_zone_to_list_opp
    zone_to_rid_to_list_opp = intersector.exchangePointLists(rid_to_zones, zonerank, Cmpi.rank, Cmpi.size, zone_to_rid_to_list_owned, com)

    if zone_to_rid_to_list_opp == {}: return # single block
    #
    for z in zones:
        zid = CD.getProperty(z, 'zid')
        if zid not in zone_to_rid_to_list_opp: continue

        rid_to_list_opp = zone_to_rid_to_list_opp[zid]
        if rid_to_list_opp == {}: continue

        XOR.updateJoinsPointLists3(z, zidDict, rid_to_list_opp, 'PointListDonor')


#------------------------------------------------------------------------------
# closeCells MPI
#------------------------------------------------------------------------------
def closeCells(t, procDict, zidDict, com = MPI.COMM_WORLD):
    """Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
    Usage: closeCells(t, com, procDict, zidDict)"""
    tp = Internal.copyRef(t)
    _closeCells(tp, procDict, zidDict, com)
    return tp

def _closeCells(t, procDict, zidDict, com = MPI.COMM_WORLD):
    """Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
    Usage: closeCells(t, com, procDict, zidDict)"""
    if procDict is {}:
        print ('INPUT ERROR: you must also give as an argument the processors dictionary')
        return

    if zidDict is {}:
        print ('INPUT ERROR: you must also give as an argument the zone id dictionary')
        return

    #
    zid_to_rid_to_list_owned = XOR.getJoinsPtLists(t, zidDict)
    #
    zonerank = getZonesRanks(zidDict, procDict)

    rid_to_zones = XOR.getRidToZones(t, zidDict)

    zs = Internal.getZones(t)

    meshes = []
    zids = []
    for z in zs:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        meshes.append(m)
        zid = CD.getProperty(z, 'zid')
        zids.append(zid)

    meshes = intersector.closeCells_mpi(meshes, zids, zid_to_rid_to_list_owned, zonerank, rid_to_zones, com) #fxme agglo

    # associate zid and closed mesh
    zid_to_m = {}
    for i in range(len(zids)):
        zid = zids[i]
        m = meshes[i]
        zid_to_m[zid] = m

    for z in zs : 
        zid = CD.getProperty(z, 'zid')
        m = zid_to_m[zid]
        # MAJ du maillage de la zone
        C.setFields([m], z, 'nodes')

    return t
