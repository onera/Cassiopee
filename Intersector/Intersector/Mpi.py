
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

def _adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, sensor_metric_policy = 0, subdiv_type=0, hmesh=None, sensor=None, com = MPI.COMM_WORLD, procDict=None, zidDict=None):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""

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

    if hmesh is None :
      #print("create hm : ", NBZ)
      hmesh = createHMesh(t, subdiv_type)
      if hmesh == [None] : # no basic elts in t (assumed because hmesh creation as done [None]
        return
      owesHmesh=1

    owesSensor=0
    if sensor is None : 
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
    zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)

    #
    zonerank = getZonesRanks(zidDict, procDict)
    
    rid_to_zones = getRidToZones(t, zidDict)

    #print('adaptCells..')
    intersector.adaptCells_mpi(hmesh, sensor, zone_to_rid_to_list_owned, zonerank, rid_to_zones)#, com) #fxme agglo
    
    if owesHmesh == 1 : #and owesSensor == 1 :
      #print("_conformizeHMesh")
      _conformizeHMesh(t, hmesh, zidDict, procDict, rid_to_zones, zonerank, zone_to_rid_to_list_owned, com)
    
    if owesHmesh == 1 :
    #   #print('delete owned hmesh')
      XOR.deleteHMesh(hmesh)
    if owesSensor == 1 : 
      #print('delete owned sensor')
      XOR.deleteSensor(sensor)
    

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
# setZonesAndJoinsUId : Assigns a zid and  the dictionary zid <-> rank
# IN: zidDict : dictionary zname <-> zid
# IN: procDict : dictionary zname <-> rank
# OUT: Returns the dictionary zid <-> rank
#==============================================================================
def setZonesAndJoinsUId(t):

  zs = Internal.getZones(t)
  iz = -1
  ir = -1

  # ---------------------------------------------------------------------
  # dict_VD = {z1 : {z2 : {0 : {[idx_min_1, idx_min_2], [ir_1, ir_2]}, 1 : {[idx_min_1, idx_min_2], [ir_3, ir_4]}}}
  #                 {z4 : {0 : {[idx_min], [ir_5]}, 1 : {[idx_min], [ir_6]}, 1 : {[idx_min], [ir_7]}, 1 : {[idx_min], [ir_8]}}}
  #           {z2 : {z3 : {0 : {[idx_min, ...], [ir, ...]}}}
  #           {z3 : {z4 : {0 : {[idx_min, ...], [ir, ...]}}}}
  # ---------------------------------------------------------------------

  dict_VD = {}
  zidDict = {}

  for z in zs: # loop on blocks
    iz += 1

    zidDict[z[0]] = iz #build zidDict

    l = Internal.newIntegralData(name='zid', parent=z) #ceate a zid for each block    
    l[1] = iz

  for z in zs: # loop on blocks

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t') # check if rac exists for block z
    nb_racs  = len(raccords)
    
    for rac in raccords: # loop on racs

      # GET ZONES INDICES
      z1 = CD.getProperty(z, 'zid')

      donnorName = "".join(Internal.getValue(rac))
      z2 = zidDict[donnorName]

      # re-arrange
      if z1 > z2:
        c  = z2
        z2 = z1
        z1 = c            

      ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]
      ptList_donor = Internal.getNodeFromName1(rac, 'PointListDonor')[1][0]

      is_perio = Internal.getNodeFromName1(rac, 'GridConnectivityProperty')
      idx_min  = numpy.minimum(numpy.min(ptList), numpy.min(ptList_donor))

      # 0==match connec; 1==perio connec
      if is_perio:
        idx_perio = 1
      else:
        idx_perio = 0

      # -----

      if z1 in dict_VD:
        if z2 in dict_VD[z1]:
          # if 'idx_min' in dict_VD[z1][z2]:
          if idx_perio in dict_VD[z1][z2]:
            if idx_min in dict_VD[z1][z2][idx_perio]:
              r    = Internal.newIntegralData(name='rid', parent=rac)
              r[1] = dict_VD[z1][z2][idx_perio][idx_min]
            else:
              ir += 1
              dict_VD[z1][z2][idx_perio][idx_min]=ir
              
              r    = Internal.newIntegralData(name='rid', parent=rac) 
              r[1] = ir
          else:
            ir += 1
            dict_VD[z1][z2][idx_perio] = {}
            dict_VD[z1][z2][idx_perio][idx_min]=ir

            r    = Internal.newIntegralData(name='rid', parent=rac) 
            r[1] = ir
        else:
          ir += 1
          dict_VD[z1][z2] = {}
          dict_VD[z1][z2][idx_perio] = {}
          dict_VD[z1][z2][idx_perio][idx_min]=ir

          r    = Internal.newIntegralData(name='rid', parent=rac) 
          r[1] = ir
      else:
        ir += 1
        dict_VD[z1] = {}
        dict_VD[z1][z2] = {}
        dict_VD[z1][z2][idx_perio] = {}
        dict_VD[z1][z2][idx_perio][idx_min]=ir

        r    = Internal.newIntegralData(name='rid', parent=rac) 
        r[1] = ir


#==============================================================================
# createHMesh : Returns a hierarchical mesh hook per zone
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook 
#==============================================================================
def createHMesh(t, subdiv_type= 0):
  """Returns a hierarchical mesh hook per zone.
  Usage: createHMesh(t, subdiv_type= 0)"""

  zones = Internal.getZones(t)
  
  hmeshs = []
  for z in zones:

    m = C.getFields(Internal.__GridCoordinates__, z)[0]

    # ----- 

    dico_rotation_to_ptList = {}

    key1 = [0.]*6
    key1 = tuple(key1)

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')

    for rac in raccords:
      rac_match_perio = Internal.getNodesFromType2(rac, 'GridConnectivityProperty_t')

      ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]

      if rac_match_perio: #perio BC
        center_rotat  = Internal.getNodesFromName3(rac, 'RotationCenter')[0][1]
        angle_rotat   = Internal.getNodesFromName3(rac, 'RotationAngle' )[0][1]

        # Condition to set negative angle as a translation like dict

        key2 = numpy.concatenate((center_rotat, angle_rotat))
        key2 = tuple(key2)

        if numpy.sum(angle_rotat) < 0: #negative rotation
          if key1 in dico_rotation_to_ptList:
            dico_rotation_to_ptList[key1] = numpy.concatenate((dico_rotation_to_ptList[key1], ptList))
          else:
            dico_rotation_to_ptList[key1] = ptList
        else: #translation or positive rotation
          if key2 in dico_rotation_to_ptList:
            dico_rotation_to_ptList[key2] = numpy.concatenate((dico_rotation_to_ptList[key2], ptList))
          else:
            dico_rotation_to_ptList[key2] = ptList
      else: #match BC
        if key1 in dico_rotation_to_ptList: #following match BC
          dico_rotation_to_ptList[key1] = numpy.concatenate((dico_rotation_to_ptList[key1], ptList))
        else: #first match BC
          dico_rotation_to_ptList[key1] = ptList

    # --- 
    
    m = intersector.initForAdaptCells(m, dico_rotation_to_ptList)

    zid = CD.getProperty(z, 'zid')
    hm = intersector.createHMesh2(m, subdiv_type, zid)
    hmeshs.append(hm)

  return hmeshs

#==============================================================================
# getRidToZones : Computes a dictionary rid <-> (zid1,zid2)
# IN: t : 3D NGON PyTree
# IN: zidDict : dictionary zname <-> zid
# OUT: Returns the dictionary rid <-> (zid1,zid2)
#==============================================================================
def getRidToZones(t, zidDict):
  """ Function returning ridDict (zones concerned 
  by rac rid, given tree t as well as zidDict.
  """

  ridDict = {}

  zones = Internal.getZones(t)
  for z in zones:
    zid = CD.getProperty(z, 'zid')

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for rac in raccords:
      rid = CD.getProperty(rac, 'rid')

      z1 = CD.getProperty(z, 'zid')

      donnorName = "".join(Internal.getValue(rac))
      z2 = zidDict[donnorName]

      ridDict[rid] = (z1,z2)

  return ridDict


#==============================================================================
# getJoinsPtList : XXX
# IN: t : 3D NGON PyTree
# IN: XXX
# OUT: XXX
#==============================================================================
def getJoinsPtLists(t, zidDict):

  zone_to_rid_to_list_owned = {}

  zones = Internal.getZones(t)

  for z in zones:
    zid = CD.getProperty(z, 'zid')

    rid_to_list_owned = {}

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for rac in raccords:
      
      rid = CD.getProperty(rac, 'rid')

      donnorName = "".join(Internal.getValue(rac))
      ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]

      jzid = zidDict[donnorName]
      if jzid==zid: # auto-match case
        if rid in rid_to_list_owned:
          rid_to_list_owned[rid] = numpy.concatenate((rid_to_list_owned[rid], ptList))
        else:
          rid_to_list_owned[rid] = ptList
      else:         # other cases
        rid_to_list_owned[rid] = ptList

    zone_to_rid_to_list_owned[zid] = rid_to_list_owned

  return zone_to_rid_to_list_owned


#==============================================================================
# getBCsPtLists : XXX
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook 
#==============================================================================
def getBCsPtLists(t):

  zone_to_bcptlists = {}
  zones = Internal.getZones(t)
  #
  for z in zones:
    zid = CD.getProperty(z, 'zid')
    # BC and Joins point list (owned side)
    bcptlists = XOR.getBCPtList(z)
    zone_to_bcptlists[zid]=bcptlists

  return zone_to_bcptlists

#==============================================================================
# _conformizeHMesh : Converts the basic element leaves of a hierarchical mesh (hooks is a list of hooks to hiearchical zones) to a conformal polyhedral mesh.
#                   Each hiearchcial zone is referring to a zone in the original mesh t. So the mesh is replaced in the tree and the BCs/Joins/Fields are transferred.
# IN: t : PyTree before adaptation
# IN: hook : list of hooks to hiearchical zones (same size as nb of zones in t).
# OUT: Nothing 
#==============================================================================
def _conformizeHMesh(t, hooks, zidDict, procDict, rid_to_zones = None, zonerank = None, zone_to_rid_to_list_owned = None, com = MPI.COMM_WORLD):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: _conformizeHMesh(t, hooks)"""

    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    #tt = time.time()
    
    if zonerank == None:
      zonerank = getZonesRanks(zidDict, procDict)
      #print(zonerank)

    # if Cmpi.rank == 0:
    #   print('_conformizeHMesh CPU get Dicts : ' + str(time.time() -tt))

    #tt = time.time()

    ### 1. UPDATE BC & JOIN POINTLISTS WITH CURRENT ENABLING STATUS AND MPI EXCHANGES
    if zone_to_rid_to_list_owned == None :
      zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)

    # if Cmpi.rank == 0:
    #   print('_conformizeHMesh CPU get Joins PTLists : ' + str(time.time() -tt))

    #tt = time.time()

    zone_to_bcptlists = getBCsPtLists(t)

    # if Cmpi.rank == 0:
    #   print('_conformizeHMesh CPU get BC PTLists : ' + str(time.time() -tt))

    if rid_to_zones == None:
      rid_to_zones = getRidToZones(t, zidDict)

    #dtconf=0
    #dtfield = 0

    i=-1
    for z in zones:

      i +=1
      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      if m == []: continue
      if hooks[i] == None : continue

      zid = CD.getProperty(z, 'zid')

      # BC and Joins point list (owned side)
      bcptlists = zone_to_bcptlists[zid]
      rid_to_ptlist = zone_to_rid_to_list_owned[zid]

      fieldsC = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
      if fieldsC == [] : fieldsC = None

      fieldsN = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
      if fieldsN == [] : fieldsN = None

      try:
        pfieldsF = Internal.getNodeFromName(z, 'CADData')
        fieldsF = Internal.getChildFromName(pfieldsF, 'fcadid')
        fieldsF = [fieldsF[1]]
        if fieldsF == [] : fieldsF = None # Unnecessary ?
        
      except TypeError:
        fieldsF = None

      #tconf=time.time()
      res = intersector.conformizeHMesh2(hooks[i], bcptlists, rid_to_ptlist, rid_to_zones, fieldsC, fieldsN, fieldsF)
      #dtconf += time.time() - tconf

      # res[0] : mesh
      # res[1] : List : updated bc pointlist
      # res[2] : Dict : updated 'jzone to PtList'
      # res[3] : List : center fields
      # res[4] : List : node fields
      # res[5] : List : face fields

      #tfield = time.time()

      mesh = res[0]

      # MAJ du maillage de la zone
      C.setFields([mesh], z, 'nodes')

      if len(res) < 2 : continue

      # MAJ BCs
      bcptlists = res[1]
      #print(bcptlists)

      if bcptlists != [] :
        XOR.updateBCPointLists2(z, bcptlists)
      else:
        C._deleteZoneBC__(z)

      # MAJ Joins
      rid_to_ptlist = res[2]

      zone_to_rid_to_list_owned[zid] = rid_to_ptlist #dico z_to_rid_to_ptlist_owned est mis à jour avant de le donner à exchangePointList

      if rid_to_ptlist != {}:
        updateJoinsPointLists3(z, zidDict, rid_to_ptlist, 'PointList')

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
      if fieldz != [None] :
        Internal.newDataArray('fcadid', value=fieldz[0], parent=Internal.getNodeFromName(z, 'CADData'))

      #dtfield += time.time() - tfield

    #texch = time.time()
    _exchangePointLists(t, hooks, zidDict, procDict, rid_to_zones, zonerank, zone_to_bcptlists, zone_to_rid_to_list_owned, com)
    #dtexch = time.time() - texch

    #if Cmpi.rank == 0:
    #print('_conformizeHMesh CPU conformization : ' + str(dtconf))
    #print('_conformizeHMesh CPU MAJ Fields : ' + str(dtfield))
    #print('_conformizeHMesh CPU MPI Exch PtLists : ' + str(dtexch))


#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------
def updateJoinsPointLists3(z, zidDict, rid_to_ptlist, ptlname): # 'PointList', 'PointListDonor'

  joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
  zid = CD.getProperty(z, 'zid')

  # update the Join pointlist and synchronize with other zones (their PointListDonor)

  processed_rid = set()

  for j in joins:

    ptl    = Internal.getNodeFromName1(j, ptlname)

    rid    = CD.getProperty(j, 'rid')
    donnorName = "".join(Internal.getValue(j))
    jzid   = zidDict[donnorName]

    if rid not in rid_to_ptlist : continue

    L1     = rid_to_ptlist[rid]

    if jzid==zid:
      half   = int(len(rid_to_ptlist[rid]) / 2)

      if ptlname == 'PointList':        # appel apres conformisation
        if rid in processed_rid:        ## deuxieme passe (i.e deuxime demi-raccord)
          L1 = L1[half:]
        else:                           ## premiere passe
          processed_rid.add(rid)
          L1 = L1[:half]
      elif ptlname == 'PointListDonor': # appel apres echange de pointlist
        if rid in processed_rid:        ## deuxieme passe (i.e deuxime demi-raccord)
          L1 = L1[:half]
        else:                           ## premiere passe
          processed_rid.add(rid)
          L1 = L1[half:]
  
    L1     = numpy.reshape(L1, (1,len(L1)))
    ptl[1] = L1


def _exchangePointLists(t, hooks, zidDict, procDict, rid_to_zones = None, zonerank=None, zone_to_bcptlists=None, zone_to_rid_to_list_owned=None, com = MPI.COMM_WORLD):

  nb_hooks = len(hooks)
  zones = Internal.getZones(t)
  nb_zones = len(zones)

  if nb_zones != nb_hooks:
    print('must give one hook per zone')
    return

  if zonerank == None:
    zonerank = getZonesRanks(zidDict, procDict)
    #print(zonerank)

  if zone_to_rid_to_list_owned == None :
    zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict) #already up-to-date if conformizeHMesh is the caller or has been called

  if rid_to_zones == None:
    rid_to_zones = getRidToZones(t, zidDict)

  # MPI exchange
  # zone_to_rid_to_list_opp or zone_to_zone_to_list_opp
  zone_to_rid_to_list_opp = intersector.exchangePointLists(rid_to_zones, zonerank, Cmpi.rank, Cmpi.size, zone_to_rid_to_list_owned)#, com)

  if zone_to_rid_to_list_opp == {} : return # single block
  #
  for z in zones:
    zid = CD.getProperty(z, 'zid')
    if zid not in zone_to_rid_to_list_opp : continue

    rid_to_list_opp = zone_to_rid_to_list_opp[zid]
    if rid_to_list_opp == {} : continue

    updateJoinsPointLists3(z, zidDict, rid_to_list_opp, 'PointListDonor')

