
import Converter.Mpi as Cmpi
from mpi4py import MPI
import Converter.Internal as Internal
import Converter.PyTree as C
import Converter.Distributed as CD

import Intersector.PyTree as XOR
from . import intersector

from collections import OrderedDict

import numpy

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
def adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None, com = MPI.COMM_WORLD):
  """Adapts an unstructured mesh a with respect to a sensor.
  Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""
  tp = Internal.copyRef(t)
  _adaptCells(tp, sensdata, sensor_type, smoothing_type, itermax, subdiv_type, hmesh, sensor, com)
  return tp

def _adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None, com = MPI.COMM_WORLD):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""

    if sensdata is None and sensor is None:
      print('INPUT ERROR : no source data to initialize a sensor')
      return

    if hmesh is None and sensor is not None:
      print ('INPUT ERROR : you must also give as an argument the hmesh targeted by the sensor')
      return

    NBZ = len(Internal.getZones(t))
    #print('NBZ : '+str(NBZ))

    procDict = Cmpi.getProcDict(t)
    zidDict = Cmpi.getPropertyDict(t, 'zid')
    if zidDict == {} :
      print ('error : zone id is not set on zones')
      return
    
    # FIXME : verifying both dicts have same keys ###############################################

    owesHmesh=0

    if hmesh is None :
      #print("create hm : ", NBZ)
      hmesh = createHMesh(t, subdiv_type)
      if hmesh == [None] : # no basic elts in t (assumed because hmesh creation as done [None]
        return
      owesHmesh=1

    owesSensor=0
    if sensor is None : 
      #print("create sensor")
      sensor = XOR.createSensor(hmesh, sensor_type, smoothing_type, itermax)
      owesSensor=1

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
    zone_to_zone_to_list_owned = getJoinsPtLists(t, zidDict, 'PointList')

    zonerank = getZonesRanks(zidDict, procDict)
    #print(zonerank)

    #print('adaptCells..')
    intersector.adaptCells_mpi(hmesh, sensor, zone_to_zone_to_list_owned, zonerank)#, com) #fxme agglo

    if owesHmesh == 1 : #and owesSensor == 1 :
    	#print("_conformizeHMesh")
    	_conformizeHMesh(t, hmesh, zidDict, procDict, zonerank, zone_to_zone_to_list_owned, com)

    if owesHmesh == 1 :
    # 	#print('delete owned hmesh')
     	XOR.deleteHMesh(hmesh)
    if owesSensor == 1 : 
    	#print('delete owned sensor')
    	XOR.deleteSensor(sensor)

def getZonesRanks(zidDict = None, procDict = None):
  zonerank = {}
  if zidDict == None : 
    zidDict = Cmpi.getPropertyDict(t, 'zid')
  if procDict == None:
    procDict = Cmpi.getProcDict(t)

  sorted_keys = sorted(zidDict, key=zidDict.get)
  #print(sorted_keys)
  for w in sorted_keys:
    zid = zidDict[w]
    rank = procDict[w]
    zonerank[zid]=rank
  #print(zonerank)
  return zonerank

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
      zid = CD.getProperty(z, 'zid')
      hm = intersector.createHMesh2(m, subdiv_type, zid)
      hmeshs.append(hm)

    return hmeshs


#==============================================================================
# getJoinsPtList : XXX
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook 
#==============================================================================
def getJoinsPtLists(t, zidDict, ptlname): # 'PointList' or 'PointListDonor'

  zone_to_zone_to_list_owned = {}

  zones = Internal.getZones(t)

  for z in zones:

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    nb_racs = len(raccords)

    zid = CD.getProperty(z, 'zid')

    zone_to_list_owned = {}

    for rac in raccords:
      rt = Internal.getNodeFromType1(rac, 'GridConnectivityType_t')
      
      donnorName = "".join(Internal.getValue(rac))
      #print(donnorName)
      jzid = zidDict[donnorName]
      #print('id is ' + str(jzid))
      ptList = Internal.getNodeFromName1(rac, ptlname)[1][0]
      #print ('PTLSZ : ' + str(len(ptList)))
      #print ('MMAXVAL : ' + str(max(ptList)))

      zone_to_list_owned[jzid] = ptList

    zone_to_zone_to_list_owned[zid] = zone_to_list_owned

  return zone_to_zone_to_list_owned


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
def _conformizeHMesh(t, hooks, zidDict = None, procDict = None, zonerank = None, zone_to_zone_to_list_owned = None, com = MPI.COMM_WORLD):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: _conformizeHMesh(t, hooks)"""
    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    if zidDict == None : 
      zidDict = Cmpi.getPropertyDict(t, 'zid')
    if procDict == None:
      procDict = Cmpi.getProcDict(t)
    if zonerank == None:
      zonerank = getZonesRanks(zidDict, procDict)
      #print(zonerank)

    ### 1. UPDATE BC & JOIN POINTLISTS WITH CURRENT ENABLING STATUS AND MPI EXCHANGES
    if zone_to_zone_to_list_owned == None :
      zone_to_zone_to_list_owned = getJoinsPtLists(t, zidDict, 'PointList')

    zone_to_bcptlists = getBCsPtLists(t)

    i=-1
    for z in zones:

        i +=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        if hooks[i] == None : continue

        zid = CD.getProperty(z, 'zid')

        # BC and Joins point list (owned side)
        bcptlists = zone_to_bcptlists[zid]
        jzone_to_ptlist = zone_to_zone_to_list_owned[zid]

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

        res = intersector.conformizeHMesh2(hooks[i], bcptlists, jzone_to_ptlist, fieldsC, fieldsN, fieldsF)

        # res[0] : mesh
        # res[1] : List : updated bc pointlist
        # res[2] : Dict : updated 'jzone to PtList'
        # res[3] : List : center fields
        # res[4] : List : node fields
        # res[5] : List : face fields

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
        jzone_to_ptlist = res[2]
        #print(jzone_to_ptlist)
        
        if jzone_to_ptlist != {}:
          updateJoinsPointLists3(z, zidDict, jzone_to_ptlist, 'PointList')


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

    _exchangePointLists(t, hooks, zidDict, procDict, zonerank, zone_to_bcptlists, zone_to_zone_to_list_owned, com)


#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------
def updateJoinsPointLists3(z, zidDict, jzone_to_ptlist, ptlname): # 'PointList', 'PointListDonor'

  joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
  #zname=z[0]
  #print(jzone_to_ptlist)

  # update the Join pointlist and synchronize with other zones (their PointListDonor)

  for j in joins:
    donnorName = "".join(Internal.getValue(j))
    donnor_zid = zidDict[donnorName]
    #print(donnorName + ' - ' + str(donnor_zid))
    ptl = Internal.getNodeFromName1(j, ptlname)
    L1 = jzone_to_ptlist[donnor_zid]
    L1 = numpy.reshape(L1, (1,len(L1)))
    ptl[1]= L1


def _exchangePointLists(t, hooks, zidDict=None, procDict=None, zonerank=None, zone_to_bcptlists=None, zone_to_zone_to_list_owned=None, com = MPI.COMM_WORLD):

  nb_hooks = len(hooks)
  zones = Internal.getZones(t)
  nb_zones = len(zones)

  if nb_zones != nb_hooks:
      print('must give one hook per zone')
      return

  if zidDict == None : 
    zidDict = Cmpi.getPropertyDict(t, 'zid')
  if procDict == None:
    procDict = Cmpi.getProcDict(t)
  if zonerank == None:
    zonerank = getZonesRanks(zidDict, procDict)
    #print(zonerank)

  if zone_to_zone_to_list_owned == None :
    zone_to_zone_to_list_owned = getJoinsPtLists(t, zidDict, 'PointList') #already up-to-date if conformizeHMesh is the caller or has been called

  # MPI exchange
  zone_to_zone_to_list_opp = intersector.exchangePointLists(zonerank, Cmpi.rank, Cmpi.size, zone_to_zone_to_list_owned)#, com)

  if zone_to_zone_to_list_opp == {} : return # single block
  #
  for z in zones:
    zid = CD.getProperty(z, 'zid')
    zone_to_list_opp = zone_to_zone_to_list_opp[zid]
    if zone_to_list_opp == {} : continue

    updateJoinsPointLists3(z, zidDict, zone_to_zone_to_list_opp[zid], 'PointListDonor')

