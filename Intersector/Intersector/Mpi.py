
import Converter.Mpi as Cmpi
from mpi4py import MPI
import Converter.Internal as Internal
import Converter.PyTree as C
import Converter.Distributed as CD

import Intersector.PyTree as XOR
from . import intersector

from collections import OrderedDict

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
  
    zonerank = {}
    sorted_keys = sorted(zidDict, key=zidDict.get)
    #print(sorted_keys)
    for w in sorted_keys:
      zid = zidDict[w]
      rank = procDict[w]
      zonerank[zid]=rank
    #print(zonerank)

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
    zone_to_zone_to_list_owned = getJoinsPtList(t, zidDict)

    #print('adaptCells..')
    intersector.adaptCells_mpi(hmesh, sensor, zone_to_zone_to_list_owned, zonerank)#, com) #fxme agglo

    if owesHmesh == 1 : #and owesSensor == 1 :
    	#print("_conformizeHMesh")
    	XOR._conformizeHMesh(t, hmesh)

    # if owesHmesh == 1 :
    # 	#print('delete owned hmesh')
    # 	XOR.deleteHMesh(hmesh)
    # if owesSensor == 1 : 
    # 	#print('delete owned sensor')
    # 	XOR.deleteSensor(sensor)


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
def getJoinsPtList(t, zidDict):

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
      ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]
      #print ('PTLSZ : ' + str(len(ptList)))
      #print ('MMAXVAL : ' + str(max(ptList)))

      zone_to_list_owned[jzid] = ptList

    zone_to_zone_to_list_owned[zid] = zone_to_list_owned

  return zone_to_zone_to_list_owned
