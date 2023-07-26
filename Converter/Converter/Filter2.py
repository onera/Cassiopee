# Part loader - taken from old maia
import Converter.PyTree as C
import Converter.Internal as I
import Converter.Filter as Filter
import Converter.Mpi as Cmpi
import numpy
from functools import partial

#============================================================
# add size to tree
#============================================================
def _addSizesToZoneTree(zone, zone_path, size_data):
  """
  Creates the MyArray#Size node using the size_data dict on the given zone
  for the following nodes:
  - ElementConnectivity array of Element_t nodes
  - PointList (or Unstr PointRange) array of BC_t
  - PointList array of GC_t, GC1to1_t, BCDataSet_t and ZoneSubRegion_t nodes
  - PointListDonor array of GC_t and GC1to1_t nodes
  """
  for elmt in I.getNodesFromType1(zone, 'Elements_t'):
    elmt_path = zone_path+"/"+elmt[0]
    ec_path   = elmt_path+"/ElementConnectivity"
    I.newIndexArray('ElementConnectivity#Size', value=size_data[ec_path][2], parent=elmt)

  for zone_bc in I.getNodesFromType1(zone, 'ZoneBC_t'):
    zone_bc_path = zone_path+"/"+zone_bc[0]
    for bc in I.getNodesFromType1(zone_bc, 'BC_t'):
      bc_path = zone_bc_path+"/"+bc[0]
      if I.getNodeFromName1(bc, 'PointList') is not None:
        pl_path = bc_path+"/PointList"
        I.newIndexArray('PointList#Size', value=size_data[pl_path][2], parent=bc)
      for bcds in I.getNodesFromType1(bc, 'BCDataSet_t'):
        if I.getNodeFromName1(bcds, 'PointList') is not None:
          pl_path = bc_path+"/"+bcds[0]+"/PointList"
          I.newIndexArray('PointList#Size', value=size_data[pl_path][2], parent=bcds)

  for zone_gc in I.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
    zone_gc_path = zone_path+"/"+zone_gc[0]
    gcs = I.getNodesFromType1(zone_gc, 'GridConnectivity_t') \
        + I.getNodesFromType1(zone_gc, 'GridConnectivity1to1_t')
    for gc in gcs:
      gc_path = zone_gc_path+"/"+gc[0]
      if I.getNodeFromName1(gc, 'PointList') is not None:
        pl_path = gc_path+"/PointList"
        I.newIndexArray('PointList#Size', value=size_data[pl_path][2], parent=gc)
      if I.getNodeFromName1(gc, 'PointListDonor') is not None:
        pld_path = gc_path+"/PointListDonor"
        I.newIndexArray('PointListDonor#Size', value=size_data[pld_path][2], parent=gc)

  for zone_subregion in I.getNodesFromType1(zone, 'ZoneSubRegion_t'):
    zone_subregion_path = zone_path+"/"+zone_subregion[0]
    if I.getNodeFromName1(zone_subregion, 'PointList') is not None:
      pl_path = zone_subregion_path+"/PointList"
      I.newIndexArray('PointList#Size', value=size_data[pl_path][2], parent=zone_subregion)

  for flow_sol in I.getNodesFromType1(zone, 'FlowSolution_t'):
    sol_path = zone_path + "/" + I.getName(flow_sol)
    if I.getNodeFromName1(flow_sol, 'PointList') is not None:
      pl_path = sol_path+"/PointList"
      I.newIndexArray('PointList#Size', value=size_data[pl_path][2], parent=flow_sol)

def _addSizesToTree(size_tree, size_data):
  """
  Convenience function which loops over zones to add size
  data in each one.
  """
  for base in I.getNodesFromType1(size_tree, 'CGNSBase_t'):
    base_path = '/'+base[0]
    for zone in I.getZones(base):
      zone_path = base_path+"/"+zone[0]
      _addSizesToZoneTree(zone, zone_path, size_data)
  return None

#=====================================================================
# load un squelette sur le rank 0 et bcast
# ajoute les noeuds size
#=====================================================================
def loadCollectiveSizeTree(filename):
  """
    Load on all ranks a "size tree"
    a size tree is a partial tree that contains only the data needed to distribute the tree:
      nb of nodes, nb of elements, size of bcs and gcs...
    Convention:
      when we load the dimensions of an array "MyArray" without loading the array,
      then the dimensions are kept in a "MyArray#Size" node,
      at the same level as the array node would be.
  """
  skeleton_depth  = 7
  skeleton_n_data = 3

  # In order to avoid filesystem overload only 1 proc reads the squeleton, then we broadcast
  if Cmpi.rank == 0:
    size_data = {}
    t = C.convertFile2PyTree(filename,
                             skeletonData=[skeleton_n_data, skeleton_depth],
                             dataShape=size_data,
                             format='bin_hdf')
    _addSizesToTree(t, size_data)
    #fix_point_ranges(t)
    #load_grid_connectivity_property(filename, t)
  else: t = None
  t = Cmpi.bcast(t, root=0)
  return t

#============================================================
# distribution
#============================================================
def newDistribution(distributions=dict(), parent=None):
  """
  Create and return a CGNSNode to be used to store distribution data
  Attach it to parent node if not None
  In addition, add distribution arrays specified in distributions dictionnary.
  distributions must be a dictionnary {DistriName : distri_array}
  """
  distriNode = I.newUserDefinedData(':CGNS#Distribution', None, parent)
  for name, value in distributions.items():
    I.newDataArray(name, value, parent=distriNode)
  return distriNode

def _createDistributionNodeFromDistrib(name, parent_node, distrib):
  distrib_ud = newDistribution(parent=parent_node)
  I.newDataArray(name, value=distrib, parent=distrib_ud)
  return None

def getDistribution(node, distri_name=None):
  return I.getNodeFromPath(node, '/'.join([':CGNS#Distribution', distri_name])) if distri_name \
      else I.getNodeFromName1(node, ':CGNS#Distribution')

def uniformDistributionAt(n_elt, i, n_interval):
  step      = n_elt // n_interval
  remainder = n_elt %  n_interval

  if i < remainder:
    inf = i * (step + 1)
    sup = inf + step + 1
  else:
    inf = i * step + remainder
    sup = inf + step

  return inf, sup

def uniformDistribution(n_elt, comm):
  int_type = type(n_elt)
  i_rank = int_type(comm.Get_rank())
  n_rank = int_type(comm.Get_size())
  u_dist = uniformDistributionAt(n_elt, i_rank, n_rank)
  proc_indices = numpy.empty(3, dtype=type(n_elt))
  proc_indices[0] = u_dist[0]
  proc_indices[1] = u_dist[1]
  proc_indices[2] = n_elt
  return proc_indices

def createDistributionNodeFromDistrib(name, parent_node, distrib):
  distrib_ud = newDistribution(parent=parent_node)
  I.newDataArray(name, value=distrib, parent=distrib_ud)

def computeElementsDistribution(zone, comm):
  if I.getZoneType(zone) == 1: pass
  else:
    elts = I.getNodesFromType1(zone, 'Elements_t')
    for elt in elts:
      er = I.getNodeFromName(elt, 'ElementRange')
      n_tot_elmt = er[1][1] - er[1][0] + 1
      createDistributionNode(n_tot_elmt, comm, 'Element', elt)

def createDistributionNode(n_elt, comm, name, parent_node):
  distrib = uniformDistribution(n_elt, comm)
  createDistributionNodeFromDistrib(name, parent_node, distrib)

def _computeZoneDistribution(zone, comm):
  n_vtx = C.getNPts(zone)
  n_cell = C.getNCells(zone)

  distrib_vtx  = createDistributionNode(n_vtx , comm, 'Vertex', zone)
  distrib_cell = createDistributionNode(n_cell, comm, 'Cell'  , zone)

  computeElementsDistribution(zone, comm)

  #for zone_subregion in I.getNodesFromType1(zone, 'ZoneSubRegion_t'):
  #  compute_plist_or_prange_distribution(zone_subregion, comm)

  #for flow_sol in I.getNodesFromType1(zone, 'FlowSolution_t'):
  #  compute_plist_or_prange_distribution(flow_sol, comm)

  #for zone_bc in I.getNodesFromType1(zone, 'ZoneBC_t'):
  #  for bc in I.getNodesFromType1(zone_bc, 'BC_t'):
  #    compute_plist_or_prange_distribution(bc, comm)
  #    for bcds in I.getNodesFromType1(bc, 'BCDataSet_t'):
  #      compute_plist_or_prange_distribution(bcds, comm)

  #for zone_gc in I.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
  #  gcs = I.getNodesFromType1(zone_gc, 'GridConnectivity_t') + I.getNodesFromType1(zone_gc, 'GridConnectivity1to1_t')
  #  for gc in gcs:
  #    compute_plist_or_prange_distribution(gc, comm)

def _addDistributionInfo(t, comm=Cmpi.KCOMM):
  for z in I.getZones(t):
    _computeZoneDistribution(z, comm)
  return None

#=============================================================
# Filters
#=============================================================
def plOrPrSize(node):
  """
  For a given node, search for a PointList of a PointRange children and return the
  size of the corresponding region:
   - directly from the value of PointRange for the PointRange case
   - using PointList#Size, if existing, or the :CGNS#Distribution otherwise for the
     PointList case.
  If both node are present, size of PointList is returned (but this should not happens!)
  Raises if both are absent.
  """
  # CBX
  index = I.getNodesFromType(node, 'IndexArray_t')
  for ind in index:
    ps = I.getNodeFromName1(ind, 'PointList#Size')
    if ps is not None: return ps
    #No PL#Size, try to get info from :CGNS#Distribution and suppose size = 1,N
    distri = I.getVal(getDistribution(node, 'Index'))
    return numpy.array([1, distri[2]])
  index = I.getNodesFromType(ind, 'IndexRange_t')
  for ind in index:
    return ind
    
def apply_dataspace_to_arrays(node, node_path, data_space, hdf_filter):
  """
  Fill the hdf_filter with the specified data_space for all the DataArray_t nodes
  below the parent node node
  """
  for data_array in I.getNodesFromType1(node, 'DataArray_t'):
    path = node_path+"/"+data_array[0]
    hdf_filter[path] = data_space

def apply_dataspace_to_pointlist(node, node_path, data_space, hdf_filter):
  """
  Fill the hdf_filter with the specified data_space for PointList and PointListDonor nodes
  (if existing) below the parent node node
  """
  if I.getNodeFromName1(node, 'PointList') is not None:
    hdf_filter[node_path + "/PointList"] = data_space
  if I.getNodeFromName1(node, 'PointListDonor') is not None:
    hdf_filter[node_path + "/PointListDonor"] = data_space

def cell_to_indexes(i_cell, plan_size, line_size):
  """ Compute the (i,j,k) indices of a cell or a node
  from its global index.
  Numbering convention is increasing i,j,k. Here global index
  and i,j,k start at 0.
  """
  k = i_cell // plan_size
  j = (i_cell - k*plan_size) // line_size
  i = i_cell - k*plan_size - j*line_size
  return i,j,k

def compute_slabs(array_shape, gnum_interval):
  """ Compute HDF HyperSlabs to be used in order to contiguously load a part
  of a structured tridimensionnal array.

  array_shape   : Number of elements in x,y,z directions
  gnum_interval : semi open interval of elements to load in global numbering
  returns    : list of hyperslabs

  Each slab is a list [[istart, iend], [jstart, jend], [kstart, kend]] of
  semi open intervals, starting at zero. The flattening order for the 3d
  array is increasing i, j, k.
  """
  hslab_list = []
  nx, ny, nz = array_shape
  line_size = nx
  plan_size = nx*ny

  ncell_to_load = gnum_interval[1] - gnum_interval[0]
  # print("{0} : cellInterval is [{1}:{2}[\n".format(iRank, gnum_interval[0], gnum_interval[1]))
  imin, jmin, kmin = cell_to_indexes(gnum_interval[0],   plan_size, line_size)
  imax, jmax, kmax = cell_to_indexes(gnum_interval[1]-1, plan_size, line_size)

  # print('toLoad : {0}  -- {1} {2} {3}  -- {4} {5} {6} \n'.format(
  #  ncell_to_load, imin, jmin, kmin, imax, jmax, kmax))

  istart = imin
  jstart = jmin
  kstart = kmin

  #If the line is full, merged it with next plan
  this_line_size = min(nx, istart+ncell_to_load) - istart
  if this_line_size != nx:
    jstart += 1
    if this_line_size > 0:
      start_line  = [[istart, min(nx, istart+ncell_to_load)], [jmin, jmin+1], [kmin, kmin+1]]
      ncell_to_load -= this_line_size
      hslab_list.append(start_line)
      # print('start_line {0}, loaded {1} elmts\n'.format(
        # start_line, start_line[0][1] - start_line[0][0]))

  #If the plan is full, merged it with the block
  this_plan_size = min(ny, jstart+(ncell_to_load // nx)) - jstart
  if this_plan_size != ny:
    kstart += 1
    if this_plan_size > 0:
      start_plane = [[0, nx], [jstart, min(ny, jstart+(ncell_to_load // nx))], [kmin, kmin+1]]
      ncell_to_load -= nx*this_plan_size
      hslab_list.append(start_plane)
      # print('start_plane {0}, loaded {1} lines ({2} elmts)\n'.format(
        # start_plane, start_plane[1][1] - start_plane[1][0], nx*(start_plane[1][1] - start_plane[1][0])))

  this_block_size = min(nz, kstart+(ncell_to_load // plan_size)) - kstart
  if this_block_size > 0:
    central_block = [[0, nx], [0, ny], [kstart, min(nz, kstart+(ncell_to_load // plan_size))]]
    ncell_to_load -= plan_size*this_block_size
    hslab_list.append(central_block)
    # print('central_block {0}, loaded {1} planes ({2} elmts)\n'.format(
      # central_block, central_block[2][1] - central_block[2][0], plan_size*(central_block[2][1] - central_block[2][0])))

  if ncell_to_load >= nx:
    end_plane = [[0, nx], [0, (ncell_to_load // nx)], [kmax, kmax+1]]
    ncell_to_load -= nx*(end_plane[1][1] - end_plane[1][0])
    hslab_list.append(end_plane)
    # print('end_plane {0}, loaded {1} lines ({2} elmts)\n'.format(
      # end_plane, end_plane[1][1] - end_plane[1][0], nx*(end_plane[1][1] - end_plane[1][0])))
  if ncell_to_load > 0:
    end_line = [[0, ncell_to_load], [jmax, jmax+1], [kmax, kmax+1]]
    ncell_to_load -= (end_line[0][1] - end_line[0][0])
    hslab_list.append(end_line)
    # print('end_line {0}, loaded {1} elmts\n'.format(
      # end_line, end_line[0][1] - end_line[0][0]))
  assert(ncell_to_load == 0)

  return hslab_list

def create_combined_dataspace(data_shape, distrib):
  """
  Create a dataspace from a flat distribution, but for arrays having a 3d (resp. 2d) stucture
  ie (Nx, Ny, Nz) (resp. (Nx, Ny)) numpy arrays.
  First, the 1d distribution is converted into slabs to load with the function compute_slabs.
  Those slabs are then combined to create the dataspace :
   for DSFile, we are expecting a list including all the slabs looking like
   [[startI_1, startJ_1, startK_1], [1,1,1], [nbI_1, nbJ_1, nbK_1], [1,1,1],
    [startI_2, startJ_2, startK_2], [1,1,1], [nbI_2, nbJ_2, nbK_2], [1,1,1], ...
    [startI_N, startJ_N, startK_N], [1,1,1], [nbI_N, nbJ_N, nbK_N], [1,1,1]]
   DSGlob me be the list of the tree dimensions sizes
   DSMmry and DSFrom have the same structure than flat / 1d dataspaces

  Mostly usefull for structured blocks.
  """
  slab_list  = compute_slabs(data_shape, distrib[0:2])
  dn_da    = distrib[1] - distrib[0]
  DSFILEDA = []
  for slab in slab_list:
    iS,iE, jS,jE, kS,kE = [item for bounds in slab for item in bounds]
    DSFILEDA.extend([[iS,jS,kS], [1,1,1], [iE-iS, jE-jS, kE-kS], [1,1,1]])
  DSMMRYDA = [[0]    , [1]    , [dn_da], [1]]
  DSFILEDA = list([list(DSFILEDA)])
  DSGLOBDA = [list(data_shape)]
  DSFORMDA = [[0]]
  return DSMMRYDA + DSFILEDA + DSGLOBDA + DSFORMDA

def create_flat_dataspace(distrib):
  """
  Create the most basic dataspace (1d / flat) for a given
  distribution.
  """
  dn_da    = distrib[1] - distrib[0]
  DSMMRYDA = [[0         ], [1], [dn_da], [1]]
  DSFILEDA = [[distrib[0]], [1], [dn_da], [1]]
  DSGLOBDA = [[distrib[2]]]
  DSFORMDA = [[0]]
  return DSMMRYDA + DSFILEDA + DSGLOBDA + DSFORMDA

def create_pe_dataspace(distrib):
  """
  Create a dataspace from a flat distribution, of elements,
  but adapted to "ParentElements" arrays ie (N,2) numpy arrays.
  """
  dn_pe    = distrib[1] - distrib[0]
  DSMMRYPE = [[0              , 0], [1, 1], [dn_pe, 2], [1, 1]]
  DSFILEPE = [[distrib[0], 0], [1, 1], [dn_pe, 2], [1, 1]]
  DSGLOBPE = [[distrib[2], 2]]
  DSFORMPE = [[1]]
  return DSMMRYPE + DSFILEPE + DSGLOBPE + DSFORMPE

def create_pointlist_dataspace(distrib):
  """
  Create a dataspace from a flat distribution, but adapted to "fake 2d" arrays
  ie (1,N) numpy arrays.
  Mostly usefull for PointList arrays and DataArray of the related BCDataSets.
  """
  dn_pl    = distrib[1] - distrib[0]
  DSMMRYPL = [[0,0          ], [1, 1], [1, dn_pl], [1, 1]]
  DSFILEPL = [[0, distrib[0]], [1, 1], [1, dn_pl], [1, 1]]
  DSGLOBPL = [[1, distrib[2]]]
  DSFORMPL = [[0]]
  return DSMMRYPL + DSFILEPL + DSGLOBPL + DSFORMPL

def create_data_array_filter(distrib, data_shape=None):
  """
  Create an hdf dataspace for the given distribution. The kind of
  dataspace depends of the data_shape optional argument, representing
  the size of the array for which the dataspace is created in each dimension:
  - If data_shape is None or a single value, dataspace is 1d/flat
  - If data_shape is a 2d list [1, N], a dataspace adpated to pointlist is created
  - In other cases (which should correspond to true 2d array or 3d array), the
    dataspace is create from combine method (flat in memory, block in file).
  """
  if data_shape is None or len(data_shape) == 1: #Unstructured
    hdf_data_space = create_flat_dataspace(distrib)
  elif len(data_shape) == 2 and data_shape[0] == 1:
    hdf_data_space = create_pointlist_dataspace(distrib)
  else: #Structured
    hdf_data_space = create_combined_dataspace(data_shape, distrib)

  return hdf_data_space

def gen_elemts(zone_tree):
  elmts_ini = I.getNodesFromType1(zone_tree, 'Elements_t')
  for elmt in elmts_ini:
    yield elmt

def getSubregionExtent(sub_region_node, zone):
  """
  Return the path of the node (starting from zone node) related to sub_region_node
  node (BC, GC or itself)
  """
  #if I.getNodeFromName1(sub_region_node, "BCRegionName") is not None:
  #  for zbc, bc in iterNodesWithParentsByMatching(zone, "ZoneBC_t/BC_t"):
  #    if I.getName(bc) == I.getValue(I.getNodeFromName1(sub_region_node, "BCRegionName")):
  #      return I.getName(zbc) + '/' + I.getName(bc)
  #elif I.getNodeFromName1(sub_region_node, "GridConnectivityRegionName") is not None:
  #  gc_pathes = ["ZoneGridConnectivity_t/GridConnectivity_t", "ZoneGridConnectivity_t/GridConnectivity1to1_t"]
  #  for gc_path in gc_pathes:
  #    for zgc, gc in iterNodesWithParentsByMatching(zone, gc_path):
  #      if I.getName(gc) == I.getValue(I.getNodeFromName1(sub_region_node, "GridConnectivityRegionName")):
  #        return I.getName(zgc) + '/' + I.getName(gc)
  #else:
  return I.getName(sub_region_node)
  
def create_zone_eso_elements_filter(elmt, zone_path, hdf_filter, mode):
  distrib_elmt = I.getVal(getDistribution(elmt, 'Element'))
  dn_elmt      = distrib_elmt[1] - distrib_elmt[0]

  # > For NGon only
  pe = I.getNodeFromName1(elmt, 'ParentElements')
  if pe:
    data_space = create_pe_dataspace(distrib_elmt)
    hdf_filter[f"{zone_path}/{I.getName(elmt)}/ParentElements"] = data_space
    if I.getNodeFromName1(elmt, 'ParentElementsPosition'):
      hdf_filter[f"{zone_path}/{I.getName(elmt)}/ParentElementsPosition"] = data_space

  eso = I.getNodeFromName1(elmt, 'ElementStartOffset')
  eso_path = None
  if eso:
    # Distribution for NGon/NFace -> ElementStartOffset is the same than DistrbutionFace, except
    # that the last proc have one more element
    n_elmt = distrib_elmt[2]
    if mode == 'read':
      dn_elmt_idx = dn_elmt + 1 # + int(distrib_elmt[1] == n_elmt)
    elif mode == 'write':
      dn_elmt_idx = dn_elmt + int((distrib_elmt[1] == n_elmt) and (distrib_elmt[0] != distrib_elmt[1]))
    DSMMRYESO = [[0              ], [1], [dn_elmt_idx], [1]]
    DSFILEESO = [[distrib_elmt[0]], [1], [dn_elmt_idx], [1]]
    DSGLOBESO = [[n_elmt+1]]
    DSFORMESO = [[0]]

    eso_path = zone_path+"/"+elmt[0]+"/ElementStartOffset"
    hdf_filter[eso_path] = DSMMRYESO + DSFILEESO + DSGLOBESO + DSFORMESO

  ec = I.getNodeFromName1(elmt, 'ElementConnectivity')
  if ec:
    if eso_path is None:
      raise RuntimeError("In order to load ElementConnectivity, the ElementStartOffset is mandatory")
    ec_path = zone_path+"/"+elmt[0]+"/ElementConnectivity"
    hdf_filter[ec_path] = partial(load_element_connectivity_from_eso, elmt, zone_path)

def load_element_connectivity_from_eso(elmt, zone_path, hdf_filter):
  distrib_ud   = getDistribution(elmt)
  distrib_elmt = I.getNodeFromName1(distrib_ud, 'Element')[1]
  dn_elmt      = distrib_elmt[1] - distrib_elmt[0]

  eso_n = I.getNodeFromName1(elmt, 'ElementStartOffset') # Maintenant il est chargÃ©
  if eso_n[1] is None: raise RuntimeError
  eso = eso_n[1]

  beg_face_vtx = eso[0]
  end_face_vtx = eso[eso.shape[0]-1]
  dn_face_vtx  = end_face_vtx - beg_face_vtx

  # print("beg_face_vtx::", beg_face_vtx)
  # print("end_face_vtx::", end_face_vtx)
  distrib_n  = None
  ec_size_n  = I.getNodeFromName1(elmt, 'ElementConnectivity#Size')
  if ec_size_n is not None:
    n_face_vtx = numpy.prod(ec_size_n[1])
  else:
    distrib_ud = getDistribution(elmt)
    distrib_n  = I.getNodeFromName1(distrib_ud, "ElementConnectivity")
    n_face_vtx = distrib_n[1][2]

  # print("n_face_vtx::", n_face_vtx)

  n_face      = distrib_elmt[2]
  dn_face_idx = dn_elmt + int(distrib_elmt[1] == n_face)
  DSMMRYEC = [[0           ], [1], [dn_face_vtx], [1]]
  DSFILEEC = [[beg_face_vtx], [1], [dn_face_vtx], [1]]
  DSGLOBEC = [[n_face_vtx ]]
  DSFORMEC = [[0]]

  ec_path = zone_path+"/"+elmt[0]+"/ElementConnectivity"
  hdf_filter[ec_path] = DSMMRYEC + DSFILEEC + DSGLOBEC + DSFORMEC

  if distrib_n is None:
    distrib = numpy.empty(3, dtype=eso.dtype)
    distrib[0] = beg_face_vtx
    distrib[1] = end_face_vtx
    distrib[2] = n_face_vtx
    I.newDataArray("ElementConnectivity", value=distrib, parent=distrib_ud)

def create_zone_std_elements_filter(elmt, zone_path, hdf_filter):
  distrib_elmt = I.getVal(getDistribution(elmt, 'Element'))
  dn_elmt      = distrib_elmt[1] - distrib_elmt[0]

  elmt_npe = I.eltName2EltNo(elmt)[1] # nbre de vertex de l'element

  DSMMRYElmt = [[0                       ], [1], [dn_elmt*elmt_npe], [1]]
  DSFILEElmt = [[distrib_elmt[0]*elmt_npe], [1], [dn_elmt*elmt_npe], [1]]
  DSGLOBElmt = [[distrib_elmt[2]*elmt_npe]]
  DSFORMElmt = [[0]]

  path = zone_path+"/"+elmt[0]+"/ElementConnectivity"
  hdf_filter[path] = DSMMRYElmt + DSFILEElmt + DSGLOBElmt + DSFORMElmt

  pe = I.getNodeFromName1(elmt, 'ParentElements')
  if pe:
    data_space = create_pe_dataspace(distrib_elmt)
    hdf_filter[f"{zone_path}/{I.getName(elmt)}/ParentElements"] = data_space
    if I.getNodeFromName1(elmt, 'ParentElementsPosition'):
      hdf_filter[f"{zone_path}/{I.getName(elmt)}/ParentElementsPosition"] = data_space

def create_zone_elements_filter(zone_tree, zone_path, hdf_filter, mode):
  """
  Prepare the hdf_filter for all the Element_t nodes found in the zone.
  """
  zone_elmts = gen_elemts(zone_tree)
  for elmt in zone_elmts:
    if elmt[1][0] == 22 or elmt[1][0] == 23:
      create_zone_eso_elements_filter(elmt, zone_path, hdf_filter, mode)
    elif elmt[1][0] == 20:
      raise ValueError('MIXED elements not implemented.')
    else:
      create_zone_std_elements_filter(elmt, zone_path, hdf_filter)

def create_zone_bc_filter(zone, zone_path, hdf_filter):
  """
  Fill up the hdf filter for the BC_t nodes present in
  the zone.
  Filter is created for the following nodes:
   - PointList (if present = unstruct. only)
   - All arrays founds in BCDataSets. Those arrays are supposed
     to be shaped as the PointList array. If a BCDataSet contains
     no PointList/PointRange node, the data is assumed to be consistent
     with the PointList/PointRange of the BC. Otherwise, the PointList/
     PointRange node of the BCDataSet is used to set the size of the BCData
     arrays. In this case, the PointList (if any) of the BCDataSet is
     written in the filter as well.
  """
  for zone_bc in I.getNodesFromType1(zone, 'ZoneBC_t'):
    zone_bc_path = zone_path+"/"+zone_bc[0]
    for bc in I.getNodesFromType1(zone_bc, 'BC_t'):
      bc_path = zone_bc_path+"/"+bc[0]

      distrib_bc = I.getVal(getDistribution(bc, 'Index'))

      bc_shape = plOrPrSize(bc)
      data_space = create_data_array_filter(distrib_bc, bc_shape)
      apply_dataspace_to_pointlist(bc, bc_path, data_space, hdf_filter)

      for bcds in I.getNodesFromType1(bc, "BCDataSet_t"):
        bcds_path = bc_path + "/" + bcds[0]
        distrib_bcds_n = getDistribution(bcds)

        if distrib_bcds_n is None: #BCDS uses BC distribution
          distrib_data = distrib_bc
          data_shape   = bc_shape
        else: #BCDS has its own distribution
          distrib_data = I.getNodeFromName1(distrib_bcds_n, 'Index')[1]
          data_shape = plOrPrSize(bcds)

        data_space_pl = create_data_array_filter(distrib_data, data_shape)
        #BCDataSet always use flat data array
        data_space_array = create_data_array_filter(distrib_data, [data_shape.prod()])
        apply_dataspace_to_pointlist(bcds, bcds_path, data_space_pl, hdf_filter)
        for bcdata in I.getNodesFromType1(bcds, 'BCData_t'):
          bcdata_path = bcds_path + "/" + bcdata[0]
          apply_dataspace_to_arrays(bcdata, bcdata_path, data_space_array, hdf_filter)


def create_zone_grid_connectivity_filter(zone, zone_path, hdf_filter):
  """
  Fill up the hdf filter for the GC_t nodes present in the zone.
  For unstructured GC (GridConnectivity_t), the filter is set up for
  the PointList and PointListDonor arrays.
  Structured GC (GridConnectivity1to1_t) are skipped since there is
  no data to load for these nodes.
  """
  for zone_gc in I.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
    zone_gc_path = zone_path+"/"+zone_gc[0]
    for gc in I.getNodesFromType1(zone_gc, 'GridConnectivity_t'):
      gc_path = zone_gc_path+"/"+gc[0]
      distrib_ia = I.getVal(getDistribution(gc, 'Index'))

      gc_shape   = plOrPrSize(gc)
      data_space = create_data_array_filter(distrib_ia, gc_shape)
      apply_dataspace_to_pointlist(gc, gc_path, data_space, hdf_filter)

def create_flow_solution_filter(zone, zone_path, hdf_filter):
  """
  Fill up the hdf filter for the FlowSolution_t nodes present in the
  zone. The size of the dataspace are computed from the pointList node
  if present, or using allCells / allVertex if no pointList is present.
  Filter is created for the arrays and for the PointList if present
  """
  distrib_vtx  = I.getVal(getDistribution(zone, 'Vertex'))
  distrib_cell = I.getVal(getDistribution(zone, 'Cell'))
  for flow_solution in I.getNodesFromType1(zone, 'FlowSolution_t'):
    flow_solution_path = zone_path + "/" + I.getName(flow_solution)
    grid_location = I.getNodeFromType1(flow_solution, 'GridLocation_t')
    if grid_location is None: grid_location = 'Vertex'
    else: grid_location = I.GetValue(grid_location)
    distrib_ud_n = getDistribution(flow_solution)
    if distrib_ud_n:
      distrib_data = I.getNodeFromName1(distrib_ud_n, 'Index')[1]
      data_shape = plOrPrSize(flow_solution)
      data_space_pl = create_data_array_filter(distrib_data, data_shape)
      data_space = create_data_array_filter(distrib_data, [data_shape.prod()])
      apply_dataspace_to_pointlist(flow_solution, flow_solution_path, data_space_pl, hdf_filter)
    elif grid_location == 'CellCenter':
      data_space = create_data_array_filter(distrib_cell, zone[1][:,1])
    elif grid_location == 'Vertex':
      data_space = create_data_array_filter(distrib_vtx, zone[1][:,0])
    else:
      raise RuntimeError(f"GridLocation {grid_location} is not allowed without PL")
    apply_dataspace_to_arrays(flow_solution, flow_solution_path, data_space, hdf_filter)

def create_zone_subregion_filter(zone, zone_path, hdf_filter):
  """
  Fill up the hdf filter for the ZoneSubRegion_t nodes present in
  the zone.
  The size of the dataspace are computed from
  - the corresponding BC or GC if the subregion is related to one of them
    (this information is given by a BCRegionName or GridConnectivityRegionName
    node)
  - the PointList / PointSize node of the subregion otherwise.
  Filter is created for the following nodes:
   - All arrays present in ZoneSubRegion;
   - PointList array if the zone is unstructured and if the subregion
     is not related to a BC/GC.
  """
  for zone_subregion in I.getNodesFromType1(zone, 'ZoneSubRegion_t'):
    zone_subregion_path = zone_path+"/"+zone_subregion[0]

    # Search matching region
    matching_region_path = getSubregionExtent(zone_subregion, zone)
    matching_region = I.getNodeFromPath(zone, matching_region_path)
    
    distrib_ud_n = getDistribution(matching_region)
    if not distrib_ud_n:
      raise RuntimeError("ZoneSubRegion {0} is not well defined".format(zone_subregion[0]))
    distrib_data = I.getNodeFromName1(distrib_ud_n, 'Index')[1]

    data_shape = plOrPrSize(matching_region)
    data_space_pl = create_data_array_filter(distrib_data, data_shape)
    data_space_ar = create_data_array_filter(distrib_data, [data_shape.prod()])

    apply_dataspace_to_pointlist(zone_subregion, zone_subregion_path, data_space_pl, hdf_filter)
    apply_dataspace_to_arrays(zone_subregion, zone_subregion_path, data_space_ar, hdf_filter)

def createZoneFilter(zone, zone_path, hdf_filter, mode):
  """
  Fill up the hdf filter for the following elements of the zone:
  Coordinates, Elements (NGon / NFace, Standards), FlowSolution
  (vertex & cells only), ZoneSubRegion, ZoneBC (including BCDataSet)
  and ZoneGridConnectivity.

  The bounds of the filter are determined by the :CGNS#Distribution
  node and, for the structured zones, by the size of the blocks.
  """
  # Coords
  distrib_vtx  = I.getVal(getDistribution(zone, 'Vertex'))
  all_vtx_dataspace = create_data_array_filter(distrib_vtx, zone[1][:,0])
  for grid_c in I.getNodesFromType1(zone, 'GridCoordinates_t'):
    grid_coord_path = zone_path + "/" + I.getName(grid_c)
    apply_dataspace_to_arrays(grid_c, grid_coord_path, all_vtx_dataspace, hdf_filter)
  create_zone_elements_filter(zone, zone_path, hdf_filter, mode)
  create_zone_bc_filter(zone, zone_path, hdf_filter)
  create_zone_grid_connectivity_filter(zone, zone_path, hdf_filter)
  create_flow_solution_filter(zone, zone_path, hdf_filter)
  create_zone_subregion_filter(zone, zone_path, hdf_filter)

def createTreeHdfFilter(dist_tree, hdf_filter, mode='read'):
  for base in I.getNodesFromType1(dist_tree, 'CGNSBase_t'):
    for zone in I.getNodesFromType1(base, 'Zone_t'):
      zone_path = "/"+I.getName(base)+"/"+I.getName(zone)
      createZoneFilter(zone, zone_path, hdf_filter, mode)

#==========================================================
# load
#==========================================================
def update_tree_with_partial_load_dict(dist_tree, partial_dict_load):
  for path, data in partial_dict_load.items():
    Node = I.getNodeFromPath(dist_tree, path)
    Node[1] = data

def loadTreeFromFilter(filename, dist_tree, comm, hdf_filter):
  # print("load_tree_from_filter")
  hdf_filter_with_dim = {key: value for (key, value) in hdf_filter.items() \
      if isinstance(value, (list, tuple))}

  partial_dict_load = C.convertFile2PartialPyTreeFromPath(filename, hdf_filter_with_dim, comm)
  update_tree_with_partial_load_dict(dist_tree, partial_dict_load)

  # > Match with callable
  hdf_filter_with_func = {key: value for (key, value) in hdf_filter.items() \
      if not isinstance(value, (list, tuple))}
  unlock_at_least_one = True
  while (len(hdf_filter_with_func) > 0 and unlock_at_least_one):
    # Update if you can
    next_hdf_filter = dict()
    unlock_at_least_one = False
    for key, f in hdf_filter_with_func.items():
      try:
        f(next_hdf_filter)
        unlock_at_least_one = True
      except RuntimeError: # Not ready yet
        pass
    partial_dict_load = C.convertFile2PartialPyTreeFromPath(filename, next_hdf_filter, comm)

    update_tree_with_partial_load_dict(dist_tree, partial_dict_load)
    hdf_filter_with_func = {key: value for (key, value) in next_hdf_filter.items() \
        if not isinstance(value, (list, tuple))}

  if unlock_at_least_one is False:
    raise RuntimeError("Something strange in the loading process")

def clean_distribution_info(dist_tree):
  for base in I.getNodesFromType1(dist_tree, 'CGNSBase_t'):
    for zone in I.getNodesFromType1(base, 'Zone_t'):
      I._rmNodesByName1(zone, ':CGNS#Distribution')
      for elmt in I.getNodesFromType1(zone, 'Elements_t'):
        I._rmNodesByName1(elmt, ':CGNS#Distribution')
        I._rmNodesByName1(elmt, 'ElementConnectivity#Size')
      for zone_bc in I.getNodesFromType1(zone, 'ZoneBC_t'):
        for bc in I.getNodesFromType1(zone_bc, 'BC_t'):
          I._rmNodesByName2(bc, ':CGNS#Distribution')
          I._rmNodesByName2(bc, 'PointList#Size')
      for zone_gc in I.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
        for gc in I.getNodesFromType1(zone_gc, 'GridConnectivity_t') + \
                  I.getNodesFromType1(zone_gc, 'GridConnectivity1to1_t'):
          I._rmNodesByName1(gc, ':CGNS#Distribution')
          I._rmNodesByName1(gc, 'PointList#Size')
      for zone_subregion in I.getNodesFromType1(zone, 'ZoneSubRegion_t'):
        I._rmNodesByName1(zone_subregion, ':CGNS#Distribution')
        I._rmNodesByName1(zone_subregion, 'PointList#Size')
      for zone_sol in I.getNodesFromType1(zone, 'FlowSolution_t'):
        I._rmNodesByName1(zone_sol, ':CGNS#Distribution')
        I._rmNodesByName1(zone_sol, 'PointList#Size')

def saveTreeFromFilter(filename, dist_tree, comm, hdf_filter):
  hdf_filter_with_dim  = {key: value for (key, value) in hdf_filter.items() if isinstance(value, list)}
  hdf_filter_with_func = {key: value for (key, value) in hdf_filter.items() if not isinstance(value, list)}

  next_hdf_filter = dict()
  for key, f in hdf_filter_with_func.items():
    f(hdf_filter_with_dim)

  #Dont save distribution info, but work on a copy to keep it for further use
  saving_dist_tree = I.copyRef(dist_tree)
  clean_distribution_info(saving_dist_tree)

  C.convertPyTree2FilePartial(saving_dist_tree, filename, comm, hdf_filter_with_dim, ParallelHDF=True)

#========================================================
# resume
#========================================================
#distTree = Part.loadCollectiveSizeTree(inputfile)
#Part._addDistributionInfo(distTree)
#hdf_filter = {}
#Part.createTreeHdfFilter(distTree, hdf_filter)
##skip_type_ancestors = [["Zone_t", "FlowSolution#EndOfRun", "Momentum*"],
##                       ["Zone_t", "ZoneSubRegion_t", "Velocity*"]]
#Part.loadTreeFromFilter(inputfile, distTree, Cmpi.KCOMM, hdf_filter)
