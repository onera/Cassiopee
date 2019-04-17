# Conservation des fonctions de lecture/ecriture partielle / MPI

# -- convertPyTree2FileMPI
def convertPyTree2FileMPI(t, fileName, comm, SkeletonTree, ParallelHDF=False,
                             format=None, isize=4, rsize=8,
                             endian='big', colormap=0, dataFormat='%.9e '):
  """Convert a pyTree to a file.
  Usage: convertPyTree2File(t, fileName, format, options)"""

  # > GardeFou
  if t == []: print('Warning: convertPyTree2File: nothing to write.'); return
  format = 'bin_hdf'

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > First step : Prepare a dictonary of Filter and a dictionnary of Property
  #   in order to prepare for each procs the data to write ...
  Filter = dict()
  Proper = dict()
  Elmts  = dict()
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  import CGNS.PAT.cgnsutils as CGU
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  printTree(t)

  for Zone in Internal.getZones(t):
    pathsToArray = Internal.getPathsFromType(Zone, 'IndexRange_t')
    print(pathsToArray)
    pathsToArray  = CGU.getPathsByTypeSet(Zone, 'IndexRange_t')
    print(pathsToArray)
    
  # > The path who wants to effectivly write is DataArray_t
  pathsToArray = Internal.getPathsFromType(t, 'IndexRange_t')
  pathsToArray += Internal.getPathsFromType(t, 'DataArray_t')
  print(pathsToArray)

  pathsToArray  = CGU.getPathsByTypeSet(t, 'IndexRange_t')
  pathsToArray += CGU.getPathsByTypeSet(t, 'DataArray_t')
  print(pathsToArray)
  for path in pathsToArray:
     print(path)
     node = Internal.getNodeFromPath(t, path)

     # > The ideo is to make a global Dataspace full for the current proc and void for the other
     NbE = list(node[1].shape)
     Beg = [0]*len(NbE); Sti = [1]*len(NbE); Blk = [1]*len(NbE)
     DataSpaceMMRY = [Beg, Sti, NbE, Blk]

     # > Partial filter (Voluntary not fill completely see after ...)
     Filter[path] = DataSpaceMMRY

     # > You need to setup all label to rebuild the tree
     Label    = []
     ListPath = path.split("/")[1:]
     topnode = t
     for i,l in enumerate(ListPath):
        topnode = Internal.getNodeFromName1(topnode, l)
        Label.append(topnode[3])

     # > Fill property dictionnary
     Proper[path] = {'ProcNumber' : comm.Get_rank(),
                     'DataType'   : node[1].dtype,
                     'DataShape'  : node[1].shape,
                     'Label'      : Label}
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  pathsToElmts = Internal.getPathsFromType(t, 'Elements_t')
  for p in pathsToElmts:
    node = Internal.getNodeFromPath(t, p)
    Elmts[p] = [node[0], node[1], [], node[3]]

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Merge filter among procs ...
  ListFilter = comm.gather(Filter, root=0)
  ListProper = comm.gather(Proper, root=0)
  ListElmts  = comm.gather(Elmts , root=0)

  if comm.Get_rank() == 0:
    DictFilterAll = dict()
    DictProperAll = dict()
    DictElmtsAll  = dict()
    for Proc in ListFilter:
       for Path in Proc:
         DictFilterAll[Path] = Proc[Path]
    for Proc in ListProper:
      for Path in Proc:
         DictProperAll[Path] = Proc[Path]
    for Proc in ListElmts:
      for Path in Proc:
         DictElmtsAll[Path] = Proc[Path]
  else:
    DictFilterAll = None
    DictProperAll = None
    DictElmtsAll  = None
  Filter   = comm.bcast(DictFilterAll, root=0)
  Proper   = comm.bcast(DictProperAll, root=0)
  Elmts    = comm.bcast(DictElmtsAll , root=0)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Perform sort of the receive dictionnary
  for path in Filter:
    if(Proper[path]['ProcNumber'] != comm.Get_rank()):
       # > Change the global size
       TMP    = Filter[path]
       DataSpaceGLOB = TMP[2]
       TMP[2] = [0]*len(TMP[2])
       Filter[path] = TMP+TMP+[DataSpaceGLOB]
    else:
      Filter[path] = Filter[path]+Filter[path]+[Filter[path][2]]
    # print path, Filter[path]
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > A partir de l'arbre on recreer tout les noeuds
  for path in Filter:
    # Node = Internal.getNodeFromPath(t, path)
    ListPath = path.split("/")[1:]
    topnode = t
    for i,l in enumerate(ListPath):
      if Internal.getNodeFromName1(topnode, l) is None:
        if i == len(ListPath)-1:
          shape   = [0]*len(list(Proper[path]['DataShape']))
          lvalue  = numpy.empty(shape, dtype=Proper[path]['DataType'])
          # topnode = Internal.createUniqueChild(topnode, l, 'DataArray_t', value=lvalue)
          topnode = Internal.createUniqueChild(topnode, l, Proper[path]['Label'][i], value=lvalue)
        else:
          topnode = Internal.createUniqueChild(topnode, l, Proper[path]['Label'][i])
      else:
        topnode = Internal.getNodeFromName1(topnode, l)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Make a proper Skeleton tree wit
  SkeletonTree2 = Internal.copyRef(t)
  for path in Filter:
    Node = Internal.getNodeFromPath(SkeletonTree2, path)
    Node[1] = None

  if comm.Get_rank() == 0:
    for path in Elmts:
      ListPath = path.split("/")[1:-1]
      EndPath  = path.split("/")[-1]
      topnode  = SkeletonTree2
      for l in ListPath:
        print(l)
        topnode = Internal.getNodeFromName1(topnode, l)
      print(path)
      if Internal.getNodeFromName1(topnode, EndPath) is None:
        Internal._addChild(topnode, Elmts[path])
      else:
        Node = Internal.getNodeFromName1(topnode, EndPath)
        Node[1] = Elmts[path][1]

  for Zone in Internal.getZones(SkeletonTree2):
     Zone[1] = Internal.getNodeFromName2(SkeletonTree, Zone[0])[1]

  SkeletonTree = Internal.merge([SkeletonTree2, SkeletonTree])

  for path in Filter:
    Node = Internal.getNodeFromPath(SkeletonTree, path)
    Node[1] = None

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if comm.Get_rank() == 0: printTree(SkeletonTree)

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Write skeleton
  if comm.Get_rank() == 0:
    convertPyTree2File(SkeletonTree, fileName, format)
  comm.barrier()
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Write data in filter in file (With creation of DataSpace )
  skeletonData = None  # Skeleton Data is inecfective (Normaly)
  Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  return None

# -- convertPyTree2FilePartial
def convertPyTree2FilePartial(t, fileName, comm, Filter, ParallelHDF=False,
                              format=None, isize=4, rsize=8,
                              endian='big', colormap=0, dataFormat='%.9e '):
  """Convert a pyTree to a file.
  Usage: convertPyTree2File(t, fileName, format, options)"""
  if t == []: print('Warning: convertPyTree2File: nothing to write.'); return
  format = 'bin_hdf'

  if not ParallelHDF:
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write Tree Data execpt Data in Filter
    SkeletonTree = Internal.copyRef(t)
    for path in Filter:
      print(path)
      Node = Internal.getNodeFromPath(SkeletonTree, path)
      Node[1] = None
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Si MPI Mode Off (HDF Not Parallel)
    if comm.Get_rank() == 0:
      convertPyTree2File(SkeletonTree, fileName, format)
      # > Fill up Dimension
      skeletonData = None
      Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wait for Skeleton write
    comm.barrier()

    # > Set skeletonData to Not None
    skeletonData = []

    # > Cette maniere de faire provoque de la non reproductibilite ...
    # Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > On serialize
    for lock in range(comm.Get_size()):
      if lock == comm.Get_rank():
        Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
      comm.barrier()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  else:  # > Si MPI Mode Off (HDF Not Parallel)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write Tree Data execpt Data in Filter
    SkeletonTree = Internal.copyRef(t)
    for path in Filter:
      Node = Internal.getNodeFromPath(SkeletonTree, path)
      Node[1] = None
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if comm.Get_rank() == 0:
      convertPyTree2File(SkeletonTree, fileName, format)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # > On wait l'ecriture Skelette ...
    comm.barrier()

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write data in filter in file (With creation of DataSpace )
    skeletonData = None  # Skeleton Data is inefective (Normaly)
    Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# -- convertFile2PyTree
# def convertFile2PartialPyTree(fileName, comm, format=None, nptsCurve=20, nptsLine=2,
#                                         density=-1., skeletonData=None):
#   """Convert a file to pyTree.
#   Usage: convertFile2PyTree(fileName, format, options)"""
#   if format is None:
#     format = Converter.convertExt2Format__(fileName); autoTry = True
#   else: autoTry = False
#   try: file = open(fileName, 'r')
#   except: raise IOError("convertFile2PyTree: file %s not found."%fileName)

#   t = Converter.converter.convertFile2PartialPyTree(fileName, format, skeletonData, comm)
#   t = Internal.createRootNode(children=t[2])
#   return t

def convertFile2PyTreeFromPath(fileName, Filter,
                               format=None, nptsCurve=20, nptsLine=2,
                               density=-1., skeletonData=None):
  """Convert a file to pyTree.
  Usage: convertFile2PyTree(fileName, format, options)"""
  if format is None:
    format = Converter.convertExt2Format__(fileName); autoTry = True
  else: autoTry = False
  try: file = open(fileName, 'r')
  except: raise IOError("convertFile2PartialPyTreeFromPath: file %s not found."%fileName)

  t = Converter.converter.convertFile2PyTreeFromPath(fileName, format, Filter)
  # t = Internal.createRootNode(children=t[2])
  return t
