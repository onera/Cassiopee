// Create packet hook
#define CREATEHOOK \
  PyObject* hook; \
  void** packet = new void* [6]; \
  packet[0] = NULL; \
  packet[1] = NULL; \
  packet[2] = NULL; \
  packet[3] = NULL; \
  packet[4] = NULL; \
  packet[5] = NULL; \
  hook = PyCapsule_New(packet, NULL, NULL);

// GET PACKET
#define GETPACKET \
  void** packet = (void**) PyCapsule_GetPointer(hook, NULL);

// GET PACKET+SHAPE
#define GETSHAPE \
  void** packet = (void**) PyCapsule_GetPointer(hook, NULL); \
  TopoDS_Shape* shape = (TopoDS_Shape*)packet[0];

// GET MAPS
#define GETMAPSURFACES \
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

#define GETMAPEDGES \
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];

// SET MAPS
#define SETMAPSURFACES \
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1]; \
  delete ptr; \
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape(); \
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf); \
  packet[1] = sf;

#define SETMAPEDGES \
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2]; \
  delete ptr2; \
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape(); \
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se); \
  packet[2] = se;

#define SETSHAPE(newshp) \
  packet[0] = newshp; \
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1]; \
  delete ptr; \
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape(); \
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf); \
  packet[1] = sf; \
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2]; \
  delete ptr2; \
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape(); \
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se); \
  packet[2] = se;
