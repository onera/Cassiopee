/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepCheck_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" 
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "Geom2d_Curve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepGProp.hxx"
#include "GProp_GProps.hxx"

// ultimate (best) functions

// ============================================================================
// Return the nbPoints and ue for meshing E regular with hmax
// ============================================================================
E_Int __getParamHmax(const TopoDS_Edge& E, E_Float hmax, E_Int& nbPoints, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints=2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  nbPoints = (E_Int)round(L / hmax);
  nbPoints = std::max(nbPoints, E_Int(3));
  GCPnts_UniformAbscissa param(geomAdap, nbPoints, u0, u1);
  ue = new E_Float [nbPoints];
  for (E_Int i = 0; i < nbPoints; i++) ue[i] = param.Parameter(i+1);
  printf("L=%f hmax=%f nbPoints=" SF_D_ "\n", L, hmax, nbPoints); fflush(stdout);
  return 0;
}

// ============================================================================
// Return the nbPoints and ue for meshing E with given deflection
// ============================================================================
E_Int __getParamHausd(const TopoDS_Edge& E, E_Float hausd, E_Int& nbPoints, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  if (BRep_Tool::Degenerated(E)) 
  { 
    nbPoints=2; 
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  GCPnts_UniformDeflection param(geomAdap, hausd, u0, u1);
  nbPoints = param.NbPoints();
  ue = new E_Float [nbPoints];
  for (E_Int i = 0; i < nbPoints; i++) ue[i] = param.Parameter(i+1);
  printf("hausd=%f nbPoints=" SF_D_ "\n", hausd, nbPoints); fflush(stdout);
  return 0;
}

// ============================================================================
// Return ue for meshing with given param in [0,1]
// ============================================================================
E_Int __getParamExt(const TopoDS_Edge& E, E_Int nbPoints, E_Float* uext, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  ue = new E_Float [nbPoints];
  for (E_Int i = 0; i < nbPoints; i++) ue[i] = uext[i]*(u1-u0)+u0;
  //for (E_Int i = 0; i < nbPoints; i++) printf("%d : %f %f\n", i, ue[i], uext[i]);
  return 0;
}

// ============================================================================
// Mesh an edge from param [ue]
// ============================================================================
E_Int __meshEdge(const TopoDS_Edge& E, 
                 E_Int nbPoints, E_Float* ue,
                 K_FLD::FldArrayF& coords,
                 E_Boolean reverse)
{
  BRepAdaptor_Curve C0(E);
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  E_Float* pu = coords.begin(4);

  E_Float u; gp_Pnt Pt;
  for (E_Int i = 0; i < nbPoints; i++)
  {
    u = ue[i];
    C0.D0(u, Pt);
    if (reverse)
    { px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z(); pu[nbPoints-i-1] = u; }
    else
    { px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z(); pu[i] = u; }
  }
  return 0;
}

// ============================================================================
// Mesh an edge with fe by face
// ============================================================================
E_Int __meshEdgeByFace(const TopoDS_Edge& E, const TopoDS_Face& F,
                       E_Int& nbPoints, K_FLD::FldArrayF& fe, 
                       K_FLD::FldArrayF& coords,
                       E_Boolean reverse)
{
  BRepAdaptor_Curve C0(E);
  Handle(Geom_Surface) surf = BRep_Tool::Surface(F);
  Standard_Real aFirst = C0.FirstParameter(), aEnd=C0.LastParameter();
  Standard_Real pFirst = aFirst, pEnd=aEnd;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
  
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  //Standard_Real u0 = geomAdap.FirstParameter();
  //Standard_Real u1 = geomAdap.LastParameter();
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  E_Float* pu = coords.begin(4);
  E_Float* pv = coords.begin(5);

  // degenerated
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt; gp_Pnt2d Puv; 
    E_Float u = pFirst;
    C0.D0(u, Pt);

    for (E_Int i = 0; i < nbPoints; i++)
    {
      u = i*1./(nbPoints-1);
      u = u*(pEnd-pFirst)+pFirst;
      pCurve->D0(u, Puv);
      if (reverse)
      {
        px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z();
        pu[nbPoints-i-1] = Puv.X(); pv[nbPoints-i-1] = Puv.Y();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
        pu[i] = Puv.X(); pv[i] = Puv.Y();
      }
    }
    return 1;
  }
  
  // non degenerated    
  {
    gp_Pnt Pt; gp_Pnt2d Puv; 
    for (E_Int i = 1; i <= nbPoints; i++)
    {
      E_Float u = fe(i-1,4);
      C0.D0(u, Pt);
      pCurve->D0(u, Puv);
      if (reverse)
      { 
        px[nbPoints-i] = Pt.X(); py[nbPoints-i] = Pt.Y(); pz[nbPoints-i] = Pt.Z();
        pu[nbPoints-i] = Puv.X(); pv[nbPoints-i] = Puv.Y(); 
      }
      else
      {
        px[i-1] = Pt.X(); py[i-1] = Pt.Y(); pz[i-1] = Pt.Z();
        pu[i-1] = Puv.X(); pv[i-1] = Puv.Y(); 
      }
    }
  }
  return 0;
}

// ============================================================================
/* Mesh one edge of CAD, return STRUCT
   hmax
   hausd
   hmax + hausd
   external param ue
 */
// ============================================================================
PyObject* K_OCC::meshOneEdge(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int i;
  E_Float hmax; E_Float hausd; PyObject* externalEdge;
  if (!PYPARSETUPLE_(args, O_ I_ RR_ O_, &hook, &i, &hmax, &hausd, &externalEdge)) return NULL;
    
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  // Get the ith edge (start 1)
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  const TopoDS_Edge& E = TopoDS::Edge(edges(i));
  
  E_Int nbPoints = 0; // nbre of points of discretized edge
  E_Float* ue; // edge param

  if (hmax > 0 && hausd < 0 && externalEdge == Py_None) // pure hmax
  {
    __getParamHmax(E, hmax, nbPoints, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }
  else if (hmax < 0 && hausd > 0 && externalEdge == Py_None) // pure hausd
  {
    __getParamHausd(E, hausd, nbPoints, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }
  else if (hmax > 0 && hausd > 0 && externalEdge == Py_None) // mix hmax + hausd
  {
    // pour l'instant on retourne hmax comme pour les mailleurs precedents
    __getParamHmax(E, hmax, nbPoints, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }
  else if (externalEdge != Py_None)
  {
    E_Int ni, nj, nk;
    K_FLD::FldArrayF* fe; K_FLD::FldArrayI* ce;
    char* varString; char* eltType;
    K_ARRAY::getFromArray3(externalEdge, varString, fe, ni, nj, nk, ce, eltType);
    E_Int pos = K_ARRAY::isNamePresent("s", varString);
    if (pos == -1) pos = K_ARRAY::isNamePresent("u", varString);
    if (pos == -1)
    {
      RELEASESHAREDS(externalEdge, fe);
      PyErr_SetString(PyExc_ValueError,
                      "getParam: can't find parameter field (s or u) in array.");
      return NULL;
    }
    E_Float* uext = fe->begin(pos+1);
    __getParamExt(E, ni, uext, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", ni, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, ni, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    RELEASESHAREDS(externalEdge, fe);
    return o;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "meshOneEdge: invalid input.");
    return NULL;
  }
  return NULL;
}

// ============================================================================
/* Mesh edges for one face by wires 
   sort by wire, reverse for orientation, close by wire */
// ============================================================================
PyObject* K_OCC::meshEdgesOfFace(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int noFace;
  PyObject* discretizedEdges;
  if (!PYPARSETUPLE_(args, O_ I_ O_, &hook, &noFace, &discretizedEdges)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
  // outer wire
  const TopoDS_Wire& OW = ShapeAnalysis::OuterWire(F);
  
  TopAbs_Orientation forientation = F.Orientation();
  //if (forientation == TopAbs_FORWARD) printf("face orientation=forward\n");
  //else if (forientation == TopAbs_REVERSED) printf("face orientation=reversed\n");

  PyObject* out = PyList_New(0); // sortie par wire
  
  for (expl.Init(surfaces(noFace), TopAbs_WIRE); expl.More(); expl.Next())
  {
    PyObject* le = PyList_New(0); // sortie par edges
    const TopoDS_Wire& W = TopoDS::Wire(expl.Current());
    //printf("getting a wire\n");
    //BRepCheck_Wire check = BRepCheck_Wire(W);
    BRepCheck_Wire check(W);
    
    //BRepCheck_Status status = check.Orientation(F);
    //if (status == BRepCheck_NoError) printf("wire correctly oriented\n");
    //else printf("WIRE BADLY oriented\n");
    BRepCheck_Status wclosed = check.Closed();
    //if (wclosed == BRepCheck_NoError) printf("wire correctly closed\n");
    //else printf("WIRE BADLY closed\n");
    
    //E_Boolean isOuter = false;
    //if (W == OW) { printf("is outer wire\n"); isOuter = true; }
    //status = check.CheckOuterBound();
    //if (status == BRepCheck_NoError) printf("is outer (test2)\n");
    //else printf("is inner (test2)\n");

    //E_Float surface = ShapeAnalysis::ContourArea(W);
    //printf("wire surface=%f\n", surface);	

    // Compute oriented surface to orient externe = CCW, intern = CW
    //BRepBuilderAPI_MakeFace faceMaker(F, W);
    //faceMaker.Build();
    //const TopoDS_Face& F2 = faceMaker.Face();
    //GProp_GProps SProps;
    //BRepGProp::SurfaceProperties(F2, SProps);
    //Standard_Real area = SProps.Mass();
    //gp_Pnt center = SProps.CentreOfMass();
    //printf("Signed surface = %f\n", area);
    //printf("Center of mass = %f %f %f\n", center.X(), center.Y(), center.Z());

    TopAbs_Orientation worientation = W.Orientation();
    //if (worientation == TopAbs_FORWARD) printf("wire orientation=forward\n");
    //else if (worientation == TopAbs_REVERSED) printf("wire orientation=reversed\n");
    //else if (worientation == TopAbs_INTERNAL) printf("wire orientation=internal\n");
    //else if (worientation == TopAbs_EXTERNAL) printf("wire orientation=external\n");
    worientation = TopAbs_FORWARD; // ignored

    //TopExp_Explorer expl2;
    BRepTools_WireExplorer expl2;
    std::vector<FldArrayF*> discreteWire;
    std::vector<PyObject*> discreteWire2;

    // adapt hmax from wire length
    /*
    E_Float Lw = 0.;
    for (expl2.Init(W, F); expl2.More(); expl2.Next())
    {
      const TopoDS_Edge& E = TopoDS::Edge(expl2.Current());  
      Lw += __getLength(E);
    }
    E_Float alpha = Lw / (41*hmax);
    alpha = std::max(alpha, -2.);
    alpha = std::min(alpha, 2.);
    hmaxw = hmax + alpha * hmax;
    */
   
    //for (expl2.Init(W, TopAbs_EDGE); expl2.More(); expl2.Next())
    for (expl2.Init(W, F); expl2.More(); expl2.Next())
    {
      //printf("getting edges of wire\n");
      
      const TopoDS_Edge& E = TopoDS::Edge(expl2.Current());

      // Get the edge in the list of discretized edges
      E_Int eno = 0;
      for (eno = 1; eno <= edges.Extent(); eno++)
      {
        const TopoDS_Edge& El = TopoDS::Edge(edges(eno));
        if (El.TShape() == E.TShape()) break;
      }

      // Extract discrete edge
      PyObject* ge = PyList_GetItem(discretizedEdges, eno-1);
      E_Int ni, nj, nk;
      K_FLD::FldArrayF* fe; K_FLD::FldArrayI* ce;
      char* varString; char* eltType;
      K_ARRAY::getFromArray3(ge, varString, fe, ni, nj, nk, ce, eltType);
      //E_Float* ue = fe->begin(4);
      E_Int nbPoints = ni;
      
      // get orientation      
      TopAbs_Orientation eorientation = E.Orientation();
      //if (eorientation == TopAbs_FORWARD) printf("edge orientation=forward\n");
      //else if (eorientation == TopAbs_REVERSED) printf("edge orientation=reversed\n");
      //else if (eorientation == TopAbs_INTERNAL) printf("edge orientation=internal\n");
      //else if (eorientation == TopAbs_EXTERNAL) printf("edge orientation=external\n");

      // create array
      PyObject* o = K_ARRAY::buildArray2(5, "x,y,z,u,v", nbPoints, 1, 1, 1);
      FldArrayF* f; K_ARRAY::getFromArray2(o, f);
      discreteWire.push_back(f);
      discreteWire2.push_back(o);

#define REVERSEFACE
      E_Boolean reversed = false;
      if (forientation == TopAbs_FORWARD)
      {
        if (worientation == TopAbs_FORWARD)
        {
          if (eorientation == TopAbs_FORWARD) reversed = false;
          else reversed = true; 
        }
        else
        {
          if (eorientation == TopAbs_FORWARD) reversed = true;
          else reversed = false; 
        }
      }
      else
      {
#ifdef REVERSEFACE
        if (worientation == TopAbs_FORWARD)
        {
          if (eorientation == TopAbs_FORWARD) reversed = true;
          else reversed = false;
        }
        else
        {
          if (eorientation == TopAbs_FORWARD) reversed = false;
          else reversed = true; 
        }
#else   
        if (worientation == TopAbs_FORWARD)
        {
          if (eorientation == TopAbs_FORWARD) reversed = false;
          else reversed = true;
        }
        else
        {
          if (eorientation == TopAbs_FORWARD) reversed = true;
          else reversed = false;
        }
#endif   
      }
      
      __meshEdgeByFace(E, F, nbPoints, *fe, *f, reversed);

      //RELEASESHAREDS(o, f); // done later
      if (forientation == TopAbs_FORWARD)
      {
        if (worientation == TopAbs_FORWARD) PyList_Append(le, o); 
        else PyList_Insert(le, 0, o); 
      }
      else
      {
#ifdef REVERSEFACE
        if (worientation == TopAbs_REVERSED) PyList_Append(le, o);
        else PyList_Insert(le, 0, o);
#else   
        if (worientation == TopAbs_FORWARD) PyList_Append(le, o); 
        else PyList_Insert(le, 0, o);
#endif   
      }
      Py_DECREF(o);
    }

    if (forientation == TopAbs_FORWARD)
    {
      if (worientation == TopAbs_REVERSED)
      { 
        std::reverse(discreteWire.begin(), discreteWire.end());  
        std::reverse(discreteWire2.begin(), discreteWire2.end());  
      }
    }
    else
    {
#ifdef REVERSEFACE
      if (worientation == TopAbs_FORWARD)
      {
        std::reverse(discreteWire.begin(), discreteWire.end());  
        std::reverse(discreteWire2.begin(), discreteWire2.end());  
      }
#else 
      if (worientation == TopAbs_REVERSED)
      { 
        std::reverse(discreteWire.begin(), discreteWire.end());  
        std::reverse(discreteWire2.begin(), discreteWire2.end());  
      }
#endif 
    }

    size_t nedges = discreteWire.size();
    if (wclosed == BRepCheck_NoError)
    {
      //printf("closing the discrete wire\n");
      // ferme le wire discret
      // Les edges sont normalement correctement orientes (tail to head)
      for (size_t i = 0; i < nedges; i++)
      {
        if (i < nedges-1)
        {
          FldArrayF& f1 = *(discreteWire[i]);
          FldArrayF& f2 = *(discreteWire[i+1]);
          E_Int np = f1.getSize();
          f1(np-1,1) = f2(0,1);
          f1(np-1,2) = f2(0,2);
          f1(np-1,3) = f2(0,3);
          f1(np-1,4) = f2(0,4);
          f1(np-1,5) = f2(0,5);
        }
        else
        {
          FldArrayF& f1 = *(discreteWire[i]);
          FldArrayF& f2 = *(discreteWire[0]);
          E_Int np = f1.getSize();
          f1(np-1,1) = f2(0,1);
          f1(np-1,2) = f2(0,2);
          f1(np-1,3) = f2(0,3);
          f1(np-1,4) = f2(0,4);
          f1(np-1,5) = f2(0,5);
        }
      }
    }
    for (size_t i = 0; i < nedges; i++) 
    { RELEASESHAREDS(discreteWire2[i], discreteWire[i]); }
    PyList_Append(out, le); Py_DECREF(le);
  }
  
  return out;
}

// Return face orientation in CAD
PyObject* K_OCC::getFaceOrientation(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int noFace;
  if (!PYPARSETUPLE_(args, O_ I_, &hook, &noFace)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));

  TopAbs_Orientation forientation = F.Orientation();
  if (forientation == TopAbs_FORWARD) return Py_BuildValue("l", 1);
  else if (forientation == TopAbs_REVERSED) return Py_BuildValue("l", 0);
  else return Py_BuildValue("l", -1); // unknown
}
