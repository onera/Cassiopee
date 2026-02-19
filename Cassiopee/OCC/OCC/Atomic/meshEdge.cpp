/*    
    Copyright 2013-2026 ONERA.

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

#define MAXNBPTSPEREDGE 100

// Open cascade snippets pour recuperer une pCurve
//Standard_Real aFirst, aLast, aPFirst, aPLast;
//Handle(Geom_Curve) aCurve3d = BRep_Tool::Curve(anEdge, aFirst, aLast);
//Handle(Geom2d_Curve) aPCurve = BRep_Tool::CurveOnSurface(anEdge, aFace, aPFirst, aPLast);

// ============================================================================
// Mesh an edge with [equal distance hmax] of nbPoints
// ============================================================================
E_Int __meshEdge1(const TopoDS_Edge& E, 
                  E_Int& nbPoints, K_FLD::FldArrayF& coords, E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);

  // degenerated
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt;
    C0.D0(0., Pt);
    for (E_Int i = 0; i < nbPoints; i++)
    {
      px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
    }
    return 1;
  }

  // non degenerated
  GCPnts_UniformAbscissa unifAbs(geomAdap, int(nbPoints), u0, u1);
  if (!unifAbs.IsDone()) return 1;
  if (nbPoints != unifAbs.NbPoints()) return 1;
    
  {
    gp_Pnt Pt;
    for (Standard_Integer i = 1; i <= nbPoints; i++)
    {
      C0.D0(unifAbs.Parameter(i), Pt);
      if (reverse)
      {px[nbPoints-i] = Pt.X(); py[nbPoints-i] = Pt.Y(); pz[nbPoints-i] = Pt.Z();}
      else
      {px[i-1] = Pt.X(); py[i-1] = Pt.Y(); pz[i-1] = Pt.Z();}
    }
  }
  return 0;
}

// Mesh an edge with [Equal distance] of nbPoints by face
E_Int __meshEdgeByFace1(const TopoDS_Edge& E, const TopoDS_Face& F,
                        E_Int& nbPoints, K_FLD::FldArrayF& coords,
                        E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  Handle(Geom_Surface) surf = BRep_Tool::Surface(F);
  Standard_Real aFirst = C0.FirstParameter(), aEnd=C0.LastParameter();
  Standard_Real pFirst = aFirst, pEnd=aEnd;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
  
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  E_Float* pu = coords.begin(4);
  E_Float* pv = coords.begin(5);

  // degenerated
  if (BRep_Tool::Degenerated(E))
  {
    printf("edge is degenerated\n");
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
  GCPnts_UniformAbscissa unifAbs(geomAdap, int(nbPoints), u0, u1);
  if (!unifAbs.IsDone()) return 1;
  if (nbPoints != unifAbs.NbPoints()) return 1;
    
  {
    gp_Pnt Pt; gp_Pnt2d Puv; 
    for (E_Int i = 1; i <= nbPoints; i++)
    {
      E_Float u = unifAbs.Parameter(i);
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
// Maille un edge avec NbPoints et une [reg param] 
// ============================================================================
E_Int __meshEdge2(const TopoDS_Edge& E,  
                  E_Int& nbPoints, K_FLD::FldArrayF& coords, 
                  E_Bool reverse)
{  
  BRepAdaptor_Curve C0(E);
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt; gp_Pnt2d Puv; 
    E_Float u = aFirst;
    C0.D0(u, Pt);

    for (E_Int i = 0; i < nbPoints; i++)
    {
      u = i*1./(nbPoints-1);
      u = u*(aEnd-aFirst)+aFirst;
      if (reverse)
      {
        px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
      }
    }
    return 0;
  }

  // non degenerated case
  {
    gp_Pnt Pt; gp_Pnt2d Puv; E_Float u;
    for (E_Int i = 0; i < nbPoints; i++)
    {
      u = i*1./(nbPoints-1);
      u = u*(aEnd-aFirst)+aFirst;
      C0.D0(u, Pt);
      if (reverse)
      { 
        px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
      }
    }
  }
  return 0;
}

// reg param by face
E_Int __meshEdgeByFace2(const TopoDS_Edge& E, const TopoDS_Face& F,  
                        E_Int& nbPoints, K_FLD::FldArrayF& coords,
                        E_Bool reverse)
{  
  BRepAdaptor_Curve C0(E);
  Handle(Geom_Surface) surf = BRep_Tool::Surface(F);
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  Standard_Real pFirst=aFirst, pEnd=aEnd;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  E_Float* pu = coords.begin(4);
  E_Float* pv = coords.begin(5);
  
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

  // non degenerated case
  {
    gp_Pnt Pt; gp_Pnt2d Puv; E_Float u;
    for (E_Int i = 0; i < nbPoints; i++)
    {
      u = i*1./(nbPoints-1);
      u = u*(aEnd-aFirst)+aFirst;
      C0.D0(u, Pt); pCurve->D0(u, Puv);
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
  }
  return 0;
}

// ============================================================================
// Maille un edge avec la parametrisation donnee dans [ue] input
// ============================================================================
E_Int __meshEdge3(const TopoDS_Edge& E,  
                  E_Int ni, E_Float* ue, K_FLD::FldArrayF& coords,
                  E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  Standard_Real pFirst=aFirst;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt; 
    E_Float u = pFirst;
    C0.D0(u, Pt);

    for (E_Int i = 0; i < ni; i++)
    {
      if (reverse)
      {
        px[ni-i-1] = Pt.X(); py[ni-i-1] = Pt.Y(); pz[ni-i-1] = Pt.Z();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
      }
    }
    return 0;
  }

  // non degenerated case
  {
    gp_Pnt Pt; E_Float u;
    for (E_Int i = 0; i < ni; i++)
    {
      u = ue[i];
      u = u*(aEnd-aFirst)+aFirst;
      C0.D0(u, Pt);
      if (reverse)
      {
        px[ni-i-1] = Pt.X(); py[ni-i-1] = Pt.Y(); pz[ni-i-1] = Pt.Z();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
      }
    }
  }
  return 0;
}

// ============================================================================
// Maille un edge avec la parametrisation donnee dans [ue] input sur une face F
// retourne les coords des edges et le u,v sur la face
// ============================================================================
E_Int __meshEdgeByFace3(const TopoDS_Edge& E, const TopoDS_Face& F,  
                        E_Int ni, E_Float* ue, K_FLD::FldArrayF& coords,
                        E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  Standard_Real pFirst=aFirst, pEnd=aEnd;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  Handle(Geom_Surface) surf = BRep_Tool::Surface(F);
  Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);
  E_Float* pu = coords.begin(4);
  E_Float* pv = coords.begin(5);
  
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt; gp_Pnt2d Puv; 
    E_Float u = pFirst;
    C0.D0(u, Pt); pCurve->D0(u, Puv);
    for (E_Int i = 0; i < ni; i++)
    {
      if (reverse)
      {
        px[ni-i-1] = Pt.X(); py[ni-i-1] = Pt.Y(); pz[ni-i-1] = Pt.Z();
        pu[ni-i-1] = Puv.X(); pv[ni-i-1] = Puv.Y();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
        pu[i] = Puv.X(); pv[i] = Puv.Y();
      }
    }
    return 0;
  }

  // non degenerated case
  {
    gp_Pnt Pt; gp_Pnt2d Puv; E_Float u;
    for (E_Int i = 0; i < ni; i++)
    {
      u = ue[i];
      u = u*(aEnd-aFirst)+aFirst;
      C0.D0(u, Pt);
      pCurve->D0(u, Puv);
      if (reverse)
      {
        px[ni-i-1] = Pt.X(); py[ni-i-1] = Pt.Y(); pz[ni-i-1] = Pt.Z();
        pu[ni-i-1] = Puv.X(); pv[ni-i-1] = Puv.Y();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
        pu[i] = Puv.X(); pv[i] = Puv.Y();
      }
    }
  }
  return 0;
}

// ============================================================================
// Mesh an edge with [deflection]
// ============================================================================
E_Int __meshEdge4(const TopoDS_Edge& E, 
                  E_Float hausd, E_Float hmax,
                  E_Int& nbPoints, K_FLD::FldArrayF& coords,
                  E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  E_Float* px = coords.begin(1);
  E_Float* py = coords.begin(2);
  E_Float* pz = coords.begin(3);

  // degenerated
  if (BRep_Tool::Degenerated(E))
  {
    gp_Pnt Pt;
    E_Float u = u0;
    C0.D0(u, Pt);
    for (E_Int i = 0; i < nbPoints; i++)
    {
      if (reverse)
      {
        px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z();
      }
      else
      {
        px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
      }
    }
    return 1;
  }

  // non degenerated
  GCPnts_UniformDeflection param(geomAdap, hausd, u0, u1);
  //GCPnts_TangentialDeflection param(GeomAdap, 10., hausd, 2, 1.e-9, hmax);
  nbPoints = param.NbPoints();
  nbPoints = std::min(nbPoints, E_Int(MAXNBPTSPEREDGE)); // hard limit
  printf("hausd edge: nbPoints=" SF_D_ " coords=" SF_D_ "\n", nbPoints, coords.getSize()); fflush(stdout);
  if (nbPoints != coords.getSize()) exit(0);

  {
    gp_Pnt Pt; E_Float u;
    for (E_Int i = 1; i <= nbPoints; i++)
    {
      u = param.Parameter(i);
      C0.D0(u, Pt);
      if (reverse)
      {
        px[nbPoints-i] = Pt.X(); py[nbPoints-i] = Pt.Y(); pz[nbPoints-i] = Pt.Z();
      }
      else
      {
        px[i-1] = Pt.X(); py[i-1] = Pt.Y(); pz[i-1] = Pt.Z();
      }
    }
  }
  return 0;
}

// ============================================================================
// Mesh an edge with [deflection] on face
// ============================================================================
E_Int __meshEdgeByFace4(const TopoDS_Edge& E, const TopoDS_Face& F, 
                        E_Float hausd,
                        E_Int& nbPoints, K_FLD::FldArrayF& coords,
                        E_Bool reverse)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  //Standard_Real aFirst=geomAdap.FirstParameter(), aEnd=geomAdap.LastParameter();
  Standard_Real pFirst=aFirst, pEnd=aEnd;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);
  Handle(Geom_Surface) surf = BRep_Tool::Surface(F);
  Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
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
  GCPnts_UniformDeflection param(geomAdap, hausd, aFirst, aEnd);
  nbPoints = param.NbPoints();
  printf("hausd face: nbPoints=" SF_D_ " coords=" SF_D_ "\n", nbPoints, coords.getSize()); fflush(stdout);

  if (nbPoints != coords.getSize()) exit(0);

  {
    gp_Pnt Pt; gp_Pnt2d Puv; E_Float u;
    for (E_Int i = 1; i <= nbPoints; i++)
    {
      u = param.Parameter(i);
      C0.D0(u, Pt); pCurve->D0(u, Puv);
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
// Return the nbPoints for meshing E regular with hmax
// can modify hmax
// ============================================================================
E_Int __getNbPts1(const TopoDS_Edge& E, E_Float& hmax, E_Int& nbPoints)
{
  if (BRep_Tool::Degenerated(E)) { nbPoints=2; return 1; }
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  E_Float u0 = geomAdap.FirstParameter();
  E_Float u1 = geomAdap.LastParameter();
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  nbPoints = (E_Int)round(L / hmax);
  nbPoints = std::max(nbPoints, E_Int(2));
  if (nbPoints > MAXNBPTSPEREDGE)
  {
    nbPoints = MAXNBPTSPEREDGE; // hard limit
    hmax = L / nbPoints;
    printf("WARNING: limit edge nbPoints\n");
  } 
  printf("L=%f hmax=%f nbPoints=" SF_D_ "\n", L, hmax, nbPoints);
  return 0;
}

// Return the length of edge
E_Float __getLength(const TopoDS_Edge& E)
{
  if (BRep_Tool::Degenerated(E)) return 0.;
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  E_Float u0 = geomAdap.FirstParameter();
  E_Float u1 = geomAdap.LastParameter();
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  return L;
}

// ============================================================================
// Return the nbPoints for meshing E with given deflection
// can change hausd
// ============================================================================
E_Int __getNbPts2(const TopoDS_Edge& E, E_Float hausd, E_Int& nbPoints)
{
  if (BRep_Tool::Degenerated(E)) { nbPoints=2; return 1; }
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  E_Float u0 = geomAdap.FirstParameter();
  E_Float u1 = geomAdap.LastParameter();
  GCPnts_UniformDeflection param(geomAdap, hausd, u0, u1);
  nbPoints = param.NbPoints();

  if (nbPoints > MAXNBPTSPEREDGE)
  {
    hausd = hausd*(MAXNBPTSPEREDGE / nbPoints); // suppose proportionnel
    GCPnts_UniformDeflection param2(geomAdap, hausd, u0, u1);
    nbPoints = param2.NbPoints();
    printf("WARNING: limit edge nbPoints=" SF_D_ "\n", nbPoints);
  }
  printf("getnpts2: hausd=%f nbPoints=" SF_D_ "\n", hausd, nbPoints); fflush(stdout);
  return 0;
}

// ============================================================================
/* Mesh global edges of CAD, equal distance hmax, return STRUCT */
// ============================================================================
PyObject* K_OCC::meshGlobalEdges1(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float hmax;
  if (!PYPARSETUPLE_(args, O_ R_, &hook, &hmax)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  E_Int nbPoints;
  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
  for (E_Int i=1; i <= edges.Extent(); i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(edges(i));
    __getNbPts1(E, hmax, nbPoints);
  
    // create array
    PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    __meshEdge1(E, nbPoints, *f, false);
    RELEASESHAREDS(o, f);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh global edges of CAD, regular param, given number of points, return STRUCT */
// ============================================================================
PyObject* K_OCC::meshGlobalEdges2(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int nbPoints;
  if (!PYPARSETUPLE_(args, O_ I_, &hook, &nbPoints)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  for (E_Int i=1; i <= edges.Extent(); i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(edges(i));
    
    // create array
    PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    __meshEdge1(E, nbPoints, *f, false);
    RELEASESHAREDS(o, f);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh global edges of CAD, given edges ue, return STRUCT */
// ============================================================================
PyObject* K_OCC::meshGlobalEdges3(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* globalEdges;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &globalEdges)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  for (E_Int i=1; i <= edges.Extent(); i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(edges(i));
    
    PyObject* ge = PyList_GetItem(globalEdges, i-1);
    E_Int ni, nj, nk;
    K_FLD::FldArrayF* fe; K_FLD::FldArrayI* ce;
    char* varString; char* eltType;
    K_ARRAY::getFromArray3(ge, varString, fe, ni, nj, nk, ce, eltType);
    E_Float* ue = fe->begin(4);

    // create array
    PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", ni, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    __meshEdge3(E, ni, ue, *f, false);
    RELEASESHAREDS(o, f);
    RELEASESHAREDS(ge, fe);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh global edges of CAD given deflection, return STRUCT */
// ============================================================================
PyObject* K_OCC::meshGlobalEdges4(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float hausd;
  if (!PYPARSETUPLE_(args, O_ R_, &hook, &hausd)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  E_Int nbPoints=1;
  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
  for (E_Int i=1; i <= edges.Extent(); i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(edges(i));
    __getNbPts2(E, hausd, nbPoints);
    
    // create array
    PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    __meshEdge4(E, hausd, -1., nbPoints, *f, false);
    RELEASESHAREDS(o, f);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh and param global edges of CAD, 
   given number of points or nb of points by hmax or by deflection on faces
   return STRUCT + uv Face */
// ============================================================================
PyObject* K_OCC::meshEdgesByFace(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int nbPoints; E_Int noFace; E_Float hmax; E_Float hausd;
  if (!PYPARSETUPLE_(args, O_ II_ RR_, &hook, &noFace, &nbPoints, &hmax, &hausd)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopExp_Explorer expl;

  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));

  for (expl.Init(surfaces(noFace), TopAbs_EDGE); expl.More(); expl.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(expl.Current());

    // get number of points
    if (hausd > 0) __getNbPts2(E, hausd, nbPoints);
    else if (hmax > 0) __getNbPts1(E, hmax, nbPoints);
    
    // create array
    PyObject* o = K_ARRAY::buildArray3(5, "x,y,z,u,v", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    if (hausd > 0) __meshEdgeByFace4(E, F, hausd, nbPoints, *f, false);
    else __meshEdgeByFace2(E, F, nbPoints, *f, false);
    RELEASESHAREDS(o, f);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh and param global edges of CAD, 
   ue given in global edges list on faces
   return STRUCT + uv Face */
// ============================================================================
PyObject* K_OCC::meshEdgesByFace2(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int noFace;
  PyObject* globalEdges;
  if (!PYPARSETUPLE_(args, O_ I_ O_, &hook, &noFace, &globalEdges)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  PyObject* out = PyList_New(0);
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
  PyObject* ge;

  for (expl.Init(surfaces(noFace), TopAbs_EDGE); expl.More(); expl.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(expl.Current());

    // find no of corresponding global edge
    E_Int noedge = 0;
    for (noedge=1; noedge <= edges.Extent(); noedge++)
    {
        const TopoDS_Edge& El = TopoDS::Edge(edges(noedge));
        //if (El.TShape() == E.TShape()) break;
        if (El.IsSame(E)) break;
        
    }
    ge = PyList_GetItem(globalEdges, noedge-1);
    E_Int ni, nj, nk;
    K_FLD::FldArrayF* fe; K_FLD::FldArrayI* ce;
    char* varString; char* eltType;
    K_ARRAY::getFromArray3(ge, varString, fe, ni, nj, nk, ce, eltType);
    E_Float* ue = fe->begin(4); 

    // create array
    PyObject* o = K_ARRAY::buildArray3(5, "x,y,z,u,v", ni, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray3(o, f);

    // fill array
    __meshEdgeByFace3(E, F, ni, ue, *f, false);
    RELEASESHAREDS(o, f);
    RELEASESHAREDS(ge, fe);
    PyList_Append(out, o); Py_DECREF(o);
  }
  return out;
}

// ============================================================================
/* Mesh edge by face and by wire */
// ============================================================================
PyObject* K_OCC::meshEdgesByFace3(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int noFace;
  E_Float hmax, hausd;
  if (!PYPARSETUPLE_(args, O_ I_ RR_, &hook, &noFace, &hmax, &hausd)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
  // outer wire
  const TopoDS_Wire& OW = ShapeAnalysis::OuterWire(F);
  
  TopAbs_Orientation forientation = F.Orientation();
  if (forientation == TopAbs_FORWARD) printf("face orientation=forward\n");
  else if (forientation == TopAbs_REVERSED) printf("face orientation=reversed\n");

  PyObject* out = PyList_New(0); // sortie par wire
  
  for (expl.Init(surfaces(noFace), TopAbs_WIRE); expl.More(); expl.Next())
  {
    PyObject* le = PyList_New(0); // sortie par edges
    const TopoDS_Wire& W = TopoDS::Wire(expl.Current());
    printf("getting a wire\n");
    //BRepCheck_Wire check = BRepCheck_Wire(W);
    BRepCheck_Wire check(W);
    
    BRepCheck_Status status = check.Orientation(F);
    if (status == BRepCheck_NoError) printf("wire correctly oriented\n");
    else printf("WIRE BADLY oriented\n");
    BRepCheck_Status wclosed = check.Closed();
    if (wclosed == BRepCheck_NoError) printf("wire correctly closed\n");
    else printf("WIRE BADLY closed\n");
    
    //E_Bool isOuter = false;
    //if (W == OW) { printf("is outer wire\n"); isOuter = true; }
    //status = check.CheckOuterBound();
    //if (status == BRepCheck_NoError) printf("is outer (test2)\n");
    //else printf("is inner (test2)\n");

    E_Float surface = ShapeAnalysis::ContourArea(W);
    printf("wire surface=%f\n", surface);	

    // Compute oriented surface to orient externe = CCW, intern = CW
    BRepBuilderAPI_MakeFace faceMaker(F, W);
    faceMaker.Build();
    const TopoDS_Face& F2 = faceMaker.Face();
    GProp_GProps SProps;
    BRepGProp::SurfaceProperties(F2, SProps);
    Standard_Real area = SProps.Mass();
    gp_Pnt center = SProps.CentreOfMass();
    printf("Signed surface = %f\n", area);
    printf("Center of mass = %f %f %f\n", center.X(), center.Y(), center.Z());

    TopAbs_Orientation worientation = W.Orientation();
    if (worientation == TopAbs_FORWARD) printf("wire orientation=forward\n");
    else if (worientation == TopAbs_REVERSED) printf("wire orientation=reversed\n");
    else if (worientation == TopAbs_INTERNAL) printf("wire orientation=internal\n");
    else if (worientation == TopAbs_EXTERNAL) printf("wire orientation=external\n");
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
      printf("getting edges of wire\n");
      
      const TopoDS_Edge& E = TopoDS::Edge(expl2.Current());  
      
      TopAbs_Orientation eorientation = E.Orientation();
      if (eorientation == TopAbs_FORWARD) printf("edge orientation=forward\n");
      else if (eorientation == TopAbs_REVERSED) printf("edge orientation=reversed\n");
      else if (eorientation == TopAbs_INTERNAL) printf("edge orientation=internal\n");
      else if (eorientation == TopAbs_EXTERNAL) printf("edge orientation=external\n");

      E_Int nbPoints = (E_Int)41; // brute force is BAD!!
      E_Float hmaxe = hmax;
      E_Float hausde = hausd;
      printf("hmax=%f hausd=%f\n", hmax, hausd);
      if (hausde < 0 && hmaxe > 0) __getNbPts1(E, hmaxe, nbPoints); // nbre de points correpondant a hmaxw (possiblement limite)
      else if (hausde > 0 && hmaxe < 0) __getNbPts2(E, hausde, nbPoints);
      else 
      {
        printf("hausd and hmax not possible for now.\n");
        return NULL;
      }

      // create array
      PyObject* o = K_ARRAY::buildArray3(5, "x,y,z,u,v", nbPoints, 1, 1, 1);
      FldArrayF* f; K_ARRAY::getFromArray3(o, f);
      discreteWire.push_back(f);
      discreteWire2.push_back(o);

      E_Bool reversed = false;
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
      }
      if (hausde < 0 && hmaxe > 0)
        __meshEdgeByFace1(E, F, nbPoints, *f, reversed);
      else if (hausde > 0 && hmaxe < 0)
        __meshEdgeByFace4(E, F, hausd, nbPoints, *f, reversed);

      //RELEASESHAREDS(o, f); // done later
      if (forientation == TopAbs_FORWARD)
      {
        if (worientation == TopAbs_FORWARD) PyList_Append(le, o); 
        else PyList_Insert(le, 0, o); 
      }
      else
      {
        if (worientation == TopAbs_FORWARD) PyList_Insert(le, 0, o); 
        else PyList_Append(le, o); 
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
      if (worientation == TopAbs_FORWARD)
      {
        std::reverse(discreteWire.begin(), discreteWire.end());  
        std::reverse(discreteWire2.begin(), discreteWire2.end());  
      }
    }

    size_t nedges = discreteWire.size();
    if (wclosed == BRepCheck_NoError)
    {
      printf("closing the discrete wire\n");
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
  fflush(stdout);

  return out;
}

//========================================
/* Get global edges number list by face */
//========================================
PyObject* K_OCC::getEdgeNoByFace(PyObject* self, PyObject* args)
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
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  PyObject* out = PyList_New(0);

  const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
  TopAbs_Orientation forientation = F.Orientation();
  
  // nouvelle version par wire pour etre coherent avec meshEdgesByFace
  for (expl.Init(surfaces(noFace), TopAbs_WIRE); expl.More(); expl.Next())
  {
    const TopoDS_Wire& W = TopoDS::Wire(expl.Current());
    BRepTools_WireExplorer expl2;
    for (expl2.Init(W, F); expl2.More(); expl2.Next())
    {
      const TopoDS_Edge& E = TopoDS::Edge(expl2.Current());
      //TopAbs_Orientation eorientation = E.Orientation();
      /* find edge number in global edge list */
      E_Int i = 0;
      for (i=1; i <= edges.Extent(); i++)
      {
        const TopoDS_Edge& El = TopoDS::Edge(edges(i));
        //if (El.TShape() == E.TShape()) break;
        if (El.IsSame(E)) break;
      }

      if (forientation == TopAbs_FORWARD)
      {
        PyList_Append(out, Py_BuildValue("l", i));
      }
      else
      {
        PyList_Insert(out, 0, Py_BuildValue("l", i));
      }
    }
  }

  return out;
}
