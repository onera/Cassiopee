/*    
    Copyright 2013-2025 Onera.

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

#include "Nuga/include/DelaunayMath.h"
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
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "Geom_Surface.hxx"
#include "TopExp.hxx"
#include "GeomLProp_SLProps.hxx"
#include <Precision.hxx>

// ultimate (best) functions

// ============================================================================
// Return uniform (h constant) distribution of NbPoints on edge
// ============================================================================
E_Int __getUniform(const TopoDS_Edge& E, E_Int nbPoints, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  GCPnts_UniformAbscissa param(geomAdap, int(nbPoints), u0, u1);
  ue = new E_Float [nbPoints];
  for (E_Int i = 0; i < nbPoints; i++) ue[i] = param.Parameter(i+1);
  return 0;
}

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
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  nbPoints = (E_Int)round(L / hmax);
  nbPoints = std::max(nbPoints, E_Int(3));
  GCPnts_UniformAbscissa param(geomAdap, int(nbPoints), u0, u1);
  ue = new E_Float [nbPoints];
  try // that fails some times
  {
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = param.Parameter(i+1);
  }
  catch (const Standard_Failure& theErr)
  {
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = i*(u1-u0)/(nbPoints-1)+u0;
    printf("Warning: regular param used on edge.\n");
  }
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
    nbPoints = 2; 
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
// Geom distrib entre u0 et u1, h0 et h1 (interieurs)
void geom1(E_Float u0, E_Float u1, E_Float h0, E_Float h1, E_Int& N, E_Float*& ue)
{
  E_Float r = (u1-u0-h0)/(u1-u0-h1);
  //printf("r=%f\n", r);
  E_Float a = log(r);
  if (a > 1.e-12) // r!=1
  {
    N = round(log(h1/h0)/a)+2;
    if (N < 2) N = 2;
    h0 = (u1-u0)*(1.-r)/(1.-pow(r, N-1));
    ue = new E_Float [N];
    ue[0] = u0;
    for (E_Int i = 1; i < N; i++) ue[i] = u0 + h0*(1.-pow(r,i))/(1.-r);
  }
  else // r=1
  {
    N = round((u1-u0)/h0+1);
    if (N < 2) N = 2;
    h0 = (u1-u0)/(N-1);
    ue = new E_Float [N];
    for (E_Int i = 0; i < N; i++) ue[i] = u0 + i*h0;
  }
  //for (E_Int i = 1; i < N; i++) ue[i] = ue[i-1] + h0*pow(r, i-1);  
  //ue[N-1] = u1; // force
  //for (E_Int i = 0; i < N; i++) printf("%g\n", ue[i]);
  //printf("h0=%f real=%f\n", h0, ue[1]-ue[0]);
  //printf("h1=%f real=%f\n", h1, ue[N-1]-ue[N-2]);
}

// ============================================================================
// Geom distrib entre u0 et u1, h0 et h1/r (interieurs)
void geom2(E_Float u0, E_Float u1, E_Float h0, E_Float h1, E_Int& N, E_Float*& ue)
{
  E_Float r = (u1-u0+h1-h0)/(u1-u0);
  //printf("r=%f\n", r);
  E_Float a = log(r);
  if (a > 1.e-12) // r!=1
  {
    N = round(log(h1/h0)/a)+2;
    if (N < 2) N = 2;
    h0 = (u1-u0)*(1.-r)/(1.-pow(r, N-1));
    ue = new E_Float [N];
    ue[0] = u0;
    for (E_Int i = 1; i < N; i++) ue[i] = u0 + h0*(1.-pow(r,i))/(1.-r);
  }
  else // r=1
  {
    N = round((u1-u0)/h0+1);
    if (N < 2) N = 2;
    h0 = (u1-u0)/(N-1);
    ue = new E_Float [N];
    for (E_Int i = 0; i < N; i++) ue[i] = u0 + i*h0;
  }
  //for (E_Int i = 1; i < N; i++) ue[i] = ue[i-1] + h0*pow(r, i-1);  
  //ue[N-1] = u1; // force
  //for (E_Int i = 0; i < N; i++) printf("%d %f\n", i, ue[i]);
  //printf("h0=%f real=%f\n", h0, ue[1]-ue[0]);
  //printf("h1/r=%f real=%f\n", h1/r, ue[N-1]-ue[N-2]);
}

// ============================================================================
// Geom distrib entre u0 et u1, h0/r et h1/r (interieurs)
void geom3(E_Float u0, E_Float u1, E_Float h0, E_Float h1, E_Int& N, E_Float*& ue)
{
  E_Float delta = (u1-u0-h0)*(u1-u0-h0)+4*(u1-u0)*h1;
  if (delta >= 0) delta = sqrt(delta);

  E_Float r = ((u1-u0-h0)+delta)/(2*(u1-u0));
  //E_Float r1 = ((u1-u0-h0)-delta)/(2*(u1-u0));
  //printf("r=%f\n", r);
  E_Float a = log(r);
  if (a > 1.e-12) // r!=1
  {
    N = round(log(h1/h0)/a)+2;
    if (N < 2) N = 2;
    h0 = (u1-u0)*(1.-r)/(1.-pow(r, N-1));
    ue = new E_Float [N];
    ue[0] = u0;
    for (E_Int i = 1; i < N; i++) ue[i] = u0 + h0*(1.-pow(r,i))/(1.-r);
  }
  else // r=1
  {
    N = round((u1-u0)/h0+1);
    if (N < 2) N = 2;
    h0 = (u1-u0)/(N-1);
    ue = new E_Float [N];
    for (E_Int i = 0; i < N; i++) ue[i] = u0 + i*h0;
  }
  //for (E_Int i = 1; i < N; i++) ue[i] = ue[i-1] + h0*pow(r, i-1);  
  //ue[N-1] = u1; // force
  for (E_Int i = 0; i < N; i++) printf("%d %f\n", i, ue[i]);
  printf("h0/r=%f real=%f\n", h0/r, ue[1]-ue[0]);
  printf("h1/r=%f real=%f\n", h1/r, ue[N-1]-ue[N-2]);
}

// ============================================================================
// Geom distrib entre u0 et u1, h0/r et h1 (interieurs)
void geom4(E_Float u0, E_Float u1, E_Float h0, E_Float h1, E_Int& N, E_Float*& ue)
{
  E_Float r = (u1-u0+h1-h0)/(u1-u0);
  //printf("r=%f\n", r);
  E_Float a = log(r);
  if (a > 1.e-12) // r!=1
  {
    N = round(log(h1/h0)/a)+2;
    if (N < 2) N = 2;
    h0 = (u1-u0)*(1.-r)/(1.-pow(r, N-1));
    ue = new E_Float [N];
    ue[0] = u0;
    for (E_Int i = 1; i < N; i++) ue[i] = u0 + h0*(1.-pow(r,i))/(1.-r);
  }
  else // r=1
  {
    N = round((u1-u0)/h0+1);
    if (N < 2) N = 2;
    h0 = (u1-u0)/(N-1);
    ue = new E_Float [N];
    for (E_Int i = 0; i < N; i++) ue[i] = u0 + i*h0;
  }
  //for (E_Int i = 1; i < N; i++) ue[i] = ue[i-1] + h0*pow(r, i-1);  
  //ue[N-1] = u1; // force
  //for (E_Int i = 0; i < N; i++) printf("%d %f\n", i, ue[i]);
  printf("h0/r=%f real=%f\n", h0/r, ue[1]-ue[0]);
  printf("h1=%f real=%f\n", h1, ue[N-1]-ue[N-2]);
}

// ============================================================================
// Return the nbPoints and ue for meshing E with hmin/hmax/hausd evaluated
// only on edge using geometric progression between points
// ============================================================================
E_Int __getParamHminHmaxHausdE(const TopoDS_Edge& E, E_Float hmin, E_Float hmax, E_Float hausd, E_Int& nbPoints, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve());
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  L = K_FUNC::E_max(L, 1.e-12);
  
  // Compute local h
  Handle(Geom_Curve) curve;
  Standard_Real first, last;
  curve = BRep_Tool::Curve(E, first, last);
  gp_Pnt point;
  gp_Vec dE, d2E;
  E_Float U, dE1, dE2, dE3, d2E1, d2E2, d2E3;
  E_Float t1, t2, t3, a, b, rho;

  E_Int npts = 2;
  std::vector<E_Float> Us(npts);
  for (E_Int i = 0; i < npts; i++) Us[i] = u0+(i/(npts-1.))*(u1-u0);
  std::vector<E_Float> h(npts);

  for (E_Int i = 0; i < npts; i++)
  {
    U = Us[i];
    // Compute the derivatives
    curve->D2(U, point, dE, d2E);
  
    // compute curvature
    dE1 = dE.X();
    dE2 = dE.Y();
    dE3 = dE.Z();
    d2E1 = d2E.X();
    d2E2 = d2E.Y();
    d2E3 = d2E.Z();
    t1 = dE2*d2E3-dE3*d2E2;
    t2 = dE3*d2E1-dE1*d2E3;
    t3 = dE1*d2E2-dE2*d2E1;
  
    a = t1*t1+t2*t2+t3*t3;
    a = std::sqrt(a);
    b = d2E1*d2E1+d2E2*d2E2+d2E3*d2E3;
    b = std::pow(b, 3./2.);
    if (a > 1.e-12) rho = b / a;
    else rho = 1.e12;

    // Compute local h
    h[i] = std::sqrt(1.*hausd*rho);
    h[i] = K_FUNC::E_min(h[i], hmax);
    h[i] = K_FUNC::E_max(h[i], hmin);
    //printf("rho=%g h=%g, h=%g\n", rho, std::sqrt(8.*hausd*rho), h[i]);
    h[i] = h[i]*(u1-u0)/L;
  }

  E_Int N;
  std::vector<E_Float*> uel(npts);
  std::vector<E_Int> Np(npts);
  
  E_Int Ntot = 0;
  for (E_Int i = 0; i < npts-1; i++)
  {
    // geom2
    geom2(Us[i], Us[i+1], h[i], h[i+1], N, uel[i]);
    Np[i] = N;
    Ntot += N-1;
  }
  Ntot += 1;
  geom1(Us[npts-2], Us[npts-1], h[npts-2], h[npts-1], N, uel[npts-1]);

  // reassemble
  ue = new E_Float [Ntot];
  nbPoints = Ntot;
  N = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    E_Float* uep = uel[i];
    for (E_Int j = 0; j < Np[i]-1; j++) ue[j+N] = uep[j];
    N += Np[i]-1;
  }
  ue[N+1] = u1;
  
  //printf("final: %g %g \n", u0, u1);
  //for (E_Int i = 0; i < Ntot; i++) printf("%g ", ue[i]);
  //printf("\n");

  // verifie que la suite est bien croissante
  for (E_Int i = 0; i < Ntot-1; i++) 
  {
    if (ue[i+1] < ue[i]) printf("warning: %g ", ue[i]);
  }

  for (E_Int i = 0; i < npts; i++)
  {
    delete [] uel[i];
  }
  return 1;
}

// ============================================================================
// Return the nbPoints and ue for meshing E with hmin/hmax/hausd evaluated
// only on edge using progressive walk
// ============================================================================
E_Int __getParamHminHmaxHausdE4(const TopoDS_Edge& E, E_Float hmin, E_Float hmax, E_Float hausd, E_Int& nbPoints, E_Float*& ue)
{
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve());
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  L = K_FUNC::E_max(L, 1.e-12);

  gp_Pnt point;
  gp_Vec dE, d2E;
  E_Float ues, dE1, dE2, dE3, d2E1, d2E2, d2E3;
  E_Float h, t1, t2, t3, a, b, rho;

  Handle(Geom_Curve) curve;
  Standard_Real first, last;
  curve = BRep_Tool::Curve(E, first, last);

  E_Float Ltot = 0.;
  E_Int N = 0;
  E_Int sizeCont = 5000;
  E_Float* hi = new E_Float [sizeCont]; // a dimensionner dynamiquement
  ues = u0; h = 0.;

  while (Ltot < L)
  {
    GCPnts_AbscissaPoint Pt(1.e-10, geomAdap, h, ues);
    ues = Pt.Parameter();

    // Compute local h
    // Compute the derivatives
    curve->D2(ues, point, dE, d2E);
  
    // compute curvature
    dE1 = dE.X();
    dE2 = dE.Y();
    dE3 = dE.Z();
    d2E1 = d2E.X();
    d2E2 = d2E.Y();
    d2E3 = d2E.Z();
    t1 = dE2*d2E3-dE3*d2E2;
    t2 = dE3*d2E1-dE1*d2E3;
    t3 = dE1*d2E2-dE2*d2E1;
  
    a = t1*t1+t2*t2+t3*t3;
    a = std::sqrt(a);
    b = d2E1*d2E1+d2E2*d2E2+d2E3*d2E3;
    b = std::pow(b, 3./2.);
    if (a > 1.e-12) rho = b / a;
    else rho = 1.e12;

    // Compute local h
    h = std::sqrt(8.*hausd*rho);
    h = K_FUNC::E_min(h, hmax);
    h = K_FUNC::E_max(h, hmin);
    hi[N] = h;
    Ltot += h; N += 1;
    if (N >= sizeCont)
    {
      sizeCont += 5000;
      E_Float* hi2 = new E_Float [sizeCont];
      for (E_Int i = 0; i < N; i++) hi2[i] = hi[i];
      delete [] hi;
      hi = hi2;
    }
  }

  // smooth hi
  /*
  for (E_Int i = 0; i < N-1; i++)
  {
    hi[i] = 0.5*(hi[i+1]+hi[i]);
  }
  hi[N-1] = 0.5*(hi[N-2]+hi[N-1]);
  */

  E_Float r = Ltot-L;
  if (r > 0.5*hi[N-1] && N >= 3) { N = N-1; Ltot = Ltot - hi[N]; }
  
  //printf("hi = ");
  //for (E_Int i = 0; i < N; i++) printf("%g ", hi[i]);
  //printf("\n");

  // scale hi to match full L
  E_Float s = L-Ltot;
  s = s/N; 
  for (E_Int i = 0; i < N; i++) hi[i] = hi[i]+s;
  
  // DBX
  //Ltot = 0.;
  //for (E_Int i = 0; i < N; i++) Ltot += hi[i];
  //printf("Ltot=%g %g\n", Ltot, L);
  // ENDDBX

  // Get ue
  ue = new E_Float [N+1];
  ues = u0;
  ue[0] = u0;
  for (E_Int i = 0; i < N; i++)
  {
    GCPnts_AbscissaPoint Pt(1.e-10, geomAdap, hi[i], ues);
    ues = Pt.Parameter();
    ue[i+1] = ues;
  }
  ue[N] = u1; // forced

  //printf("N= %d\n", N);
  //printf("ue5= ");
  //for (E_Int i = 0; i < N; i++) printf("%g ", ue[i]);
  //printf("\n");

  delete [] hi;
  nbPoints = N+1;

  return 1;
}

// ===============================================================================
// Return the nbPoints and ue for meshing E with hmin/hmax/hausd evaluated  
// as max of faces using geometric progression between points
// ===============================================================================
E_Int __getParamHminHmaxHausdF(const TopoDS_Edge& E, E_Float hmin, E_Float hmax, E_Float hausd, E_Int& nbPoints, E_Float*& ue,
  TopoDS_Shape* shape)
{
  TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
  TopExp::MapShapesAndAncestors(*shape, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve());
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  //printf("===============edge =========================\n");
  //printf("u0=%g u1=%g\n", u0, u1);
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  L = K_FUNC::E_max(L, 1.e-12);

  gp_Pnt2d Puv;
  E_Float u,v;
  E_Float ues; 
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  //printf("afirst=%g aend=%g\n", aFirst, aEnd);
  Standard_Real pFirst = aFirst, pEnd=aEnd;
  
  gp_Pnt point;
  gp_Vec DU1, DV1, DU2, DV2, DUV;
  E_Float rho, hh, K1, K2;

  E_Int npts = 2;
  std::vector<E_Float> Us(npts);
  for (E_Int i = 0; i < npts; i++) Us[i] = i/(npts-1.);
    
  std::vector<E_Float> h(npts);
  for (E_Int i = 0; i < npts; i++) h[i] = K_CONST::E_MAX_FLOAT;
  
  // Find faces of E
  const TopTools_ListOfShape& connectedFaces = edgeFaceMap.FindFromKey(E);
  for (TopTools_ListIteratorOfListOfShape it(connectedFaces); it.More(); it.Next()) 
  {
    const TopoDS_Face& F = TopoDS::Face(it.Value());
    Handle(Geom_Surface) surface = BRep_Tool::Surface(F);
    
    Standard_Real BU1, BU2, BV1, BV2; surface->Bounds(BU1, BU2, BV1, BV2);
    //printf("face bound = %g %g and %g %g\n", BU1, BU2, BV1, BV2);

    // p curve on F
    Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
    //printf("pfirst=%g pend=%g\n", pFirst, pEnd);
    //printf("pCurve = %p\n", (void*)pCurve);
    // estimate at some points
    for (E_Int i = 0; i < npts; i++)
    {
      ues = Us[i]*(pEnd-pFirst)+pFirst;

      pCurve->D0(ues, Puv);
      u = Puv.X(); v = Puv.Y(); // u,v edge on surface
      if (u <= BU1) u = BU1+1.e-1;
      if (u >= BU2) u = BU2-1.e-1;
      if (v <= BV1) v = BV1+1.e-1;
      if (v >= BV2) v = BV2-1.e-1;
      
// must be MINE (0) or OCC (1)
#define CURVATURE 1
#if CURVATURE == 0
      E_Float nx, ny, nz, n2;
      E_Float G11, G12, G21, G22, detG;
      E_Float M11, M12, M21, M22, detM;
      E_Float det, S11, S12, S21, S22, Mi11, Mi21, Mi12, Mi22, detS;

      // estimer les derivees sur la surface
      surface->D2(u, v, point, DU1, DV1, DU2, DV2, DUV);
      //printf("point=%g %g %g\n", point.X(), point.Y(), point.Z());
      //printf("dU1=%g %g %g\n", DU1.X(), DU1.Y(), DU1.Z());
      //printf("dV1=%g %g %g\n", DV1.X(), DV1.Y(), DV1.Z());
      
      // get normal
      nx = DU1.Y()*DV1.Z()-DU1.Z()*DV1.Y();
      ny = DU1.Z()*DV1.X()-DU1.X()*DV1.Z();
      nz = DU1.X()*DV1.Y()-DU1.Y()*DV1.X();
      n2 = nx*nx+ny*ny+nz*nz;
      n2 = std::sqrt(n2);
      if (n2 < 1.e-24) n2 = 1.e12;
      else n2 = 1./n2;
      nx = nx*n2;
      ny = ny*n2;
      nz = nz*n2;
      //printf("normale= %g %g %g norm=%g\n", nx, ny, nz, nx*nx+ny*ny+nz*nz);

      // Metric
      M11 = DU1.X()*DU1.X()+DU1.Y()*DU1.Y()+DU1.Z()*DU1.Z();
      M22 = DV1.X()*DV1.X()+DV1.Y()*DV1.Y()+DV1.Z()*DV1.Z();
      M12 = M21 = DU1.X()*DV1.X()+DU1.Y()*DV1.Y()+DU1.Z()*DV1.Z();
      //printf("matrice M:\n");
      //printf("  %g %g\n", M11, M12);
      //printf("  %g %g\n", M21, M22);

      detM = M11*M22-M12*M21;
      if (std::abs(detM) < 1.e-12) det = 1.e12; // invalid metric
      else det = 1./detM;
      Mi11 = det*M22;
      Mi21 = -det*M21;
      Mi12 = -det*M12;
      Mi22 = det*M11;
      detM = std::abs(detM);
      //printf("detM=%g\n", detM);

      // matrice de courbure
      G11 = DU2.X()*nx+DU2.Y()*ny+DU2.Z()*nz;
      G22 = DV2.X()*nx+DV2.Y()*ny+DV2.Z()*nz;
      G12 = G21 = DUV.X()*nx+DUV.Y()*ny+DUV.Z()*nz;
      detG = G11*G22-G12*G21;
      detG = std::abs(detG);
      
      S11 = G11*Mi11 + G12*Mi21;
      S12 = G11*Mi12 + G12*Mi22;
      S21 = G21*Mi11 + G22*Mi21;
      S22 = G21*Mi12 + G22*Mi22;
      
      detS = S11*S22-S12*S21;
      //printf("dU2=%g %g %g\n", DU2.X(), DU2.Y(), DU2.Z());
      //printf("dV2=%g %g %g\n", DV2.X(), DV2.Y(), DV2.Z());
      //printf("dUV=%g %g %g\n", DUV.X(), DUV.Y(), DUV.Z());
      //printf("detG=%g\n", detG);

      // courbure de gauss - rayon K1*K2
      //if (detS > 1.e-12) rho = 1./detS;
      //else rho = 1.e12;
      //printf("gauss=%g ou bien %g\n", detG/detM, detS);

      // K1 et K2
      //printf("matrice G:\n");
      //printf("  %g %g\n", G11, G12);
      //printf("  %g %g\n", G21, G22);
      
      K_LINEAR::DelaunayMath sl;
      sl.eigen_values(S11, S22, S12, K2, K1);
      //printf("K1=%g K2=%g K1K2=%g\n", K1, K2, K1*K2);
      //printf("mean=%g\n", 0.5*(S11+S22));

      rho = K_FUNC::E_max(std::abs(K1), std::abs(K2));
      if (rho > 1.e-12) rho = 1./rho;
      else rho = 1.e12;

#else
      // with open cascade
      GeomLProp_SLProps props(surface, u, v, 2, Precision::Confusion());

      if (props.IsCurvatureDefined()) 
      {
        Standard_Real gaussianCurvature = props.GaussianCurvature();
        Standard_Real meanCurvature = props.MeanCurvature();

        //std::cout << "OCC: gauss: " << gaussianCurvature << std::endl;
        //std::cout << "OCC: mean: " << meanCurvature << std::endl;

        // Principal curvatures
        K1 = meanCurvature + sqrt(meanCurvature * meanCurvature - gaussianCurvature);
        K2 = meanCurvature - sqrt(meanCurvature * meanCurvature - gaussianCurvature);
        //std::cout << "OCC: k1: " << K1 << " k2: " <<K2 << std::endl;
        rho = K_FUNC::E_max(std::abs(K1), std::abs(K2));
        if (rho > 1.e-12) rho = 1./rho;
        else rho = 1.e12;
      } 
      else 
      {
        rho = 1.e12;
        //std::cout << "OCC: Curvature is not defined at the given parameters." << std::endl;
      }
#endif

      hh = std::sqrt(8.*hausd*rho);
      hh = K_FUNC::E_min(hh, hmax);
      hh = K_FUNC::E_max(hh, hmin);
      //printf("rho=%g h=%g, h=%g\n", rho, std::sqrt(1.*hausd*rho), hh);
      h[i] = K_FUNC::E_min(h[i], hh);
    }
  }

  //printf("hi: ");
  //for (E_Int i = 0; i < npts; i++) 
  //{
  //  printf("%g ", h[i]);
  //}
  //printf("\n");

  for (E_Int i = 0; i < npts; i++)
  {
    // passage en parametres
    h[i] = h[i]*(u1-u0)/L;
  }

  E_Int N;
  std::vector<E_Float*> uel(npts);
  std::vector<E_Int> Np(npts);

  E_Int Ntot = 0;
  for (E_Int i = 0; i < npts-1; i++)
  {
    // geom2
    geom2(Us[i]*(aEnd-aFirst)+aFirst, Us[i+1]*(aEnd-aFirst)+aFirst, h[i], h[i+1], N, uel[i]);
    Np[i] = N;
    Ntot += N-1;
  }
  Ntot += 1;
  geom1(Us[npts-2]*(aEnd-aFirst)+aFirst, Us[npts-1]*(aEnd-aFirst)+aFirst, h[npts-2], h[npts-1], N, uel[npts-1]);

  // reassemble
  ue = new E_Float [Ntot];
  nbPoints = Ntot;
  N = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    E_Float* uep = uel[i];
    for (E_Int j = 0; j < Np[i]-1; j++) ue[j+N] = uep[j];
    N += Np[i]-1;
  }
  ue[N+1] = u1;

  // verifie que la suite est bien croissante
  for (E_Int i = 0; i < Ntot-1; i++) 
  {
    if (ue[i+1] < ue[i]) printf("warning: %g ", ue[i]);
  }

  //printf("final: ");
  //for (E_Int i = 0; i < Ntot; i++) 
  //{ printf("%g ", ue[i]); }
  //printf("\n");

  for (E_Int i = 0; i < npts; i++)
  {
    delete [] uel[i];
  }

  return 1;

}

// ===============================================================================
// Return the nbPoints and ue for meshing E with hmin/hmax/hausd evaluated 
// as max of faces using progressive walk
// ===============================================================================
E_Int __getParamHminHmaxHausdF5(const TopoDS_Edge& E, E_Float hmin, E_Float hmax, E_Float hausd, E_Int& nbPoints, E_Float*& ue,
  TopoDS_Shape* shape)
{
  TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
  TopExp::MapShapesAndAncestors(*shape, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geomAdap(C0.Curve());
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  //printf("===============edge =========================\n");
  //printf("edge u0=%g u1=%g\n", u0, u1);
  if (BRep_Tool::Degenerated(E))
  { 
    nbPoints = 2;
    ue = new E_Float [nbPoints];
    for (E_Int i = 0; i < nbPoints; i++) ue[i] = u0;
    return 1; 
  }

  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  L = K_FUNC::E_max(L, 1.e-12);

  gp_Pnt2d Puv;
  E_Float u,v, ues;
  Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
  //printf("afirst=%g aend=%g\n", aFirst, aEnd);
  Standard_Real pFirst=aFirst, pEnd=aEnd;
  
  gp_Pnt point;
  gp_Vec DU1, DV1, DU2, DV2, DUV;
  E_Float rho, hh, K1, K2, h;


  E_Float Ltot = 0.;
  E_Int N = 0;
  E_Int sizeCont = 5000;
  E_Float* hi = new E_Float [sizeCont]; // a dimensionner dynamiquement
  ues = u0; h = 0.;

  while (Ltot < L)
  {
    GCPnts_AbscissaPoint Pt(1.e-10, geomAdap, h, ues);
    ues = Pt.Parameter();

    h = K_CONST::E_MAX_FLOAT;
    
    // Find faces of E
    const TopTools_ListOfShape& connectedFaces = edgeFaceMap.FindFromKey(E);
    for (TopTools_ListIteratorOfListOfShape it(connectedFaces); it.More(); it.Next()) 
    {
      const TopoDS_Face& F = TopoDS::Face(it.Value());
      Handle(Geom_Surface) surface = BRep_Tool::Surface(F);
    
      Standard_Real BU1, BU2, BV1, BV2; surface->Bounds(BU1, BU2, BV1, BV2);

      // p curve on F
      Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(E, F, pFirst, pEnd);
      
      //C0.D0(ues, point);
      //printf("param=%g, %g %g %g\n", ues, point.X(), point.Y(), point.Z());

      pCurve->D0(ues, Puv);
      u = Puv.X(); v = Puv.Y(); // u,v edge on surface
      if (u <= BU1) u = BU1+1.e-1;
      if (u >= BU2) u = BU2-1.e-1;
      if (v <= BV1) v = BV1+1.e-1;
      if (v >= BV2) v = BV2-1.e-1;

// must be MINE (0) or OCC (1)
#define CURVATURE 1
#if CURVATURE == 0
      E_Float nx, ny, nz, n2;
      E_Float G11, G12, G21, G22, detG;
      E_Float M11, M12, M21, M22, detM;
      E_Float det, S11, S12, S21, S22, Mi11, Mi21, Mi12, Mi22, detS;

      //printf("surf: u=%g, v=%g\n", u, v);
      // estimer les derivees sur la surface
      surface->D2(u, v, point, DU1, DV1, DU2, DV2, DUV);
      //printf("point=%g %g %g\n", point.X(), point.Y(), point.Z());
      //printf("dU1=%g %g %g\n", DU1.X(), DU1.Y(), DU1.Z());
      //printf("dV1=%g %g %g\n", DV1.X(), DV1.Y(), DV1.Z());
      
      // get normal
      nx = DU1.Y()*DV1.Z()-DU1.Z()*DV1.Y();
      ny = DU1.Z()*DV1.X()-DU1.X()*DV1.Z();
      nz = DU1.X()*DV1.Y()-DU1.Y()*DV1.X();
      n2 = nx*nx+ny*ny+nz*nz;
      n2 = std::sqrt(n2);
      if (n2 < 1.e-24) n2 = 1.e12;
      else n2 = 1./n2;
      nx = nx*n2;
      ny = ny*n2;
      nz = nz*n2;
      //printf("normale= %g %g %g norm=%g\n", nx, ny, nz, nx*nx+ny*ny+nz*nz);

      // Metric
      M11 = DU1.X()*DU1.X()+DU1.Y()*DU1.Y()+DU1.Z()*DU1.Z();
      M22 = DV1.X()*DV1.X()+DV1.Y()*DV1.Y()+DV1.Z()*DV1.Z();
      M12 = M21 = DU1.X()*DV1.X()+DU1.Y()*DV1.Y()+DU1.Z()*DV1.Z();
      //printf("matrice M:\n");
      //printf("  %g %g\n", M11, M12);
      //printf("  %g %g\n", M21, M22);

      detM = M11*M22-M12*M21;
      if (std::abs(detM) < 1.e-12) det = 1.e12; // invalid metric
      else det = 1./detM;
      Mi11 = det*M22;
      Mi21 = -det*M21;
      Mi12 = -det*M12;
      Mi22 = det*M11;
      detM = std::abs(detM);
      //printf("detM=%g\n", detM);

      // matrice de courbure
      G11 = DU2.X()*nx+DU2.Y()*ny+DU2.Z()*nz;
      G22 = DV2.X()*nx+DV2.Y()*ny+DV2.Z()*nz;
      G12 = G21 = DUV.X()*nx+DUV.Y()*ny+DUV.Z()*nz;
      detG = G11*G22-G12*G21;
      detG = std::abs(detG);
      
      S11 = G11*Mi11 + G12*Mi21;
      S12 = G11*Mi12 + G12*Mi22;
      S21 = G21*Mi11 + G22*Mi21;
      S22 = G21*Mi12 + G22*Mi22;
      
      detS = S11*S22-S12*S21;
      //printf("dU2=%g %g %g\n", DU2.X(), DU2.Y(), DU2.Z());
      //printf("dV2=%g %g %g\n", DV2.X(), DV2.Y(), DV2.Z());
      //printf("dUV=%g %g %g\n", DUV.X(), DUV.Y(), DUV.Z());
      //printf("detG=%g\n", detG);

      // courbure de gauss - rayon K1*K2
      //if (detS > 1.e-12) rho = 1./detS;
      //else rho = 1.e12;
      //printf("gauss=%g ou bien %g\n", detG/detM, detS);

      // K1 et K2
      //printf("matrice G:\n");
      //printf("  %g %g\n", G11, G12);
      //printf("  %g %g\n", G21, G22);
      
      K_LINEAR::DelaunayMath sl;
      sl.eigen_values(S11, S22, S12, K2, K1);
      //printf("K1=%g K2=%g K1K2=%g\n", K1, K2, K1*K2);
      //printf("mean=%g\n", 0.5*(S11+S22));

      rho = K_FUNC::E_max(std::abs(K1), std::abs(K2));
      if (rho > 1.e-12) rho = 1./rho;
      else rho = 1.e12;

#else
      // with open cascade
      GeomLProp_SLProps props(surface, u, v, 2, Precision::Confusion());

      if (props.IsCurvatureDefined()) 
      {
        Standard_Real gaussianCurvature = props.GaussianCurvature();
        Standard_Real meanCurvature = props.MeanCurvature();

        //std::cout << "OCC: gauss: " << gaussianCurvature << std::endl;
        //std::cout << "OCC: mean: " << meanCurvature << std::endl;

        // Principal curvatures
        K1 = meanCurvature + sqrt(meanCurvature * meanCurvature - gaussianCurvature);
        K2 = meanCurvature - sqrt(meanCurvature * meanCurvature - gaussianCurvature);
        //std::cout << "OCC: k1: " << K1 << " k2: " <<K2 << std::endl;
        rho = K_FUNC::E_max(std::abs(K1), std::abs(K2));
        if (rho > 1.e-12) rho = 1./rho;
        else rho = 1.e12;
      } 
      else 
      {
        rho = 1.e12;
        //std::cout << "OCC: Curvature is not defined at the given parameters." << std::endl;
      }
#endif

      hh = std::sqrt(8.*hausd*rho);
      hh = K_FUNC::E_min(hh, hmax);
      hh = K_FUNC::E_max(hh, hmin);
      //printf("rho=%g h=%g, h=%g\n", rho, std::sqrt(1.*hausd*rho), hh); fflush(stdout);
      h = K_FUNC::E_min(h, hh);
    }
    hi[N] = h;
    Ltot += h; N += 1;
    if (N >= sizeCont)
    {
      sizeCont += 5000;
      E_Float* hi2 = new E_Float [sizeCont];
      for (E_Int i = 0; i < N; i++) hi2[i] = hi[i];
      delete [] hi;
      hi = hi2;
    }
  }

  // smooth hi
  /*
  for (E_Int i = 0; i < N-1; i++)
  {
    hi[i] = 0.5*(hi[i+1]+hi[i]);
  }
  hi[N-1] = 0.5*(hi[N-2]+hi[N-1]);
  */

  E_Float r = Ltot-L;
  if (r > 0.5*hi[N-1] && N >= 3) { N = N-1; Ltot = Ltot - hi[N]; }
  
  //printf("hi = ");
  //for (E_Int i = 0; i < N; i++) printf("%g ", hi[i]);
  //printf("\n");

  // scale hi to match full L
  E_Float s = L-Ltot;
  s = s/N; 
  for (E_Int i = 0; i < N; i++) hi[i] = hi[i]+s;
  
  // DBX
  //Ltot = 0.;
  //for (E_Int i = 0; i < N; i++) Ltot += hi[i];
  //printf("Ltot=%g %g\n", Ltot, L);
  // ENDDBX

  // Get ue
  ue = new E_Float [N+1];
  ues = u0;
  ue[0] = u0;
  for (E_Int i = 0; i < N; i++)
  {
    GCPnts_AbscissaPoint Pt(1.e-10, geomAdap, hi[i], ues);
    ues = Pt.Parameter();
    ue[i+1] = ues;
  }
  ue[N] = u1; // forced

  //printf("N= %d\n", N);
  //printf("ue5= ");
  //for (E_Int i = 0; i < N; i++) printf("%g ", ue[i]);
  //printf("\n");

  delete [] hi;
  nbPoints = N+1;
  return 1;
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
  
  E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geomAdap, u0, u1);
  E_Float abscissa = 0.;
  for (E_Int i = 0; i < nbPoints; i++)
  {
    abscissa = uext[i]*L;
    GCPnts_AbscissaPoint Pt(1.e-12, geomAdap, abscissa, u0);
    ue[i] = Pt.Parameter();
    
    // Maybe faster but less accurate
    //h = uext[i]*L - abscissa;
    //GCPnts_AbscissaPoint Pt(1.e-12, geomAdap, h, us);
    //us = Pt.Parameter();
    //ue[i] = us;
    //abscissa = uext[i]*L;
  }
  //for (E_Int i = 0; i < nbPoints; i++) printf("ue: %d : %f %f\n", i, ue[i], uext[i]);
  
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

  E_Float* pex = fe.begin(1);
  E_Float* pey = fe.begin(2);
  E_Float* pez = fe.begin(3);
  E_Float* peu = fe.begin(4);

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
        //px[nbPoints-i-1] = Pt.X(); py[nbPoints-i-1] = Pt.Y(); pz[nbPoints-i-1] = Pt.Z();
        px[nbPoints-i-1] = pex[i]; py[nbPoints-i-1] = pey[i]; pz[nbPoints-i-1] = pez[i];
        pu[nbPoints-i-1] = Puv.X(); pv[nbPoints-i-1] = Puv.Y();
      }
      else
      {
        //px[i] = Pt.X(); py[i] = Pt.Y(); pz[i] = Pt.Z();
        px[i] = pex[i]; py[i] = pey[i]; pz[i] = pez[i];
        pu[i] = Puv.X(); pv[i] = Puv.Y();
      }
    }
    return 1;
  }
  
  // non degenerated    
  {
    gp_Pnt Pt; gp_Pnt2d Puv; E_Float u;
    for (E_Int i = 1; i <= nbPoints; i++)
    {
      u = peu[i-1];
      C0.D0(u, Pt);
      pCurve->D0(u, Puv);
      if (reverse)
      { 
        //px[nbPoints-i] = Pt.X(); py[nbPoints-i] = Pt.Y(); pz[nbPoints-i] = Pt.Z();
        px[nbPoints-i] = pex[i-1]; py[nbPoints-i] = pey[i-1]; pz[nbPoints-i] = pez[i-1];
        pu[nbPoints-i] = Puv.X(); pv[nbPoints-i] = Puv.Y(); 
      }
      else
      {
        //px[i-1] = Pt.X(); py[i-1] = Pt.Y(); pz[i-1] = Pt.Z();
        px[i-1] = pex[i-1]; py[i-1] = pey[i-1]; pz[i-1] = pez[i-1];
        pu[i-1] = Puv.X(); pv[i-1] = Puv.Y(); 
      }
    }
  }
  return 0;
}

// ============================================================================
/* Mesh one edge of CAD, return STRUCT
      hmax
   or hausd
   or hmax + hausd
   or N
   or external param ue
   not used if hmax=-1, hausd=-1, N=-1, ue=None
 */
// ============================================================================
PyObject* K_OCC::meshOneEdge(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int i;
  E_Float hmin; E_Float hmax; E_Float hausd; 
  E_Int N; PyObject* externalEdge;
  if (!PYPARSETUPLE_(args, O_ I_ RRR_ I_ O_, &hook, &i, &hmin, &hmax, &hausd, &N, &externalEdge)) return NULL;
    
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  // Get the ith edge (start 1)
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  const TopoDS_Edge& E = TopoDS::Edge(edges(i));
  
  // use neighbouring faces to evaluate curvature of edge
  bool useFaces = true;
  TopoDS_Shape* shape = (TopoDS_Shape*)packet[0];

  E_Int nbPoints = 0; // nbre of points of discretized edge
  E_Float* ue; // edge param

  if (N > 0 && externalEdge == Py_None)
  {
    nbPoints = N;
    __getUniform(E, nbPoints, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }
  else if ( ((hausd > 0 && std::abs(hmax-hmin) < 1.e-12 && hmax > 0) ||
            (hausd < 0 && hmax > 0)) &&
            externalEdge == Py_None ) // pure hmax
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
  else if (not useFaces && hmax > 0 && hausd > 0 && externalEdge == Py_None) // mix hmax + hausd
  {
    // courbure uniquement sur edge
    __getParamHminHmaxHausdE4(E, hmin, hmax, hausd, nbPoints, ue);
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }
  else if (useFaces && hmax > 0 && hausd > 0 && externalEdge == Py_None) // mix hmax + hausd
  {
    // courbure sur le max des faces
    __getParamHminHmaxHausdF5(E, hmin, hmax, hausd, nbPoints, ue, shape);    
    PyObject* o = K_ARRAY::buildArray2(4, "x,y,z,u", nbPoints, 1, 1, 1);
    FldArrayF* f; K_ARRAY::getFromArray2(o, f);
    __meshEdge(E, nbPoints, ue, *f, false);
    delete [] ue;
    RELEASESHAREDS(o, f);
    return o;
  }

  else if (externalEdge != Py_None) // external parametrization
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
