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

#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <StdFail_NotDone.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <BRep_Tool.hxx>

// Return properties of face
void getShapeProperties(TopoDS_Face& F,
    E_Float& U0, E_Float& U1, E_Float& V0, E_Float& V1,
    E_Int& isUClosed, E_Int& isVClosed, 
    E_Int& isUPeriodic, E_Int& isVPeriodic,
    E_Float& uPeriod, E_Float& vPeriod)
{
  Handle(Geom_Surface) face = BRep_Tool::Surface(F);
  isUClosed = isVClosed = false;
  ShapeAnalysis::GetFaceUVBounds(F, U0, U1, V0, V1);
  printf("Face bounds U=" SF_F2_ " - V=" SF_F2_ "\n",U0,U1,V0,V1);
  ShapeAnalysis_Surface o(face);
  isUClosed = o.IsUClosed();
  isVClosed = o.IsVClosed();
  BRepAdaptor_Surface h(F);
  isUPeriodic = h.IsUPeriodic();
  if (isUPeriodic) uPeriod = h.UPeriod();
  else uPeriod = U1-U0;
  isVPeriodic = h.IsVPeriodic();
  if (isVPeriodic) vPeriod = h.VPeriod();
  else vPeriod = V1-V0;
  // Comme on dessine uniquement sur la premiere periode:
  if (uPeriod > U1-U0+1.e-2) isUClosed = false;
  if (vPeriod > V1-V0+1.e-2) isVClosed = false;
  printf("isUClosed=" SF_D_ " isVClosed=" SF_D_ " isUPeriodic=" SF_D_ " isVPeriodic=" SF_D_ "\n",
         isUClosed, isVClosed, isUPeriodic, isVPeriodic);
}

// Compute u,v du point P sur la face face
E_Int parameters2(const E_Float* pt, TopoDS_Face& F,
    E_Float& u, E_Float& v, 
    E_Int index, E_Float& up, E_Float& vp, E_Float& upp, E_Float& vpp)
{
  Handle(Geom_Surface) face = BRep_Tool::Surface(F);
  E_Float U0, U1, V0, V1, uPeriod, vPeriod;
  E_Int isUClosed, isVClosed, isUPeriodic, isVPeriodic;
  getShapeProperties(F, U0, U1, V0, V1, isUClosed, isVClosed, 
        isUPeriodic, isVPeriodic, uPeriod, vPeriod);
  
  u=v=-1;

  gp_Pnt Point;
  Point.SetCoord(pt[0], pt[1], pt[2]);
  
  ShapeAnalysis_Surface s(face);
  gp_Pnt2d uv, uvp;
  
  if (up == -K_CONST::E_MAX_FLOAT) // if starting
  {
    uv = s.ValueOfUV(Point, 1.e-6);
    u = uv.X(); v = uv.Y();
    
    if (isUPeriodic == true)
    {
      // periodic shift
      E_Float per = uPeriod;
      E_Float b0 = (U0-u)/per; E_Float b1 = (U1-u)/per;
      E_Int N = floor(b1);
      if (N < b0-1.e-3) N += 1;
      if (u <= U0-1.e-2 || u >= U1+1.e-2) u = u+N*per;
    }
    if (isVPeriodic == true)
    {
      // periodic shift
      E_Float per = vPeriod;
      E_Float b0 = (V0-v)/per; E_Float b1 = (V1-v)/per;
      E_Int N = floor(b1);
      if (N < b0-1.e-3) N += 1;
      if (v <= V0-1.e-2 || v >= V1+1.e-2) v = v+N*per;
    }
  }
  else
  {
    uvp.SetCoord(up,vp);
    uv = s.NextValueOfUV (uvp, Point, 1.e-6, -1.0);
    u = uv.X(); v = uv.Y();
    
    if (upp == -K_CONST::E_MAX_FLOAT) { upp = up; vpp = vp; }
    
    // Gestion de la periodicite
    // Avoid backsteping
    E_Float du,dv,dup,dvp,p,inv;
    
    if (isUPeriodic == true)
    {
      // periodic shift
      E_Float per = uPeriod;
      E_Float b0 = (U0-u)/per; E_Float b1 = (U1-u)/per;
      E_Int N = floor(b1);
      if (N < b0-1.e-3) N += 1;
      if (u <= U0-1.e-2 || u >= U1+1.e-2) u = u+N*per;
    }
    if (isVPeriodic == true)
    {
      // periodic shift
      E_Float per = vPeriod;
      E_Float b0 = (V0-v)/per; E_Float b1 = (V1-v)/per;
      E_Int N = floor(b1);
      if (N < b0-1.e-3) N += 1;
      if (v <= V0-1.e-2 || v >= V1+1.e-2) v = v+N*per;
    }
    
    if (isUClosed == true)
    { 
      E_Float per = uPeriod;
      // on the bound
      if (std::abs(u-U0) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);

        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        if (p < -0.2*inv) u = U1;
        else if (std::abs(du) > std::abs(U1-up)) u = U1;
      }
      else if (std::abs(u-U1) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
        
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        if (p < -0.2*inv) u = U0;
        else if (std::abs(du) > std::abs(U0-up)) u = U0;
      }
      else if (std::abs(u-U0-per) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
        
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        if (p < -0.2*inv) u = U0;
        else if (std::abs(du) > std::abs(U0-up)) u = U0;
      }
    }
    if (isVClosed == true)
    { 
      E_Float per = vPeriod;
      if (std::abs(v-V0) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          

        if (p < -0.2*inv) v = V1;
        else if (std::abs(dv) > std::abs(V1-vp)) v = V1;
      }
      else if (std::abs(v-V1) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          
        
        if (p < -0.2*inv) v = V0;
        else if (std::abs(dv) > std::abs(V0-vp)) v = V0;
      }
      else if (std::abs(v-V0-per) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          
        
        if (p < -0.2*inv) v = V0;
        else if (std::abs(dv) > std::abs(V0-vp)) v = V0;
      }
    }
  }
  
  if (up == -K_CONST::E_MAX_FLOAT) { up = u; vp = v; }
  upp = up; vpp = vp; // save for next time
  up = u; vp = v; // save for next time
  
  //__normalize(u,v);
  
  return 0;
}  
