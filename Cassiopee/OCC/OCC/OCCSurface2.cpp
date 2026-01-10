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

#include "OCCSurface.h"
#include "BRep_Tool.hxx"
#include "GeomLib_Tool.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
# include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/IdTool.h"
#include "TopExp_Explorer.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS.hxx"
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <StdFail_NotDone.hxx>
#include "String/kstring.h"

//#define DEBUG_CAD_READER

// Calcule une discretisation de la surface avec connect QUAD
void K_OCC::OCCSurface::discretize(K_FLD::FloatArray& coord3D, K_FLD::IntArray& connect, E_Int ni, E_Int nj)
{
  gp_Pnt Pt;
  coord3D.resize(3, ni*nj);
  connect.resize(4, (ni-1)*(nj-1));
  // Get UV bounds
  E_Float U0, V0, U1, V1, U, V;
  U0 = _U0; U1 = _U1; V0 = _V0; V1 = _V1;
  //ShapeAnalysis::GetFaceUVBounds(_F, U0, U1, V0, V1); // why core?
  for (E_Int j = 0; j < nj; j++)
  for (E_Int i = 0; i < ni; i++)
  {
    U = U0+i*(U1-U0)/(ni-1);
    V = V0+j*(V1-V0)/(nj-1);
    _surface->D0(U, V, Pt);
    coord3D(0,i+j*ni) = Pt.X(); coord3D(1,i+j*ni) = Pt.Y(); coord3D(2,i+j*ni) = Pt.Z();
  }
  for (E_Int j = 0; j < nj-1; j++)
  for (E_Int i = 0; i < ni-1; i++)
  {
    connect(0,i+j*(ni-1)) = i+j*ni;
    connect(1,i+j*(ni-1)) = i+1+j*ni;
    connect(2,i+j*(ni-1)) = i+1+(j+1)*ni;
    connect(3,i+j*(ni-1)) = i+(j+1)*ni;
  }
}

// Projete coord3D sur la surface
void K_OCC::OCCSurface::project(K_FLD::FloatArray& coord3D) const
{
  for (E_Int i=0; i < coord3D.cols(); ++i)
  {
    E_Float x,y,z;
    x = coord3D(0,i); y = coord3D(1,i); z = coord3D(2,i);
    gp_Pnt Point;
    Point.SetCoord(x, y, z);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    /*
    E_Int nbsol = o.NbPoints();
    printf("Nb solution=%d\n",o.NbPoints());
    for (E_Int i = 1; i <= nbsol; i++)
    {
      gp_Pnt Pj = o.Point(i);
      printf("Pt %g %g %g -> %g %g %g\n",pt[0],pt[1],pt[2],Pj.X(),Pj.Y(),Pj.Z());
      o.Parameters(i,u3,v3);
      printf("proj    %d: %g %g\n",i,u3,v3);
    }
    */
    gp_Pnt Pj = o.NearestPoint();
    //printf("projection %f %f %f -> %f %f %f\n",x,y,z,Pj.X(),Pj.Y(),Pj.Z());
    coord3D(0,i) = Pj.X(); coord3D(1,i) = Pj.Y(); coord3D(2,i) = Pj.Z();
  }
}

// Trouve l'element suivant pour un parcours en elements de la BAR
E_Int K_OCC::OCCSurface::findNextElement(E_Int e, K_FLD::IntArray& found, 
  K_FLD::IntArray& connectB, std::vector< std::vector<E_Int> >& node2Elt) const
{
  std::vector< E_Int > next;
  std::set< E_Int > nexts;
  E_Int ind0 = connectB(0,e);
  E_Int ind1 = connectB(1,e);
  for (size_t v = 0; v < node2Elt[ind0].size(); v++)
  {
    E_Int ne = node2Elt[ind0][v];
    if (found[ne] == 0) 
      //next.push_back(ne);
      nexts.insert(ne);
  }
  for (size_t v = 0; v < node2Elt[ind1].size(); v++)
  {
    E_Int ne = node2Elt[ind1][v];
    if (found[ne] == 0) 
      //next.push_back(ne);
      nexts.insert(ne);
  }
  
  // Cherche les seams si plusieurs next
  // Je ne suis pas sur de ne pas louper certains trucs
  E_Int i = 0;
  for (std::set<E_Int>::iterator it=nexts.begin(); it!=nexts.end(); ++it)
  { next.push_back(*it); i++; }
      
  size_t s = next.size();
  if (s == 0) return -1; // pas de choix
  if (s == 1) return next[0]; // un seul choix
  
  std::vector< E_Int > score(s);
  
  for (size_t i = 0; i < s; i++)
  {
    score[i] = 0.;
    // Si on redescend sur l'element precedent
    if ( connectB(0,next[i]) == connectB(0,e) && 
         connectB(1,next[i]) == connectB(1,e) ) score[i] -= 1;
    if ( connectB(1,next[i]) == connectB(0,e) && 
         connectB(0,next[i]) == connectB(1,e) ) score[i] -= 1;
    
    for (size_t j = 0; i < s; i++)
    {
      if (i != j)
      {
        // Si c'est une seam
        if ( connectB(0,next[i]) == connectB(0,next[j]) && 
             connectB(1,next[i]) == connectB(1,next[j]) ) score[i] += 2;
        if ( connectB(0,next[i]) == connectB(1,next[j]) && 
             connectB(1,next[i]) == connectB(0,next[j]) ) score[i] += 2;
      }
    }
  }
  
  E_Int bestScore = -10;
  for (size_t i = 0; i < s; i++) bestScore = std::max(bestScore, score[i]);
  
  printf("score: ");
  for (size_t i = 0; i < s; i++)
    printf(SF_D_ " (" SF_D2_ ") ", score[i], connectB(0,next[i]), connectB(1,next[i]));
  printf("\n");
    
  for (size_t i = 0; i < s; i++) 
    if (score[i] == bestScore) return next[i];
  return next[0];
}

// Parcours la BAR en elements
// Trouve l'enchainement des elements
// Ajoute ensuite des noeuds si besoin
void K_OCC::OCCSurface::parcoursBAR(K_FLD::FloatArray& pos3D, K_FLD::IntArray& connectB)
{
  E_Int npts = pos3D.cols();
  E_Int nelts = connectB.cols();
  K_FLD::IntArray found(E_Int(1), nelts, E_Int(0)); 
  
  // Node2Elt
  E_Int ind0, ind1;
  std::vector< std::vector<E_Int> > node2Elt(npts);
  for (E_Int e = 0; e < connectB.cols(); e++)
  {
    ind0 = connectB(0,e);
    ind1 = connectB(1,e);
    node2Elt[ind0].push_back(e);
    node2Elt[ind1].push_back(e);
  }
  printf("Ordering elements\n"); fflush(stdout);
  
  // Recherche un enchainement des BARS
  // si possible en une seule boucle
  std::vector<E_Int> eltChain(nelts);
  for (E_Int i = 0; i < nelts; i++) eltChain[i] = -1;
    
  E_Int elt = 0; E_Int next;
  for (E_Int i = 0; i < nelts; i++)
  {
    found[elt] = 1;
    eltChain[i] = elt;
    next = findNextElement(elt, found, connectB, node2Elt);
    if (next == -1) break;
    elt = next;
  }
  
  for (E_Int i = 0; i < nelts; i++)
    printf("chain " SF_D_ ": " SF_D_ " (" SF_D2_ ")\n",i,eltChain[i],connectB(0,eltChain[i]),connectB(1,eltChain[i]));
  
  // Ajoute les noeuds de valence > 2
  /*
  K_FLD::IntArray valence(1,npts,0); 
  E_Int add = npts;
  for (E_Int e = 0; e < nelts; e++)
  {
    E_Int ind0 = connectB(e,0);
    E_Int ind1 = connectB(e,1);
    valence[ind0] += 1;
    valence[ind1] += 1;
    if (valence[ind0] > 2) { connectB(e,0) = add; add++; }// ajoute pts dans pos3D}
    if (valence[ind1] > 2) { connectB(e,1) = add; add++; }// ajoute point dans pos3D}
  }
  */
}

// En entree : les elements sont doubles sur le backbone
// Tentative de duplication brutale
// Les pts de valence =3 ou 4 sont dupliques
// Ils sont remplaces dans deux elements au pif suivant switch
void K_OCC::OCCSurface::dupBAR(K_FLD::FloatArray& pos3D, K_FLD::IntArray& connectB,
  K_FLD::IntArray& switcha, std::map< E_Int, E_Int >& mirror)
{
  E_Int npts = pos3D.cols();
  //E_Int nelts = connectB.cols();
  
  // Node2Elt
  E_Int ind0, ind1;
  std::vector< std::vector<E_Int> > node2Elt(npts);
  for (E_Int e = 0; e < connectB.cols(); e++)
  {
    ind0 = connectB(0,e);
    ind1 = connectB(1,e);
    node2Elt[ind0].push_back(e);
    node2Elt[ind1].push_back(e);
  }
  
  // save connectB
  K_FLD::IntArray connectB2 = connectB;
  
  E_Int added = npts;
  E_Float* Pt; E_Int e;
  
  // Valence
  for (E_Int i = 0; i < npts; i++)
  {
    E_Int valence = node2Elt[i].size();
    //printf("%d: valence=%d\n", i, valence);
    if (valence == 3)
    {
      Pt = pos3D.col(i);
      pos3D.pushBack(Pt, Pt+3);
      E_Int e0 = node2Elt[i][0];
      E_Int e1 = node2Elt[i][1];
      //E_Int e2 = node2Elt[i][2];
      E_Int imodif = 2;
      if (  (connectB2(0,e0) == connectB2(1,e1) && connectB2(1,e0) == connectB2(0,e1))
         || (connectB2(0,e0) == connectB2(0,e1) && connectB2(1,e0) == connectB2(1,e1)) )
      { imodif = 1; }
      e = node2Elt[i][imodif];
      if (connectB(0,e) == i) connectB(0,e) = added;
      else if (connectB(1,e) == i) connectB(1,e) = added;
      mirror[added] = i;
      added += 1;
    }
    else if (valence == 4)
    {
      Pt = pos3D.col(i);
      pos3D.pushBack(Pt, Pt+3);
      E_Int e0 = node2Elt[i][0];
      E_Int e1 = node2Elt[i][1];
      E_Int e2 = node2Elt[i][2];
      E_Int e3 = node2Elt[i][3];
      
      //if (i == 40)
      //{
      //  printf("e0: %d %d\n", connectB2(0,e0), connectB2(1,e0));
      //  printf("e1: %d %d\n", connectB2(0,e1), connectB2(1,e1));
      //  printf("e2: %d %d\n", connectB2(0,e2), connectB2(1,e2));
      //  printf("e3: %d %d\n", connectB2(0,e3), connectB2(1,e3));
      //}
      E_Int imodif1 = 2; E_Int imodif2 = 3;
      if (  (connectB2(0,e0) == connectB2(1,e1) && connectB2(1,e0) == connectB2(0,e1))
         || (connectB2(0,e0) == connectB2(0,e1) && connectB2(1,e0) == connectB2(1,e1)) )
      {
        
        if (switcha[i] == 0) { imodif1 = 1; imodif2 = 2; }
        else { imodif1 = 1; imodif2 = 3; }
      }
      else if (  (connectB2(0,e0) == connectB2(1,e2) && connectB2(1,e0) == connectB2(0,e2))
              || (connectB2(0,e0) == connectB2(0,e2) && connectB2(1,e0) == connectB2(1,e2)) )
      {
        if (switcha[i] == 0) { imodif1 = 2; imodif2 = 1; }
        else { imodif1 = 2; imodif2 = 3; }
      }
      else if (  (connectB2(0,e0) == connectB2(1,e3) && connectB2(1,e0) == connectB2(0,e3))
              || (connectB2(0,e0) == connectB2(0,e3) && connectB2(1,e0) == connectB2(1,e3)) )
      {
        if (switcha[i] == 0) { imodif1 = 3; imodif2 = 1; }
        else { imodif1 = 3; imodif2 = 2; }
      }
      else if (  (connectB2(0,e1) == connectB2(1,e2) && connectB2(1,e1) == connectB2(0,e2))
              || (connectB2(0,e1) == connectB2(0,e2) && connectB2(1,e1) == connectB2(1,e2)) )
      {
        if (switcha[i] == 0) { imodif1 = 2; imodif2 = 3; }
        else { imodif1 = 2; imodif2 = 0; }
      }
      else if (  (connectB2(0,e1) == connectB2(1,e3) && connectB2(1,e1) == connectB2(0,e3))
              || (connectB2(0,e1) == connectB2(0,e3) && connectB2(1,e1) == connectB2(1,e3)) )
      {
        if (switcha[i] == 0) { imodif1 = 3; imodif2 = 2; }
        else { imodif1 = 3; imodif2 = 0; }
      }
      else if (  (connectB2(0,e2) == connectB2(1,e3) && connectB2(1,e2) == connectB2(0,e3))
              || (connectB2(0,e2) == connectB2(0,e3) && connectB2(1,e2) == connectB2(1,e3)) )
      {
        if (switcha[i] == 0) { imodif1 = 3; imodif2 = 1; }
        else { imodif1 = 3; imodif2 = 0; }
      }
      //if (i == 294)
      //  printf("%d: imodif = %d et %d\n", i, imodif1, imodif2);
      e = node2Elt[i][imodif1];
      if (connectB(0,e) == i) connectB(0,e) = added;
      else if (connectB(1,e) == i) connectB(1,e) = added;
      e = node2Elt[i][imodif2];
      if (connectB(0,e) == i) connectB(0,e) = added;
      else if (connectB(1,e) == i) connectB(1,e) = added;
      mirror[added] = i;
      added += 1; 
    }
  }
}

// findNextPoint. En priorite les 1, sinon les 2
E_Int K_OCC::OCCSurface::findNextPoint(K_FLD::IntArray& found, 
  std::vector< std::vector<E_Int> >& node2Elt) const
{
  E_Int i0 = -1; E_Int i1 = -1; E_Int i2 = -1;
  for (E_Int i = 0; i < found.cols(); i++)
  {
    if (found[i] == 0)
    {
      if (i0 == -1) i0 = i;
      if (node2Elt[i].size() == 1) { i1 = i; break; }
      if (node2Elt[i].size() == 2 && i2 == -1) i2 = i;
    }
  }
  if (i1 >= 0) return i1;
  if (i2 >= 0) return i2;
  return i0;
}

// Find non ambiguous start
E_Int K_OCC::OCCSurface::findNonAmbStart(E_Int npts, K_FLD::FloatArray& coord3D) const
{
  bool amb;
  E_Float u,v,up,vp,upp,vpp;
  for (E_Int i = 0; i < npts; i++)
  {
    up = -K_CONST::E_MAX_FLOAT; upp = -K_CONST::E_MAX_FLOAT;
    parameters2(coord3D.col(i), u, v, i, up, vp, upp, vpp);
    amb = false;
    if (_isUClosed)
    {
      if (std::fabs(u-0.) < 1.e-2) amb = true;
      if (std::fabs(u-1.) < 1.e-2) amb = true;
    }
    if (_isVClosed)
    {
      if (std::fabs(v-0.) < 1.e-2) amb = true;
      if (std::fabs(v-1.) < 1.e-2) amb = true;
    }
    if (amb == false) { return i; }
  }
  return 0;
}


// order the BAR
void K_OCC::OCCSurface::orderBAR(E_Int npts, 
  K_FLD::FloatArray& coord3D, 
  K_FLD::IntArray& connectB, 
  K_FLD::IntArray& index, K_FLD::IntArray& start) const
{
  // get node->elts connect
  E_Int ind0, ind1, ind2, ind3;
  std::vector< std::vector<E_Int> > node2Elt(npts);
  for (E_Int e = 0; e < connectB.cols(); e++)
  {
    ind0 = connectB(0,e);
    ind1 = connectB(1,e);
    node2Elt[ind0].push_back(e);
    node2Elt[ind1].push_back(e);
  }
  //printf("after connect\n"); fflush(stdout);
  /*
  for (E_Int i = 0; i <npts; i++)
  {
    if (node2Elt[i].size() == 0)
      printf("node2Elt %d: rien!\n", i);
    else if (node2Elt[i].size() == 1)
      printf("node2Elt %d : %d\n", i, node2Elt[i][0]);
    else if  (node2Elt[i].size() == 2)
      printf("node2Elt %d : %d %d\n", i, node2Elt[i][0], node2Elt[i][1]);
    else printf("node2Elt %d : plus de 3\n", i);
  }
  */
  // check open BAR, branch
#ifdef DEBUG_CAD_READER
  for (E_Int i = 0; i < npts; i++)
  {
    if (node2Elt[i].size() == 0) // isolated
    {
      // pas tres grave
      //printf("Warning: %d: isolated point\n", i);
    }
    else if (node2Elt[i].size() == 1) // open BAR
    {
      printf("Warning: opened BAR.\n");
      printf(SF_D_ ": elt=" SF_D_ "\n", i, node2Elt[i][0]);
    }
    else if (node2Elt[i].size() > 2) // Branch
    {
      printf("Warning: branch.\n");
      printf(SF_D_ ": ", i);
      for (size_t j = 0; j < node2Elt[i].size(); j++) printf("elt=" SF_D_ " ", node2Elt[i][j]);
      printf("\n");
    }
  }
  //printf("after check\n"); fflush(stdout);
#endif
  
  // order
  E_Int s, e0, e1;
  index.resize(1, npts);
  start.resize(E_Int(1), npts, E_Int(0));
  K_FLD::IntArray found(E_Int(1), npts, E_Int(0));
   
  // NEW ALGO
  //E_Int indCur = findNextPoint(found, node2Elt);
  E_Int indCur = findNonAmbStart(npts, coord3D);
#ifdef DEBUG_CAD_READER
  printf("starting non ambiguous index=" SF_D_ "\n", indCur);
#endif
  start[0] = indCur;
  
  for (E_Int i = 0; i < npts; i++)
  {
    // Traitement indCur
    index[i] = indCur;
    found[indCur] = 1;
    
    // find elements
    s = node2Elt[indCur].size();
    if (s == 0) 
    { 
      indCur = findNextPoint(found, node2Elt); 
      if (i < npts-1) start[i+1] = 1; 
    }
    else if (s == 1)
    {
      e0 = node2Elt[indCur][0];
      ind0 = connectB(0, e0);
      ind1 = connectB(1, e0);
      if (ind0 != indCur && found[ind0] == 0) indCur = ind0;
      else if (ind1 != indCur && found[ind1] == 0) indCur = ind1;
      else 
      { 
        indCur = findNextPoint(found, node2Elt); 
        if (i < npts-1) start[i+1] = 1; 
      }
    }
    else
    {
      e0 = node2Elt[indCur][0];
      e1 = node2Elt[indCur][1];
      ind0 = connectB(0, e0);
      ind1 = connectB(1, e0);
      ind2 = connectB(0, e1);
      ind3 = connectB(1, e1);
      if (ind0 != indCur && found[ind0] == 0) indCur = ind0;
      else if (ind1 != indCur && found[ind1] == 0) indCur = ind1;
      else if (ind2 != indCur && found[ind2] == 0) indCur = ind2;
      else if (ind3 != indCur && found[ind3] == 0) indCur = ind3;
      else 
      { 
        indCur = findNextPoint(found, node2Elt); 
        if (i < npts-1) start[i+1] = 1; 
      }
    }
  }
  
  // Check order new->old
#ifdef DEBUG_CAD_READER
  for (E_Int i = 0; i < npts; i++)
    printf("order new=" SF_D_ " -> old=" SF_D_ " (found=" SF_D_ ", start=" SF_D_ ") \n", i, index[i], found[index[i]], start[i]);
#endif
  
  // Check start
  E_Int nstart = 0;
  E_Int length = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    if ( (start[i] == 1 && i != 0) || i == npts-1)
    {
      if (node2Elt[index[i]].size() > 0) {  nstart += 1;  }
      length = 0;
    }
    else length += 1;
  }
#ifdef DEBUG_CAD_READER
  if (nstart > 1) printf("Warning: nstart = " SF_D_ "\n", nstart);
#endif
}


// parametre le contour (coord3D, connectB) avec les parametres de la surface
E_Int K_OCC::OCCSurface::parameters2
(K_FLD::FloatArray& coord3D, K_FLD::IntArray& connectB, K_FLD::FloatArray& UVs) const
{
  UVs.clear();
  
  E_Int sz(coord3D.cols());
  UVs.resize(2, sz, K_CONST::E_MAX_FLOAT);

  K_FLD::IntArray index;
  K_FLD::IntArray start;
  
  E_Int npts = coord3D.cols();
  //printf("before orderBAR\n"); fflush(stdout);
  orderBAR(npts, coord3D, connectB, index, start);
  //printf("after orderBAR\n"); fflush(stdout);
#ifdef DEBUG_CAD_READER
  printf("bounds %f %f - %f %f \n",_U0,_U1,_V0,_V1);
  printf("isClosedU = " SF_D_ " isClosedV = " SF_D_ "\n",_isUClosed,_isVClosed);
  printf("isUPeriodic = " SF_D_ " isVPeriodic = " SF_D_ "\n", _isUPeriodic, _isVPeriodic);
  printf("UPeriod=%f, VPeriod=%f\n", _uPeriod,_vPeriod);
#endif
  
  E_Int n;
  E_Float Up=-K_CONST::E_MAX_FLOAT; E_Float Vp=-K_CONST::E_MAX_FLOAT;
  E_Float Upp=-K_CONST::E_MAX_FLOAT; E_Float Vpp=-K_CONST::E_MAX_FLOAT;
  
  for (E_Int i = 0; i < npts; i++)
  {
    if (start[i] == 1) { Up=-K_CONST::E_MAX_FLOAT; Vp=-K_CONST::E_MAX_FLOAT; Upp=-K_CONST::E_MAX_FLOAT; Vpp=-K_CONST::E_MAX_FLOAT; }
    n = index[i];
    parameters2(coord3D.col(n), UVs(0,n), UVs(1,n), n, Up, Vp, Upp, Vpp);
    if (_isRevol == true && (std::fabs(Up-Upp) > 0.7*(_U1-_U0) || std::fabs(Vp-Vpp) > 0.7*(_V1-_V0)))
    {
#ifdef DEBUG_CAD_READER
      printf("Warning: %f %f | %f %f Jump detected in " SF_D_ ".\n",Up,Upp,Vp,Vpp,n);
#endif
      return index[i-1]+1;
    }
    
    //printf("" SF_D_ "/" SF_D_ ": %f %f \n",n,npts,UVs(0,n),UVs(1,n));
  }
  return 0;
  
  // Detect jumps
  /*
  for (E_Int i=0; i < connectB.cols(); ++i)
  {
    E_Int Ni = connectB(0,i);
    E_Int Nj = connectB(1,i);
    E_Float Ui = UVs(0, Ni);
    E_Float Vi = UVs(1, Ni);
    E_Float Uj = UVs(0, Nj);
    E_Float Vj = UVs(1, Nj);
    //printf("la " SF_D_ ": %g %g -> %g %g\n", i,Ui,Vi,Uj,Vj);
    if (::fabs(Ui-Uj) > 0.7 || ::fabs(Vi-Vj) > 0.7) printf("switch elt=" SF_D_ ", n1=" SF_D_ ", n2=" SF_D_ ": %g %g -> %g %g\n", i,Ni,Nj,Ui,Vi,Uj,Vj);
  }
  */
  
  // ??
  /*
  for (E_Int i=0; i < connectB.cols(); ++i)
  {
    E_Int Ni = connectB(0,i);
    E_Int Nj = connectB(1,i);
    if (_isUClosed && ::fabs(UVs(0,Ni) - UVs(0,Nj)) > K_CONST::E_PI) err = 1;
    if (_isVClosed && ::fabs(UVs(1,Ni) - UVs(1,Nj)) > K_CONST::E_PI) err = 1;
  }
  */
  
  // Correction de robustesse
  /*
  for (E_Int i = 0; i < npts; i++)
  {
    if (start[i] == 1) { Up=-1; Vp=-1; Upp=-1; Vpp=-1; }
    n = index[i];
    u = UVs(0,n); v = UVs(1,n);
    //printf("" SF_D_ ": %f %f \n",n,UVs(0,n),UVs(1,n));
    Upp = Up; Vpp = Vp;
    Up = UVs(0,n); Vp = UVs(1,n);
  }
  */
  // Ultimate security
  /*
  for (E_Int k=0; k < UVs.cols(); ++k)
    if (UVs(0,k) == K_CONST::E_MAX_FLOAT)
      UVs(0,k) = UVs(1,k) = 0.;
  */
  
}

// parameters2
E_Int
K_OCC::OCCSurface::parameters2(const E_Float* pt, E_Float& u, E_Float& v, 
  E_Int index, E_Float& up, E_Float& vp, E_Float& upp, E_Float& vpp) const 
{
  u=v=-1;
    
  gp_Pnt Point;
  Point.SetCoord(pt[0], pt[1], pt[2]);
  
  /* Another cunning development */
  ShapeAnalysis_Surface s(_surface);
  gp_Pnt2d uv, uvp;
  
  if (up == -K_CONST::E_MAX_FLOAT) // starting
  {
    uv = s.ValueOfUV(Point, 1.e-6);
    u = uv.X(); v = uv.Y();
    
    if (_isUPeriodic == true)
    {
      // periodic shift
      E_Float per = _uPeriod;
      E_Float b0 = (_U0-u)/per; E_Float b1 = (_U1-u)/per;
      E_Int N = floor(b1);
      //if (N < b0-1.e-2) printf("Warning: danger %d %f!\n",N,b0);
      if (N < b0-1.e-3) N += 1;
      if (u <= _U0-1.e-2 || u >= _U1+1.e-2) u = u+N*per;
      //E_Int N = floor(std::abs(u-_U0)/per);
      //if (u > _U1+1.e-2) u = u - N*per;
      //if (u < _U0-1.e-2) u = u + (N+1)*per;
    }
    if (_isVPeriodic == true)
    {
      // periodic shift
      E_Float per = _vPeriod;
      E_Float b0 = (_V0-v)/per; E_Float b1 = (_V1-v)/per;
      E_Int N = floor(b1);
      //if (N < b0-1.e-2) printf("Warning: danger %d %f!\n",N,b0);
      if (N < b0-1.e-3) N += 1;
      if (v <= _V0-1.e-2 || v >= _V1+1.e-2) v = v+N*per;
      //E_Int N = floor(std::abs(v-_V0)/per);
      //if (v > _V1+1.e-2) v = v - N*per;
      //if (v < _V0-1.e-2) v = v + (N+1)*per;
    }
    //E_Float u1,v1;
    //GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    //o.LowerDistanceParameters(u1,v1);
    //if (::fabs(u-u1) > 1.e-6 || ::fabs(v-v1) > 1.e-6) printf("erreur: %f %f versus %f %f\n",u,v,u1,v1);
#ifdef DEBUG_CAD_READER
    printf(SF_D_ ": startuv: %f %f (orig=%f %f)\n",index,u,v,uv.X(),uv.Y());
#endif
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
    
    if (_isUPeriodic == true)
    {
      // periodic shift
      E_Float per = _uPeriod;
      E_Float b0 = (_U0-u)/per; E_Float b1 = (_U1-u)/per;
      E_Int N = floor(b1);
      //if (N < b0-1.e-3) printf("Warning: danger N=%d < %f!\n",N,b0);
      if (N < b0-1.e-3) N += 1;
      if (u <= _U0-1.e-2 || u >= _U1+1.e-2) u = u+N*per;
      //E_Int N = floor(std::abs(u-_U0)/per);
      //if (u > _U1+1.e-2) u = u - N*per;
      //if (u < _U0-1.e-2) u = u + (N+1)*per;
      //if (_isUPeriodic) u = u + ShapeAnalysis::AdjustToPeriod(u, _U0, _U1);
    }
    if (_isVPeriodic == true)
    {
      // periodic shift
      E_Float per = _vPeriod;
      E_Float b0 = (_V0-v)/per; E_Float b1 = (_V1-v)/per;
      E_Int N = floor(b1);
      //if (N < b0-1.e-3) printf("Warning: danger N=%d %f!\n",N,b0);
      if (N < b0-1.e-3) N += 1;
      if (v <= _V0-1.e-2 || v >= _V1+1.e-2) v = v+N*per;
      //E_Int N = floor(std::abs(v-_V0)/per);
      //if (v > _V1+1.e-2) v = v - N*per;
      //if (v < _V0-1.e-2) v = v + (N+1)*per;
      //if (_isVPeriodic) v = v + ShapeAnalysis::AdjustToPeriod(v, _V0, _V1);
    }
    
    if (_isUClosed == true)
    { 
      E_Float per = _uPeriod;
      // on the bound
      if (std::fabs(u-_U0) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);

        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        //printf("switchU0 %f %f\n", p, 0.2*inv);
        //if (index == 97) printf("p=%f,du=%f,dv =%f,dup=%f,dvp=%f\n",p,du,dv,dup,dvp);
        if (p < -0.2*inv) u = _U1;
        else if (std::fabs(du) > std::fabs(_U1-up)) u = _U1;
      }
      else if (std::fabs(u-_U1) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
        
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        //printf("switchU1 %f %f\n", p, _U0-up);
        if (p < -0.2*inv) u = _U0;
        else if (std::fabs(du) > std::fabs(_U0-up)) u = _U0;
      }
      else if (std::fabs(u-_U0-per) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
        
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;
        //printf("switchU1 %f %f\n", p, _U0-up);
        if (p < -0.2*inv) u = _U0;
        else if (std::fabs(du) > std::fabs(_U0-up)) u = _U0;
      }
    }
    if (_isVClosed == true)
    { 
      E_Float per = _vPeriod;
      if (std::fabs(v-_V0) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          
          
        if (p < -0.2*inv) v = _V1;
        else if (std::fabs(dv) > std::fabs(_V1-vp)) v = _V1;
      }
      else if (std::fabs(v-_V1) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          
        
        if (p < -0.2*inv) v = _V0;
        else if (std::fabs(dv) > std::fabs(_V0-vp)) v = _V0;
      }
      else if (std::fabs(v-_V0-per) < 1.e-2)
      {
        du = (u-up); dv = (v-vp);
        dup = (up-upp); dvp = (vp-vpp);
      
        inv = std::sqrt(du*du+dv*dv)*std::sqrt(dup*dup+dvp*dvp)+1.e-10;
        p = du*dup+dv*dvp;          
        
        if (p < -0.2*inv) v = _V0;
        else if (std::fabs(dv) > std::fabs(_V0-vp)) v = _V0;
      }
    }
#ifdef DEBUG_CAD_READER
    printf(SF_D_ ": nextuv: %f %f (orig=%f %f)\n",index, u,v,uv.X(),uv.Y());
#endif
  }
  
  /*
  uv = s.ValueOfUV(Point, 1.e-6);
  u = uv.X(); v = uv.Y();
  printf("%f %f - %f %f %f %f\n", u,v,_U0,_U1,_V0,_V1);
  */
  if (up == -K_CONST::E_MAX_FLOAT) { up = u; vp = v; }
  upp = up; vpp = vp; // save for next time
  up = u; vp = v; // save for next time
  
  __normalize(u,v);
  
  return 0;
  
  // cherche a projeter
  bool fail;
  E_Float du, dv;
  
  // State start
  fail = false;
  E_Float un,vn;
  E_Float step1,stepA,stepB;
  GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
  try { o.LowerDistanceParameters(u,v); }
  catch( StdFail_NotDone& e ) { fail = true; }
  
  //printf("%f %f %f -> %f %f\n",pt[0],pt[1],pt[2],u,v);
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0]+1.e-6, pt[1], pt[2]);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0]-1.e-6, pt[1], pt[2]);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0], pt[1]+1.e-6, pt[2]);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0], pt[1]-1.e-6, pt[2]);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0], pt[1], pt[2]+1.e-6);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  if (fail == true)
  {
    // try perturbations
    Point.SetCoord(pt[0], pt[1], pt[2]-1.e-6);
    GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
    try { o.LowerDistanceParameters(u,v); }
    catch( StdFail_NotDone& e ) { fail = true; }
  }
  
  
  if (fail == true)
  {
    u = 0; v = 0; return 1;
  }
  
  // normalize in [0,1]
  __normalize(u,v);
  
  if (up < 0 || vp < 0) return 0;
  
  if (_isUClosed == false && _isVClosed == false) return 0;
  
  du = std::abs(u-up);
  dv = std::abs(v-vp);
  step1 = du+dv;
  
  // Recherche pour le periodique
  un = u; vn = v;
  
  if (_isUClosed)
  {
    if (::std::fabs(u - 0.)<2.e-2) un = 1.;
    if (::std::fabs(u - 1.)<2.e-2) un = 0.;
    du = std::fabs(un-up);
    dv = std::fabs(v-vp);
    stepA = du+dv;
    if (stepA < 0.9*step1) { u = un; }
  }
  if (_isVClosed)
  {
    if (::std::fabs(v - 0.)<2.e-2) vn = 1.;
    if (::std::fabs(v - 1.)<2.e-2) vn = 0.;
    du = std::fabs(u-up);
    dv = std::fabs(vn-vp);
    stepB = du+dv;
    if (stepB < 0.9*step1) { v = vn; }
  }
  
  return 0;
  
  
  //E_Int ok = GeomLib_Tool::Parameters(_surface, Point, 1.e+100/*dummytol*/, u, v);
  
  
  // Optimize search depending on previous Up,Vp
  
  //if (!ok || u==-1 || v==-1) return 1;
  
  //E_Float u1,v1;
  //GeomLib_Tool::Parameters(_surface, Point, 1.e+100/*dummytol*/, u1, v1);
  
  //GeomAPI_ProjectPointOnSurf o(Point, _surface, Extrema_ExtAlgo_Tree);
  //E_Float u2,v2;
  //o.LowerDistanceParameters(u2,v2);
  
  //printf("Point index=%d\n",index);
  //printf("geomlib:   %g %g\n",u1,v1);
  //printf("lowerproj: %g %g\n",u2,v2);
  //printf("=============\n");
  
  //if (index == 4)
  //{
  //  printf("Point index=%d\n",index);
  //  printf("geomlib:   %g %g\n",u1,v1);
  //  printf("lowerproj: %g %g\n",u2,v2);
  //  
  //  E_Float u3,v3;
  //  E_Int nbsol = o.NbPoints();
  //  printf("Nb solution=%d\n",o.NbPoints());
  //  for (E_Int i = 1; i <= nbsol; i++)
  //  {
  //    gp_Pnt Pj = o.Point(i);
  //    printf("Pt %g %g %g -> %g %g %g\n",pt[0],pt[1],pt[2],Pj.X(),Pj.Y(),Pj.Z());
  //    o.Parameters(i,u3,v3);
  //    printf("proj    %d: %g %g\n",i,u3,v3);
  //  }
  //  //printf("=============\n");
  //}
  //u = u2; v = v2;
    
}
