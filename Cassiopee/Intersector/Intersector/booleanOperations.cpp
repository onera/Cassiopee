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

// Some boolean operations like intersection, union, minus...
//#define FLAG_STEP
//#define DEBUG_W_PYTHON_LAYER
//#define DEBUG_MESHER
//#define DEBUG_BOOLEAN

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/include/TRI_BooleanOperator.h"
# include "Nuga/include/BAR_BooleanOperator.h"
# include "Nuga/include/NGON_BooleanOperator.h"
#include <memory>

using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
#include "Nuga/include/chrono.h"
#endif
#if defined(DEBUG_TRIANGULATOR)
    bool DELAUNAY::Triangulator::dbg_enabled = false;
#endif
#ifdef DEBUG_BOOLEAN
    std::string medith::wdir = "./";
#endif

//=============================================================================
/* Factorisation de la récupération et de la validation des arguments. */
//=============================================================================

enum eOperation {INTERSECTION, UNION, MINUS12, INTERSECTION_BORDER, MODIFIED_SOLID, DIFFSURF};

std::string getName(eOperation oper)
{
  switch (oper)
  {
    case INTERSECTION:
      return "booleanIntersection";
    case UNION:
      return "booleanUnion";
    case MINUS12:
      return "booleanMinus";
    case INTERSECTION_BORDER:
      return "booleanIntersectionBorder";
    case MODIFIED_SOLID:
      return "booleanModifiedSolid";
    case DIFFSURF:
      return "Diffsurf";
    default:
      return "Unknown";
  }
}

//==============================================================================
// Return true is eltType is ok for given oper
bool is_valid(char* eltType, eOperation oper)
{
  switch (oper)
  {
    case INTERSECTION:
    case UNION:
    case MINUS12:
    case INTERSECTION_BORDER:
    case MODIFIED_SOLID:
      if ((strcmp(eltType, "TRI") == 0) || (strcmp(eltType, "BAR") == 0) || (strcmp(eltType, "NGON") == 0))
        return true;
      else return false;
    case DIFFSURF:
    {
      if (strcmp(eltType, "NGON") == 0)
        return true;
      else return false;
    }
    default: return false;
  }
}

//==============================================================================
// Return false si error in args, true otherwise
bool getArgs(PyObject* args, eOperation oper,
             K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
             K_FLD::FloatArray& pos2, K_FLD::IntArray& connect2,
             E_Float& tolerance,
             E_Int& preserve_right, E_Int& solid_right, E_Int& agg_mode, bool& improve_conformal_cloud_qual, bool& outward_surf, E_Int& itermax, char*& eltType, char*& varString)
{
  PyObject *arrS[2];
  E_Float tol = 0.;
  E_Int preserv_r, solid_r, imp_qual(0), out_sur(0);
  std::ostringstream o;
  std::string opername = getName(oper);

  if (!PYPARSETUPLE_(args, OO_ R_ IIII_ II_,
                    &arrS[0], &arrS[1], &tol, &preserv_r, &solid_r, &agg_mode, &imp_qual, &out_sur, &itermax))
  {
    o << opername << ": wrong arguments.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    return false;
  }

  improve_conformal_cloud_qual = bool(imp_qual);
  outward_surf = bool(out_sur);

  E_Int ni, nj, nk, posx[2], posy[2], posz[2], err(0);
  K_FLD::FloatArray *fS[2];
  K_FLD::IntArray* cn[2];

  // Check the Arguments.
  for (E_Int n = 0; n < 2; n++)
  {
    // The input surface must be a unstructured triangles surface.
    E_Int res = K_ARRAY::getFromArray(arrS[n], varString, fS[n], ni, nj, nk, cn[n], eltType);
    
    if ((res != 2) || (!is_valid(eltType, oper)))
    {
      o << opername << ": invalid array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }

    // Check coordinates.
    posx[n] = K_ARRAY::isCoordinateXPresent(varString);
    posy[n] = K_ARRAY::isCoordinateYPresent(varString);
    posz[n] = K_ARRAY::isCoordinateZPresent(varString);

    if ((posx[n] == -1) || (posy[n] == -1) || (posz[n] == -1))
    {
      o << opername << ": can't find coordinates in array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }
  }
  
  if (!err)
  {
    pos1 = *fS[0];
    if (posx[0] !=0 || posy[0] != 1 || posz[0] != 2)
    {
      for (E_Int i=0; i < pos1.cols(); ++i)
      {
        pos1(0,i) = (*fS[0])(posx[0],i);
        pos1(1,i) = (*fS[0])(posy[0],i);
        pos1(2,i) = (*fS[0])(posz[0],i);    
      }
    }

    pos2 = *fS[1];
    if (posx[1] !=0 || posy[1] != 1 || posz[1] != 2)
    {
      for (E_Int i=0; i < pos2.cols(); ++i)
      {
        pos2(0,i) = (*fS[1])(posx[1],i);
        pos2(1,i) = (*fS[1])(posy[1],i);
        pos2(2,i) = (*fS[1])(posz[1],i);    
      }
    }
    
    connect1 = *cn[0];
    connect2 = *cn[1];
    tolerance = tol;
    preserve_right = preserv_r;
    solid_right=solid_r;
  }
  
  // Final cleaning.
  for (E_Int n = 0; n < 2; n++)
  {
    delete fS[n]; delete cn[n];
  }
  return (err == 0);
}

//==============================================================================
// Return true if args are ok
bool getBorderArgs(PyObject* args,
             K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
             K_FLD::FloatArray& pos2, K_FLD::IntArray& connect2,
             E_Float& tolerance, int& itermax,
             char*& eltType, char*& varString)
{
  PyObject *arrS[2];
  tolerance = 0.;
  std::ostringstream o;

  if (!PYPARSETUPLE_(args, OO_ R_ I_, 
                    &arrS[0], &arrS[1], &tolerance, &itermax))
  {
    o << "booleanIntersectionBorder" << ": wrong arguments.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    return false;
  }

  E_Int ni, nj, nk, posx[2], posy[2], posz[2], err(0);
  K_FLD::FloatArray *fS[2];
  K_FLD::IntArray* cn[2];

  // Check the Arguments.
  for (E_Int n = 0; n < 2; n++)
  {
    // The input surface must be a unstructured triangles surface.
    E_Int res = K_ARRAY::getFromArray(arrS[n], varString, fS[n], ni, nj, nk, cn[n], eltType);
    
    if ((res != 2) || (!is_valid(eltType, INTERSECTION_BORDER)))
    {
      o << "booleanIntersectionBorder" << ": invalid array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }

    // Check coordinates.
    posx[n] = K_ARRAY::isCoordinateXPresent(varString);
    posy[n] = K_ARRAY::isCoordinateYPresent(varString);
    posz[n] = K_ARRAY::isCoordinateZPresent(varString);

    if ((posx[n] == -1) || (posy[n] == -1) || (posz[n] == -1))
    {
      o << "booleanIntersectionBorder" << ": can't find coordinates in array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }
  }
  
  if (!err)
  {
    pos1 = *fS[0];
    if (posx[0] !=0 || posy[0] != 1 || posz[0] != 2)
    {
      for (E_Int i=0; i < pos1.cols(); ++i)
      {
        pos1(0,i) = (*fS[0])(posx[0],i);
        pos1(1,i) = (*fS[0])(posy[0],i);
        pos1(2,i) = (*fS[0])(posz[0],i);    
      }
    }

    pos2 = *fS[1];
    if (posx[1] !=0 || posy[1] != 1 || posz[1] != 2)
    {
      for (E_Int i=0; i < pos2.cols(); ++i)
      {
        pos2(0,i) = (*fS[1])(posx[1],i);
        pos2(1,i) = (*fS[1])(posy[1],i);
        pos2(2,i) = (*fS[1])(posz[1],i);    
      }
    }
    
    connect1 = *cn[0];
    connect2 = *cn[1];
  }
  
  // Final cleaning.
  for (E_Int n = 0; n < 2; n++)
  {
    delete fS[n]; delete cn[n];
  }
  return (err == 0);
}

//==============================================================
// Return true if args are ok
bool getUnionArgs(PyObject* args,
             K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
             K_FLD::FloatArray& pos2, K_FLD::IntArray& connect2,
             E_Float& tolerance, E_Int& preserve_right, E_Int& solid_right, E_Int& agg_mode, bool& improve_conformal_cloud_qual,
             std::vector<E_Int>& pgsList, E_Int& simplify_pgs, E_Int& hard_mode, E_Int &itermax,
             char*& eltType, char*& varString)
{
  PyObject *arrS[2], *pgs;
  E_Float tol = 0.;
  E_Int preserv_r, solid_r, imp_qual(0);
  std::ostringstream o;
  std::string opername = getName(UNION);

  pgsList.clear();

  if (!PYPARSETUPLE_(args, OO_ R_ IIII_ O_ III_, 
                    &arrS[0], &arrS[1], &tol, &preserv_r, &solid_r, &agg_mode, &imp_qual, &pgs, &simplify_pgs, &hard_mode, &itermax))
  {
    o << opername << ": wrong arguments.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    return false;
  }

  improve_conformal_cloud_qual = bool(imp_qual);

  E_Int ni, nj, nk, posx[2], posy[2], posz[2], err(0);
  K_FLD::FloatArray *fS[2];
  K_FLD::IntArray* cn[2];

  // Check the Arguments.
  for (E_Int n = 0; n < 2; n++)
  {
    // The input surface must be a unstructured triangles surface.
    E_Int res = K_ARRAY::getFromArray(arrS[n], varString, fS[n], ni, nj, nk, cn[n], eltType);
    
    if ((res != 2) || (!is_valid(eltType, UNION)))
    {
      o << opername << ": invalid array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }

    // Check coordinates.
    posx[n] = K_ARRAY::isCoordinateXPresent(varString);
    posy[n] = K_ARRAY::isCoordinateYPresent(varString);
    posz[n] = K_ARRAY::isCoordinateZPresent(varString);

    if ((posx[n] == -1) || (posy[n] == -1) || (posz[n] == -1))
    {
      o << opername << ": can't find coordinates in array.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      err = 1;
    }
  }
  
  if (!err)
  {
    pos1 = *fS[0];
    if (posx[0] !=0 || posy[0] != 1 || posz[0] != 2)
    {
      // std::cout << "posx 0 : " << posx[0] << std::endl;
      // std::cout << "posy 0 : " << posy[0] << std::endl;
      // std::cout << "posz 0 : " << posz[0] << std::endl;

      for (E_Int i=0; i < pos1.cols(); ++i)
      {
        pos1(0,i) = (*fS[0])(posx[0],i);
        pos1(1,i) = (*fS[0])(posy[0],i);
        pos1(2,i) = (*fS[0])(posz[0],i);    
      }
    }

    pos2 = *fS[1];
    if (posx[1] !=0 || posy[1] != 1 || posz[1] != 2)
    {
      // std::cout << "posx 1 : " << posx[1] << std::endl;
      // std::cout << "posy 1 : " << posy[1] << std::endl;
      // std::cout << "posz 1 : " << posz[1] << std::endl;

      for (E_Int i=0; i < pos2.cols(); ++i)
      {
        pos2(0,i) = (*fS[1])(posx[1],i);
        pos2(1,i) = (*fS[1])(posy[1],i);
        pos2(2,i) = (*fS[1])(posz[1],i);    
      }
    }
    
    connect1 = *cn[0];
    connect2 = *cn[1];
    tolerance = tol;
    preserve_right = preserv_r;
    solid_right=solid_r;
  }

  // Passing the specified pgs to the boolean to extrude polygons over contacts
  {
    E_Int res=0;
    K_FLD::FldArrayI* inds=NULL;
    if (pgs != Py_None)
      res = K_NUMPY::getFromNumpyArray(pgs, inds);

    std::unique_ptr<K_FLD::FldArrayI> pL(inds); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY is : " << res << std::endl;
    if ((res == 1) && (inds != NULL)  && (inds->getSize() != 0))
    {
      E_Int nb_ghost_pgs = inds->getSize();
      //E_Int minid(INT_MAX), maxid(-1);
      pgsList.resize(nb_ghost_pgs);
      for (E_Int i = 0; i < nb_ghost_pgs; ++i) 
      {
        pgsList[i]=(*inds)[i]-1;
        //std::cout << pgsList[i] << std::endl;
        //minid = std::min(minid, pgsList[i]);
        //maxid = std::max(maxid, pgsList[i]);
      }
    }
  }
  
  // Final cleaning.
  for (E_Int n = 0; n < 2; n++)
  {
    delete fS[n]; delete cn[n];
  }
  return (err == 0);
}

//=============================================================================
/* Factorisation des opérations booléennes. */
//=============================================================================
typedef E_Int (BooleanOperator::* operation)(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors);

template <typename OperType>
operation getOperation(eOperation oper)
{
  switch (oper)
  {
    case INTERSECTION:        
      return reinterpret_cast<operation>(&OperType::getIntersection);
    case UNION :              
      return reinterpret_cast<operation>(&OperType::getUnion);
    case MINUS12:             
      return reinterpret_cast<operation>(&OperType::get_1_minus_2);
    default: return NULL;
  }
}

//==============================================================================
// Main function that perform all kind of operations
PyObject* call_operation(PyObject* args, eOperation oper)
{
  K_FLD::FloatArray pos1, pos2, pos;
  K_FLD::IntArray connect1, connect2, connect;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors;
  E_Int preserve_right, solid_right, agg_mode, itermax;
  bool improve_conformal_cloud_qual(false), outward_surf(true);

  bool ok = getArgs(args, oper, pos1, connect1, pos2, connect2, tolerance,
		    preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, outward_surf, itermax, eltType, varString);
  if (!ok) return NULL;

  PyObject* tpl = NULL;
  E_Int err(0), et=-1;

  char eltType2[20];
  strcpy(eltType2, eltType);

  if ((strcmp(eltType, "TRI") == 0) || (strcmp(eltType, "BAR") == 0))
  {
    BooleanOperator* BO = NULL;
    operation op;
    
    if (strcmp(eltType, "TRI") == 0)
    {
      BO = new TRI_BooleanOperator (pos1, connect1, pos2, connect2, tolerance, itermax);
      op = getOperation<TRI_BooleanOperator>(oper);
    }
    else if (strcmp(eltType, "BAR") == 0)
    {
      BO = new BAR_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
      op = getOperation<BAR_BooleanOperator>(oper);
    }
    
    if (oper != INTERSECTION_BORDER)
      err = (BO->*op)(pos, connect, colors);
    else
    {
      err = BO->getIntersectionBorder(pos, connect);
      if (strcmp(eltType, "TRI") == 0) strcpy(eltType2, "BAR");
      else if (strcmp(eltType, "BAR") == 0) strcpy(eltType2, "NODE");
    }
    
    delete BO;
    et=-1;
  }
  else if (strcmp(eltType, "NGON") == 0)
  {
    typedef NGON_BooleanOperator<K_FLD::FloatArray, K_FLD::IntArray> bo_t;
    
    bo_t::eAggregation AGG=bo_t::CONVEX;
    if (agg_mode == 0) // NONE
    {
      //AGG = bo_t::NONE; fixme : not working because of new_skin PG triangulation : diconnexion between modified and unlodified parts
    }
    else if (agg_mode == 2) AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual) BO.setConformizerParams(true);

    if (oper == DIFFSURF && outward_surf == false) BO._outward = false;

    switch (oper)
    {
    case INTERSECTION:        
      /*err=*/BO.Intersection(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case UNION :              
      /*err=*/BO.Union(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case MINUS12:             
      /*err=*/BO.Diff(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case MODIFIED_SOLID:             
      /*err=*/BO.Modified_Solid(pos, connect, (bo_t::eMergePolicy)preserve_right); break;
    case DIFFSURF:             
      /*err=*/BO.Diffsurf(pos, connect); break;
    default: return NULL;
    }
    et = 8;
  }
  
  if (err != 1)
    tpl = K_ARRAY::buildArray(pos, varString, connect, et, eltType2, false);
  else
  {
    std::ostringstream o;
    o << getName(oper) << ": failed to proceed.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    return NULL;
  }
  return tpl;
}

PyObject* call_xborder(PyObject* args)
{
  K_FLD::FloatArray pos1, pos2, pos;
  K_FLD::IntArray connect1, connect2, connect;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors;
  int itermax=10;
  
  bool ok = getBorderArgs(args, pos1, connect1, pos2, connect2, tolerance, itermax, eltType, varString);
  if (!ok) return NULL;
  PyObject* tpl = NULL;
  E_Int err(0), et=-1;

  char eltType2[20];
  strcpy(eltType2, eltType);

  if ((strcmp(eltType, "TRI") == 0) || (strcmp(eltType, "BAR") == 0))
  {
    BooleanOperator* BO = NULL;
           
    if (strcmp(eltType, "TRI") == 0)
      BO = new TRI_BooleanOperator (pos1, connect1, pos2, connect2, tolerance, itermax);
    else if (strcmp(eltType, "BAR") == 0)
      BO = new BAR_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
    
    err = BO->getIntersectionBorder(pos, connect);
    if (strcmp(eltType, "TRI") == 0) strcpy(eltType2, "BAR");
    else if (strcmp(eltType, "BAR") == 0) strcpy(eltType2, "NODE");
    
    delete BO;
    et=-1;
  }
  
  if (err != 1)
    tpl = K_ARRAY::buildArray(pos, varString, connect, et, eltType2, false);
  else
  {
    std::ostringstream o;
    o << "getIntersectionBorder" << ": failed to proceed.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
  }
  return tpl;
}

//==============================================================
// Fonction generale d'union
//==============================================================
PyObject* call_union(PyObject* args)
{
  K_FLD::FloatArray pos1, pos2;
  K_FLD::IntArray connect1, connect2;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors, ghost_pgs;
  E_Int preserve_right, solid_right, agg_mode;
  bool improve_conformal_cloud_qual(false);
  E_Int simplify_pgs{1};
  E_Int hard_mode{0}, itermax{10};
  
  bool ok = getUnionArgs(args, pos1, connect1, pos2, connect2, tolerance, 
        preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, ghost_pgs, simplify_pgs, hard_mode, itermax, eltType, varString);
  if (!ok) return NULL;
  
  PyObject* tpl  = NULL;
  PyObject* tplph0 = NULL, *tplph1 = NULL;
  PyObject* tplpg0 = NULL, *tplpg1 = NULL;
  E_Int err(0), et=-1;

#ifdef DEBUG_W_PYTHON_LAYER
  PyObject *l(PyList_New(0));
#endif

  char eltType2[20];
  strcpy(eltType2, eltType);

  std::vector<K_FLD::FloatArray> crds;
  std::vector<K_FLD::IntArray> cnts;
  std::vector<std::string> znames;

  if ((strcmp(eltType, "TRI") == 0) || (strcmp(eltType, "BAR") == 0))
  {
    BooleanOperator* BO = NULL;
    operation op;
    
    if (strcmp(eltType, "TRI") == 0)
    {
      BO = new TRI_BooleanOperator (pos1, connect1, pos2, connect2, tolerance, itermax);
      op = getOperation<TRI_BooleanOperator>(UNION);
    }
    else if (strcmp(eltType, "BAR") == 0)
    {
      BO = new BAR_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
      op = getOperation<BAR_BooleanOperator>(UNION);
    }
    crds.resize(1); cnts.resize(1);
    err = (BO->*op)(crds[0], cnts[0], colors);
    delete BO;
    et = -1;
  }
  else if (strcmp(eltType, "NGON") == 0)
  {
    typedef NGON_BooleanOperator<K_FLD::FloatArray, K_FLD::IntArray> bo_t;
    
    bo_t::eAggregation AGG=bo_t::CONVEX;
    if (agg_mode == 0) // NONE
    {
      //AGG = bo_t::NONE; fixme : not working because of new_skin PG triangulation : diconnexion between modified and unlodified parts
    }
    else if (agg_mode == 2)
      AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual) BO.setConformizerParams(true);

    if (!ghost_pgs.empty())
      BO.passPGs(1, ghost_pgs);

    if (improve_conformal_cloud_qual) BO.setTriangulatorParams(/*do not shuffle : */false, /*improve quality:*/true);

    BO.simplify_pgs = (simplify_pgs == 1);

    BO.hard_mode = hard_mode;

    K_FLD::IntArray cnt;
    K_FLD::FloatArray crd;
    err = BO.Union(crd, cnt, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right);
    et = 8;

#ifndef DEBUG_W_PYTHON_LAYER
    crds.push_back(crd);
    cnts.push_back(cnt);

    ngon_type & ngo = *BO._ngoper;
       
    size_t nb_phs = ngo.PHs.size();
    size_t nb_pgs = ngo.PGs.size();

    std::vector<E_Int> phoids0(nb_phs);
    std::vector<E_Int> phoids1(nb_phs);
    std::ostringstream o;
    
    for (size_t i=0; i < nb_phs; ++i)
    {
	
      E_Int ancPH1 = ngo.PHs._ancEs(0,i);
      E_Int ancPH2 = ngo.PHs._ancEs(1,i);
      
      assert (i >=0 && i < (size_t)ngo.PHs._ancEs.cols());
      
      phoids0[i] = ancPH1;
      phoids1[i] = ancPH2;
      
      assert (i >=0 && i < phoids0.size());
      assert (i >=0 && i < phoids1.size());
      
      if (ancPH1 == E_IDX_NONE && ancPH2 == E_IDX_NONE)
      {
    	o << "Error: Union: Cell " << i << " without ancestor." << std::endl; 
    	PyErr_SetString(PyExc_TypeError, o.str().c_str());
    	err = 1;
      }
      
      if (ancPH1 != E_IDX_NONE && ancPH2 != E_IDX_NONE)
      {
    	o << "Error: Union: 2 ancestors for one cell, not handled yet.";
    	PyErr_SetString(PyExc_TypeError, o.str().c_str());
    	err = 1;
      }
      
    }

    // Reverse indirection
    std::vector<E_Int> phnids0;
    std::vector<E_Int> phnids1;
    
    K_CONNECT::IdTool::reverse_indirection(phoids0, phnids0);

    for (size_t i=0; i < phnids0.size(); ++i)
    {
      if (phnids0[i] == E_IDX_NONE){ phnids0[i] = -1 ;}
    }
    
    K_CONNECT::IdTool::reverse_indirection(phoids1, phnids1);

    for (size_t i=0; i < phnids1.size(); ++i)
    {
      if (phnids1[i] == E_IDX_NONE){ phnids1[i] = -1 ;}
    }
    
    tplph0 = K_NUMPY::buildNumpyArray(&phnids0[0], phnids0.size(), 1, 0);
    tplph1 = K_NUMPY::buildNumpyArray(&phnids1[0], phnids1.size(), 1, 0);

    
    std::vector<E_Int> pgoids0(nb_pgs);
    std::vector<E_Int> pgoids1(nb_pgs);

    K_FLD::IntArray F2E ; 
    ngo.build_noF2E(F2E);
        
    for (size_t i=0; i < nb_pgs; ++i)
    {      
      assert (i >=0 && i < (size_t)ngo.PGs._ancEs.cols());
      
      E_Int ancPG1 = ngo.PGs._ancEs(0,i);
      E_Int ancPG2 = ngo.PGs._ancEs(1,i);
      
      assert (i >=0 && i < pgoids0.size());
      assert (i >=0 && i < pgoids1.size());
      
      pgoids0[i] = ancPG1;
      pgoids1[i] = ancPG2;

      if (ancPG1 == E_IDX_NONE && ancPG2 == E_IDX_NONE)
      {
    	  std::ostringstream o;
    	  o << "Error: Union: Face " << i << " without ancestor." << std::endl;
    	  PyErr_SetString(PyExc_TypeError, o.str().c_str());
    	  return NULL;	  
      }
      
      if (ancPG1 != E_IDX_NONE && ancPG2 != E_IDX_NONE)
      {
    	  std::ostringstream o;
    	  o << "Error: Union: Face " << i << " has 2 ancestors. ancPG1 = " << ancPG1 << " and ancPG2 =" <<  ancPG2 << std::endl;
    	  PyErr_SetString(PyExc_TypeError, o.str().c_str());
    	  return NULL;
      }

      // Detect interior faces 
      E_Int elt1 = F2E(0,i);
      E_Int elt2 = F2E(1,i);
      
      if (elt1 != E_IDX_NONE && elt2 != E_IDX_NONE)
      {
	      pgoids0[i] = E_IDX_NONE; 
	      pgoids1[i] = E_IDX_NONE; 
      }
 
    }

    tplpg0 = K_NUMPY::buildNumpyArray(&pgoids0[0], pgoids0.size(), 1, 0);
    tplpg1 = K_NUMPY::buildNumpyArray(&pgoids1[0], pgoids1.size(), 1, 0);
    

#else //DEBUG_W_PYTHON_LAYER
    {
     
      // boolean answer
      if (err != 0)
      {
        crds.push_back(crd);
        cnts.push_back(cnt);
        znames.push_back("union");
      }

      // Conformized
      if (BO.conform_T3s.ng.PGs.size())
      {
        crds.push_back(BO.conform_T3s.crd);
        BO.conform_T3s.ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("conformized");
      }

      // PHT3s init (SOFT)
      if (BO.PHT3s_begin[0].ng.PGs.size())
      {
        crds.push_back(BO.PHT3s_begin[0].crd);
        BO.PHT3s_begin[0].ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_begin");
      }

      // parasite PHT3s (SOFT)
      if (BO.PHT3s_parasite[0].ng.PGs.size())
      {
        crds.push_back(BO.PHT3s_parasite[0].crd);
        BO.PHT3s_parasite[0].ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_parasites");
      }

      // irrelevant PHT3s (SOFT only)
      if (BO.PHT3s_irrelevant.ng.PGs.size())
      {
        crds.push_back(BO.PHT3s_irrelevant.crd);
        BO.PHT3s_irrelevant.ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_irrelevant");
      }

      // unclosed PHT3s (SOFT)
      if (BO.PHT3s_unclosed[0].ng.PGs.size())
      {
        crds.push_back(BO.PHT3s_unclosed[0].crd);
        BO.PHT3s_unclosed[0].ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_unclosed");
      }

      // histo failure PHT3s (SOFT)
      if (BO.PHT3s_history[0].ng.PGs.size())
      {
        //std::cout << "histo size : " << BO.PHT3s_history[0].ng.PHs.size() << std::endl;
        crds.push_back(BO.PHT3s_history[0].crd);
        BO.PHT3s_history[0].ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_histo_failure");
      }

      // PHT3s final (SOFT)
      if (BO.PHT3s_end[0].ng.PGs.size())
      {
        crds.push_back(BO.PHT3s_end[0].crd);
        BO.PHT3s_end[0].ng.export_to_array(cnt);
        cnts.push_back(cnt);
        znames.push_back("soft_PHT3s_end");
      }

    }
#endif

  } // strcmp(eltType, "NGON") == 0 
  
#ifndef DEBUG_W_PYTHON_LAYER
  {
    if (err != 1)
    {
      tpl = K_ARRAY::buildArray(crds[0], varString, cnts[0], et, eltType2, false);
    }
    else
    {
      std::ostringstream o;
      o << "Union: failed to proceed.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
      return NULL;
    }

    PyObject* l = PyList_New(0);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    if (strcmp(eltType, "NGON") == 0)
    {  
      PyList_Append(l, tplph0);
      Py_DECREF(tplph0);
    
      PyList_Append(l, tplph1);
      Py_DECREF(tplph1);
    
      PyList_Append(l, tplpg0);
      Py_DECREF(tplpg0);
    
      PyList_Append(l, tplpg1);
      Py_DECREF(tplpg1);
    }
    
    return l;
  }
#else
  {
    {
      size_t nb_zones = cnts.size();
      //std::cout << "nb zones: " << nb_zones << std::endl;
      for (size_t i=0; i < nb_zones; ++i)
      {
        if (crds[i].cols() == 0) continue;

      	tpl = K_ARRAY::buildArray(crds[i], varString, cnts[i], et, "NGON", false);

      	PyObject* o = Py_BuildValue("(Os)", tpl, znames[i].c_str());
        PyList_Append(l, o);
        Py_DECREF(tpl);
      }
    }
  }

  return l;

#endif
}

//=============================================================================
/* Intersection de 2 surfaces fermees. */
//=============================================================================
PyObject* K_INTERSECTOR::booleanIntersection(PyObject* self, PyObject* args)
{
  return call_operation(args, INTERSECTION);
}

//=============================================================================
/* Union de 2 surfaces fermees (enveloppe externe). */
//=============================================================================
PyObject* K_INTERSECTOR::booleanUnion(PyObject* self, PyObject* args)
{
  return call_union(args);
}

PyObject* K_INTERSECTOR::booleanUnionMZ(PyObject* self, PyObject* args)
{
  E_Int agg_mode{2}, improve_qual{0}, simplify_pgs{1}, hard_mode{0};
  E_Float xtol(0.), closetol(0.);
  PyObject* arr1s, *arr2s;

  if (!PYPARSETUPLE_(args, OO_ RR_ IIII_, 
                    &arr1s, &arr2s, &xtol, &closetol, &agg_mode, &improve_qual, &simplify_pgs, &hard_mode))
  {
    PyErr_SetString(PyExc_TypeError, "booleanUnion2: wrong args");
    return NULL;
  }

  E_Int nb_zones1 = PyList_Size(arr1s);
  E_Int nb_zones2 = PyList_Size(arr2s);

  std::vector<K_FLD::FloatArray*> crd1s(nb_zones1, nullptr), crd2s(nb_zones2, nullptr);
  std::vector<K_FLD::IntArray*>   cnt1s(nb_zones1, nullptr), cnt2s(nb_zones2, nullptr);
  char* varString{ nullptr }, *eltType{ nullptr };
  PyObject *l(PyList_New(0));

  // get the zones
  for (E_Int i=0; i < nb_zones1; ++i)
  {
    PyObject* py_zone = PyList_GetItem(arr1s, i);
    
    E_Int err = check_is_NGON(py_zone, crd1s[i], cnt1s[i], varString, eltType);
    if (err)
    {
      for (E_Int i=0; i < nb_zones1; ++i)
      {
        delete crd1s[i];
        delete cnt1s[i];
      }
      PyErr_SetString(PyExc_TypeError, "booleanUnionMZ: not NGON elts.");
      return NULL;
    }
  }
  for (E_Int i=0; i < nb_zones2; ++i)
  {
    PyObject* py_zone = PyList_GetItem(arr2s, i);
    
    E_Int err = check_is_NGON(py_zone, crd2s[i], cnt2s[i], varString, eltType);
    if (err)
    {
      for (E_Int i=0; i < nb_zones2; ++i)
      {
        delete crd2s[i];
        delete cnt2s[i];
      }
      PyErr_SetString(PyExc_TypeError, "booleanUnionMZ : not NGON elts.");
      return NULL;
    }
  }

  // join and close as 2 operands but keeping track of zones
  std::vector<E_Int> zonePHids1, zonePHids2;
  std::vector<E_Int> zoneshiftPH1(nb_zones1), zoneshiftPH2(nb_zones2); 
  std::vector<E_Int> zoneshiftPG1(nb_zones1), zoneshiftPG2(nb_zones2);
  
  ngon_type ng1;
  K_FLD::FloatArray crd1;
  for (E_Int i=0; i < nb_zones1; ++i)
  {
    ngon_type ng(*cnt1s[i]);
    ng.PGs.shift(crd1.cols());
    ng1.append(ng);
    zoneshiftPH1[i] = ng1.PHs.size();
    zoneshiftPG1[i] = ng1.PGs.size();
    zonePHids1.resize(ng1.PHs.size(),i);
    crd1.pushBack(*crd1s[i]);
  }
  
  ngon_type ng2;
  K_FLD::FloatArray crd2;
  for (E_Int i=0; i < nb_zones2; ++i)
  {
    ngon_type ng(*cnt2s[i]);
    ng.PGs.shift(crd2.cols());
    ng2.append(ng);
    zoneshiftPH2[i] = ng2.PHs.size();
    zoneshiftPG2[i] = ng2.PGs.size();
    zonePHids2.resize(ng2.PHs.size(),i);
    crd2.pushBack(*crd2s[i]);
  }

  // Clean connectivity and keeping track of histories 
  std::vector<E_Int> clean_pgnids1;
  std::vector<E_Int> clean_pgnids2;
  std::vector<E_Int> clean_phnids1;
  std::vector<E_Int> clean_phnids2;
    
  ngon_type::clean_connectivity(ng1, crd1, -1/*ngon_dim*/, closetol/*tolerance*/, true/*remove_dup_phs*/, true/*do_omp*/, &clean_pgnids1, &clean_phnids1);
  ngon_type::clean_connectivity(ng2, crd2, -1/*ngon_dim*/, closetol/*tolerance*/, true/*remove_dup_phs*/, true/*do_omp*/, &clean_pgnids2, &clean_phnids2);

  // Reverse indirection (preserving information when 2 ancestors exist) 
  map<E_Int,std::vector<E_Int>> vect_pgoids1;
  map<E_Int,std::vector<E_Int>> vect_pgoids2;

  for (size_t i=0; i < clean_pgnids1.size(); ++i)
  {
    E_Int nids = clean_pgnids1[i]; 
    vect_pgoids1[nids].push_back(i); 
  }
  
  for (size_t i=0; i < clean_pgnids2.size(); ++i)
  {
    E_Int nids         = clean_pgnids2[i]; 
    vect_pgoids2[nids].push_back(i); 
  }

  // Union
  K_FLD::IntArray cnt1, cnt2;
  ng1.export_to_array(cnt1);
  ng2.export_to_array(cnt2);

  using bo_t = NGON_BooleanOperator<K_FLD::FloatArray, K_FLD::IntArray>;
  bo_t::eAggregation AGG=bo_t::CONVEX;
  if (agg_mode == 0) // NONE
  {
    //AGG = bo_t::NONE; fixme : not working because of new_skin PG triangulation : diconnexion between modified and unlodified parts
  }
  else if (agg_mode == 2) AGG = bo_t::FULL;
  
  bo_t BO(crd1, cnt1, crd2, cnt2, xtol, AGG);

  if (improve_qual)
  {
    BO.setConformizerParams(true);
    BO.setTriangulatorParams(/*do not shuffle : */false, /*improve quality:*/true);
  }

  BO.simplify_pgs = (simplify_pgs == 1);

  BO.hard_mode = hard_mode;
  
  K_FLD::IntArray cnt;
  K_FLD::FloatArray crd;
  E_Int err=BO.Union(crd, cnt, bo_t::eInterPolicy::SOLID_RIGHT, bo_t::eMergePolicy::PRESERVE_RIGHT);

  if (err == 1) return NULL;

  if (err == -7) // empty X
  {
#ifdef E_DOUBLEINT
   return Py_BuildValue("l", long(err));
#else
  return Py_BuildValue("i", err);
#endif 
  }

  ngon_type & ngo = *BO._ngoper;

  std::vector<K_FLD::FloatArray> crd1so(nb_zones1), crd2so(nb_zones2);
  std::vector<ngon_type> ng1so(nb_zones1), ng2so(nb_zones2);

  // PH history
  std::vector<std::vector<E_Int>> phoids1(nb_zones1), phoids2(nb_zones2); // phoids[zid][i] = k <=> the i-th cell in zid-th zone either had k as id, or was a piece of k (in same zone) 

  // Structure for keeping new matches information 
  std::map<E_Int, std::map<E_Int, std::set<E_Int>>> z1_jz_to_ptl1, z2_jz_to_ptl2;
   
  //
  K_FLD::IntArray F2E; 
  ngo.build_noF2E(F2E);
  //

  // Retrieve zones and PH oids
  E_Int nb_phs = ngo.PHs.size();
  
  for (E_Int i=0; i < nb_phs; ++i)
  {
    E_Int ancPH1 = ngo.PHs._ancEs(0,i);
    E_Int ancPH2 = ngo.PHs._ancEs(1,i);

    // exclusive OR : coming from op1 OR op2
    assert ((ancPH1 != E_IDX_NONE) || (ancPH2 != E_IDX_NONE));
    assert ((ancPH1 == E_IDX_NONE) || (ancPH2 == E_IDX_NONE));

    const E_Int* faces = ngo.PHs.get_facets_ptr(i);
    E_Int nb_faces = ngo.PHs.stride(i);

    if (ancPH2 == E_IDX_NONE) // comes from op1
    {	
      E_Int zid = zonePHids1[ancPH1];

      if (zid > 0) ancPH1 = ancPH1 - zoneshiftPH1[zid-1]; // Shift to get index local to zone 

      ng1so[zid].PHs.add(nb_faces, faces);
      phoids1[zid].push_back(ancPH1);

      // on cherche les nouveaux raccords : entre op1 et op2 (uniquement)
      for (E_Int f = 0; f < nb_faces; ++f)
      {
        E_Int PGi = faces[f]-1;
        E_Int PHleft = F2E(0, PGi);
        E_Int PHright = F2E(1, PGi);
        assert (PHleft == i || PHright == i);
        
        E_Int PHother    = (PHleft == i) ? PHright : PHleft;
        E_Int ancPHother = E_IDX_NONE ;
        if (PHother != E_IDX_NONE)  ancPHother = ngo.PHs._ancEs(1,PHother);
	  
        if (ancPHother == IDX_NONE) continue;// provient de op1 => pgoids en sortie permettra de mettre à jour le ptlist
        E_Int zidother = zonePHids2[ancPHother] + nb_zones1 ; // zid dans op2 + shift (pour num. globale)

        z1_jz_to_ptl1[zid][zidother].insert(PGi);
      }
      
    }
    else  // comes from op2
    {
      E_Int zid = zonePHids2[ancPH2];
      
      if (zid > 0) ancPH2 = ancPH2 - zoneshiftPH2[zid-1]; // Shift to get index local to zone
      
      ng2so[zid].PHs.add(nb_faces, faces);
      phoids2[zid].push_back(ancPH2);
      
      // on cherche les nouveaux raccords : entre op1 et op2 (uniquement)
      for (E_Int f = 0; f < nb_faces; ++f)
      {
        E_Int PGi = faces[f]-1;
        E_Int PHleft = F2E(0, PGi);
        E_Int PHright = F2E(1, PGi);
        assert (PHleft == i || PHright == i);
        
        E_Int PHother = (PHleft == i) ? PHright : PHleft;
        E_Int ancPHother = E_IDX_NONE ;
        if (PHother != E_IDX_NONE)  ancPHother = ngo.PHs._ancEs(0,PHother);
        
        if (ancPHother == IDX_NONE) continue;// provient de op2 => pgoids en sortie permettra de mettre à jour le ptlist
        E_Int zidother = zonePHids1[ancPHother]; // dans op1

        z2_jz_to_ptl2[zid+nb_zones1][zidother].insert(PGi);
      }
    }
  }
  
  // Construction d'un historique pgoids pour l'union global (toutes les zones)
  Vector_t<E_Int> glo_pgoids;
  
  // On change de container parce qu'il ne faut plus de set ici 
  std::map<E_Int, std::map<E_Int, std::vector<E_Int>>> z1loc_jz_to_ptl1;  

  //E_Int shift = 0;
  std::vector<std::vector<E_Int>> idx_to_remove1(nb_zones1);
  
  for (size_t i=0; i < ng1so.size(); ++i)
  {
    //if (i > 0) shift = zoneshiftPG1[i-1];
    crd1so[i] = crd;
    ng1so[i].PGs = ngo.PGs;
    Vector_t<E_Int> pgnids, phnids;
    ng1so[i].remove_unreferenced_pgs(pgnids, phnids);

    // indice global => local
    for (auto& i : z1_jz_to_ptl1)
    {
      E_Int zid = i.first;
      for (auto& j : i.second)
      {
        E_Int jzid = j.first;
        auto & ptl = j.second;
        for (auto & pg : ptl)
        {
          if (pgnids[pg] !=  E_IDX_NONE)
          {
            z1loc_jz_to_ptl1[zid][jzid].push_back(pgnids[pg]);
          }
        }
      }
    }
    
    ngon_type::compact_to_used_nodes(ng1so[i].PGs, crd1so[i]);

    K_FLD::IntArray F2E_loc ; 
    ng1so[i].build_noF2E(F2E_loc);

    for (E_Int k = 0; k < ng1so[i].PGs.size(); ++k)
    {
      E_Int ancPG1 = ng1so[i].PGs._ancEs(0,k);
      glo_pgoids.push_back(ancPG1);
    }
  }
  
  // Remove faces already belonging to matches
  for (auto& ii : z1loc_jz_to_ptl1)
  {
    E_Int zid = ii.first;
    for (auto& j : ii.second)
    {
      //E_Int jzid = j.first;
      auto & ptl = j.second;
      for (auto & pg : ptl)
      {
        if (pg != E_IDX_NONE)
          idx_to_remove1[zid].push_back(pg); 
      }
    }
  }
  
  E_Int prev = 0 ;
  
  for (size_t i=0; i < ng1so.size(); ++i)
  {
    Vector_t<E_Int> pgnids, phnids;

    K_FLD::IntArray cnto;
    ng1so[i].export_to_array(cnto);
  
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(crd1so[i], varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // computing PG histo
    std::vector<E_Int> pgoids(ng1so[i].PGs.size(), IDX_NONE);
    E_Int shift = 0;
      
    for (size_t k = prev; k < prev+pgoids.size(); ++k)
    {
      if (i>0) shift = zoneshiftPG1[i-1];
      
      E_Int ids_u = glo_pgoids[k];

      if (ids_u != E_IDX_NONE)
      {
        pgoids[k-prev] = vect_pgoids1[ids_u][0]-shift ;

        if ( vect_pgoids1[ids_u].size() > 1)
          vect_pgoids1[ids_u].erase(vect_pgoids1[ids_u].begin());
      }
      else pgoids[k-prev] = -1; 
    }

    // Remove pg already flagged as match
    for (size_t kk = 0 ; kk < idx_to_remove1[i].size(); ++kk )
    {
      E_Int id_rm   = idx_to_remove1[i][kk];
      pgoids[id_rm] = -1 ;
    }

    prev = prev + pgoids.size();

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&pgoids[0], pgoids.size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PH history
    tpl = K_NUMPY::buildNumpyArray(&phoids1[i][0], phoids1[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

  }

  // Construction d'un historique pgoids pour l'union global (toutes les zones)
  Vector_t<E_Int> glo_pgoids2;
  std::vector<std::vector<E_Int>> idx_to_remove2(nb_zones2);
  
  // On change de container parce qu'il ne faut plus de set ici 
  std::map<E_Int, std::map<E_Int, std::vector<E_Int>>> z2loc_jz_to_ptl2;
  
  for (size_t i=0; i < ng2so.size(); ++i)
  {
    crd2so[i] = crd;
    ng2so[i].PGs = ngo.PGs;
    Vector_t<E_Int> pgnids, phnids;
    ng2so[i].remove_unreferenced_pgs(pgnids, phnids);

    // indice global => local
    for (auto& i : z2_jz_to_ptl2)  
    {
      E_Int zid = i.first;
      for (auto& j : i.second)
      {
        E_Int jzid = j.first;
        auto & ptl = j.second;
        for (auto & pg : ptl)
        {
          if (pgnids[pg] !=  E_IDX_NONE)
            z2loc_jz_to_ptl2[zid][jzid].push_back(pgnids[pg]);
        }
      }
    }
   
    ngon_type::compact_to_used_nodes(ng2so[i].PGs, crd2so[i]);
    
    for (E_Int k = 0; k < ng2so[i].PGs.size(); ++k)
    {
      E_Int ancPG2 = ng2so[i].PGs._ancEs(1,k);
 
      glo_pgoids2.push_back(ancPG2);
    }
    
  }
  
  // Remove faces belonging to matches
  for (auto& ii : z2loc_jz_to_ptl2)
  {
    E_Int zid = ii.first;
    for (auto& j : ii.second)
    {
      //E_Int jzid = j.first;
      auto & ptl = j.second;
      for (auto & pg : ptl)
      {
        if (pg != E_IDX_NONE) idx_to_remove2[zid-nb_zones1].push_back(pg);
      }
    }
  }
   
  prev = 0;
  
  for (size_t i=0; i < ng2so.size(); ++i)
  {
    Vector_t<E_Int> pgnids, phnids;

    K_FLD::IntArray cnto;
    ng2so[i].export_to_array(cnto);
  
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(crd2so[i], varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // computing PG histo
    std::vector<E_Int> pgoids(ng2so[i].PGs.size(), IDX_NONE);
    E_Int shift = 0;
    
    for (size_t k = prev; k < prev+pgoids.size(); ++k)
    {
      if (i>0) shift = zoneshiftPG2[i-1] ;

      E_Int ids_u = glo_pgoids2[k];

      if (ids_u != E_IDX_NONE)
      {
        pgoids[k-prev] = vect_pgoids2[ids_u][0]-shift ;

        if ( vect_pgoids2[ids_u].size() > 1)
          vect_pgoids2[ids_u].erase(vect_pgoids2[ids_u].begin());
      }
      else pgoids[k-prev] = -1 ;
    }
    
    // Remove pg already flagged as match
    for (size_t kk = 0 ; kk < idx_to_remove2[i].size(); ++kk )
    {
      E_Int id_rm   = idx_to_remove2[i][kk];
      pgoids[id_rm] = -1;
    }
    
    prev = prev + pgoids.size();

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&pgoids[0], pgoids.size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PH history
    tpl = K_NUMPY::buildNumpyArray(&phoids2[i][0], phoids2[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // pushing out map z1loc_jz_to_ptl1 (>> Python dictionary)
  PyObject * pointList1_dict = PyDict_New();
  for (auto& i : z1loc_jz_to_ptl1)
  {

    E_Int zid = i.first;
    auto& jz_to_ptl = i.second;
    
    PyObject * titi_dict = PyDict_New();
      
    for (auto& j : jz_to_ptl)
    {
      E_Int jzid = j.first;
      auto& ptl = j.second;
	
      PyObject* key = Py_BuildValue("i", jzid);
      PyObject* np = K_NUMPY::buildNumpyArray(&ptl[0], ptl.size(), 1, 0);

      PyDict_SetItem(titi_dict, key, np);
      Py_DECREF(np);
    }

    PyObject* key = Py_BuildValue("i", zid);
    PyDict_SetItem(pointList1_dict, key, titi_dict);
    
  }
  PyList_Append(l, pointList1_dict);
  Py_DECREF(pointList1_dict);


  // pushing out map z2loc_jz_to_ptl2 (>> Python dictionary)
  PyObject * pointList2_dict = PyDict_New();
  for (auto& i : z2loc_jz_to_ptl2)
  {

    E_Int zid = i.first;
    auto& jz_to_ptl = i.second;
    
    PyObject * titi_dict = PyDict_New();
      
    for (auto& j : jz_to_ptl)
    {
      E_Int jzid = j.first;
      auto& ptl = j.second;
	
      PyObject* key = Py_BuildValue("i", jzid);
      PyObject* np = K_NUMPY::buildNumpyArray(&ptl[0], ptl.size(), 1, 0);

      PyDict_SetItem(titi_dict, key, np);
      Py_DECREF(np);
    }

    PyObject* key = Py_BuildValue("i", zid);
    PyDict_SetItem(pointList2_dict, key, titi_dict);
    
  }
  PyList_Append(l, pointList2_dict);
  Py_DECREF(pointList2_dict);


  for (E_Int i=0; i < nb_zones1; ++i)
  {
    delete crd1s[i];
    delete cnt1s[i];
  }

  for (E_Int i=0; i < nb_zones2; ++i)
  {
    delete crd2s[i];
    delete cnt2s[i];
  }
 
  return l;
}

//=============================================================================
/* Difference de 2 surfaces. */
//=============================================================================
PyObject* K_INTERSECTOR::booleanMinus(PyObject* self, PyObject* args)
{
  return call_operation(args, MINUS12);
}
//=============================================================================
/* Contour de l'intersection de 2 surfaces fermées. Le resultat est de 
   type BAR */
//=============================================================================
PyObject* K_INTERSECTOR::booleanIntersectionBorder(PyObject* self, PyObject* args)
{
  return call_xborder(args);
}

//=============================================================================
/* Modification d'un maillage volumique après résolution de l'intersection de sa peau 
   avec un autre maillage volumique.*/
//=============================================================================
PyObject* K_INTERSECTOR::booleanModifiedSolid(PyObject* self, PyObject* args)
{
  return call_operation(args, MODIFIED_SOLID);
}


PyObject* K_INTERSECTOR::DiffSurf(PyObject* self, PyObject* args)
{
  K_FLD::FloatArray pos1, pos2, pos;
  K_FLD::IntArray connect1, connect2, connect;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors;
  E_Int preserve_right, solid_right, agg_mode, itermax;
  bool improve_conformal_cloud_qual(false), outward_surf(true);
  
  bool ok = getArgs(args, DIFFSURF, pos1, connect1, pos2, connect2, tolerance, 
        preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, outward_surf, itermax, eltType, varString);
  if (!ok) return NULL;
  PyObject* tpl = NULL;
  E_Int err(0)/*, et=-1*/;

  char eltType2[20];
  strcpy(eltType2, eltType);

   if (strcmp(eltType, "NGON") != 0)
   {
      PyErr_SetString(PyExc_TypeError, "only NGON");
      return NULL; //fixme mem leak ?
   }

  {
    typedef NGON_BooleanOperator<K_FLD::FloatArray, K_FLD::IntArray> bo_t;
    
    bo_t::eAggregation AGG=bo_t::CONVEX;
    if (agg_mode == 0) // NONE
    {
      //AGG = bo_t::NONE; fixme: not working because of new_skin PG triangulation: disconnexion between modified and unlodified parts
    }
    else if (agg_mode == 2) AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual) BO.setConformizerParams(true);
    BO._outward = outward_surf;

    err = BO.Diffsurf(pos, connect);
    //et = 8;
  }
  
  if (!err)
    tpl = K_ARRAY::buildArray(pos, varString, connect, 8, "NGON", false);
  else
  {
    std::ostringstream o;
    o << "Diffsurf: failed to proceed.";
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
  }
  return tpl;
}

//=======================  Generator/booleanOperations.cpp ====================
