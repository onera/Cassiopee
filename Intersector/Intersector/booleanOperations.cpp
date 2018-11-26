/*    
    Copyright 2013-2018 Onera.

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

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/Boolean/TRI_BooleanOperator.h"
# include "Nuga/Boolean/BAR_BooleanOperator.h"
# include "Nuga/Boolean/NGON_BooleanOperator.h"
#include <memory>

using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
E_Int chrono::verbose=1;
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
bool is_valid (char* eltType, eOperation oper)
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
bool getArgs(PyObject* args, eOperation oper,
             K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
             K_FLD::FloatArray& pos2, K_FLD::IntArray& connect2,
             E_Float& tolerance, E_Int& preserve_right, E_Int& solid_right, E_Int& agg_mode, bool& improve_conformal_cloud_qual, char*& eltType, char*& varString)
{
  PyObject *arrS[2];
  E_Float tol = 0.;
  E_Int preserv_r, solid_r, imp_qual(0);
  std::ostringstream o;
  std::string opername = getName(oper);

  if (!PYPARSETUPLE(args, 
                    "OOdllll", "OOdiiii", "OOfllll", "OOfiiii", 
                    &arrS[0], &arrS[1], &tol, &preserv_r, &solid_r, &agg_mode, &imp_qual))
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

bool getUnionArgs(PyObject* args,
             K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
             K_FLD::FloatArray& pos2, K_FLD::IntArray& connect2,
             E_Float& tolerance, E_Int& preserve_right, E_Int& solid_right, E_Int& agg_mode, bool& improve_conformal_cloud_qual,
             std::vector<E_Int>& pgsList,
             char*& eltType, char*& varString)
{
  PyObject *arrS[2], *pgs;
  E_Float tol = 0.;
  E_Int preserv_r, solid_r, imp_qual(0);
  std::ostringstream o;
  std::string opername = getName(UNION);

  pgsList.clear();

  if (!PYPARSETUPLE(args, 
                    "OOdllllO", "OOdiiiiO", "OOfllllO", "OOfiiiiO", 
                    &arrS[0], &arrS[1], &tol, &preserv_r, &solid_r, &agg_mode, &imp_qual, &pgs))
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
      res = K_NUMPY::getFromNumpyArray(pgs, inds, true);

    std::auto_ptr<K_FLD::FldArrayI> pL(inds); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY is : " << res << std::endl;
    if ((res == 1) && (inds != NULL)  && (inds->getSize() != 0))
    {
      E_Int nb_ghost_pgs = inds->getSize();
      //E_Int minid(INT_MAX), maxid(-1);
      pgsList.resize(nb_ghost_pgs);
      for (size_t i = 0; i < nb_ghost_pgs; ++i) 
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
PyObject* call_operation(PyObject* args, eOperation oper)
{
  K_FLD::FloatArray pos1, pos2, pos;
  K_FLD::IntArray connect1, connect2, connect;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors;
  E_Int preserve_right, solid_right, agg_mode;
  bool improve_conformal_cloud_qual(false);
  
  bool ok = getArgs(args, oper, pos1, connect1, pos2, connect2, tolerance, 
		    preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, eltType, varString);
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
      BO = new TRI_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
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
    else if (agg_mode == 2)
      AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual) BO.setConformizerParams(true);

    switch (oper)
    {
    case INTERSECTION:        
      err=BO.Intersection(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case UNION :              
      err=BO.Union(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case MINUS12:             
      err=BO.Diff(pos, connect, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right); break;
    case MODIFIED_SOLID:             
      err=BO.Modified_Solid(pos, connect, (bo_t::eMergePolicy)preserve_right); break;
    case DIFFSURF:             
      err=BO.Diffsurf(pos, connect); break;
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
  }
  return tpl;
}

PyObject* call_union(PyObject* args)
{
  K_FLD::FloatArray pos1, pos2;
  K_FLD::IntArray connect1, connect2;
  E_Float tolerance;
  char *eltType, *varString;
  std::vector<E_Int> colors, ghost_pgs;
  E_Int preserve_right, solid_right, agg_mode;
  bool improve_conformal_cloud_qual(false);
  
  bool ok = getUnionArgs(args, pos1, connect1, pos2, connect2, tolerance, 
        preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, ghost_pgs, eltType, varString);
  if (!ok) return NULL;
  PyObject* tpl = NULL;
  E_Int err(0), et=-1;
  PyObject *l(PyList_New(0));

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
      BO = new TRI_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
      op = getOperation<TRI_BooleanOperator>(UNION);
    }
    else if (strcmp(eltType, "BAR") == 0)
    {
      BO = new BAR_BooleanOperator (pos1, connect1, pos2, connect2, tolerance);
      op = getOperation<BAR_BooleanOperator>(UNION);
    }
    
    crds.resize(1);
    cnts.resize(1);

    err = (BO->*op)(crds[0], cnts[0], colors);
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
    else if (agg_mode == 2)
      AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual) BO.setConformizerParams(true);

    if (!ghost_pgs.empty())
      BO.passPGs(1, ghost_pgs);

    //BO.setTriangulatorParams(/*do not shuffle : */true, /*improve quality:*/false);

    K_FLD::IntArray cnt;
    K_FLD::FloatArray crd;
    err=BO.Union(crd, cnt, (bo_t::eInterPolicy)solid_right, (bo_t::eMergePolicy)preserve_right);
    et = 8;

#ifndef DEBUG_W_PYTHON_LAYER
    crds.push_back(crd);
    cnts.push_back(cnt);
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

  }
  
#ifndef DEBUG_W_PYTHON_LAYER
  {
  	if (err != 1)
  	  tpl = K_ARRAY::buildArray(crds[0], varString, cnts[0], et, eltType2, false);
  	else
    {
      std::ostringstream o;
      o << "Union : failed to proceed.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
    }
  	
  	return tpl;
  }
#else
  {
    {
      size_t nb_zones = cnts.size();
      //std::cout << "nb zones : " << nb_zones << std::endl;
      for (size_t i=0; i < nb_zones; ++i)
      {
        if (crds[i].cols() == 0) continue;

      	tpl = K_ARRAY::buildArray(crds[i], varString, cnts[i], et, "NGON", false);

      	PyObject* o = Py_BuildValue("(Os)", tpl, znames[i].c_str());
        PyList_Append(l, o);
        Py_DECREF(tpl);
      }
    }
    if (err != 1)
    {
      std::ostringstream o;
      o << "Union : failed to proceed.";
      PyErr_SetString(PyExc_TypeError, o.str().c_str());
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
  return call_operation(args, INTERSECTION_BORDER);
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
  E_Int preserve_right, solid_right, agg_mode;
  bool improve_conformal_cloud_qual(false);
  
  bool ok = getArgs(args, DIFFSURF, pos1, connect1, pos2, connect2, tolerance, 
        preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, eltType, varString);
  if (!ok) return NULL;
  PyObject* tpl = NULL;
  E_Int err(0), et=-1;

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
      //AGG = bo_t::NONE; fixme : not working because of new_skin PG triangulation : diconnexion between modified and unlodified parts
    }
    else if (agg_mode == 2)
      AGG = bo_t::FULL;
        
    bo_t BO(pos1, connect1, pos2, connect2, tolerance, AGG);

    if (improve_conformal_cloud_qual)BO.setConformizerParams(true);

    err=BO.Diffsurf(pos, connect);
    et = 8;
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
