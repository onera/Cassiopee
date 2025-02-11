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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef MEDITH_HXX
#define MEDITH_HXX

#ifdef VISUAL
#pragma warning (disable : 4996)
#endif

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <cstring>

#include "Nuga/include/DynArray.h"
#include "Nuga/include/MeshTool.h"


class medith
{
public:
  enum eType { NONE = -2, VERT = -1, EDGE, TRI, QUAD, TET, HEX};
  enum eColor { COL_NONE = 0, RED = 1, GREEN = 2, YELLOW = 3};

  static std::string wdir;

public:
  static E_Int read (const char* filename, K_FLD::FloatArray& pos, K_FLD::IntArray& connect)
  {
  std::string                  line, entity;
  char cline[512];
  char* words[10];
  std::ifstream                file (filename);
  FILE * fp = fopen(filename, "r");
  int                          nb_entities, nb_read;
  std::vector<K_FLD::IntArray> connects;
  double                       P[3];
  int                          nb_nodes[5], S[8], nods, dim(3);

  pos.clear();
  connect.clear();

  if (fp == NULL)
    return 1;

  nb_nodes[EDGE] = 2;
  nb_nodes[TRI] = 3;
  nb_nodes[QUAD] = nb_nodes[TET] = 4;
  nb_nodes[HEX] = 8;

  eType curType(NONE);

  connects.resize(5);
  for (size_t i = 0; i < 5; ++i)
    connects[i].set_alloc(connect);

  /*keys.insert("Dimension\n");
  keys.insert("Vertices\n");
  keys.insert("Edges\n");
  keys.insert("Triangles\n");
  keys.insert("Quadrilaterals\n");
  keys.insert("Tetrahedra\n");
  keys.insert("Hexahedra\n");
  keys.insert("End\n");*/

  curType = NONE;
  
  while (fgets(cline, 512, fp) != NULL)
  {
    line = cline;
    if (line.empty())
      continue;
            
    get_words(line, ' ', &words[0]);

    if ((curType == VERT) && (nb_read < nb_entities))
    {
      for (int i = 0; i < dim; ++i)
        P[i] = fast_atof(words[i]);

      pos.pushBack(P, P+dim);
      ++nb_read;
      continue;
    }

    if ((curType != NONE) && (nb_read < nb_entities))
    {
      nods = nb_nodes[curType];
      for (int i = 0; i < nods; ++i)
        S[i] = fast_atoindex(words[i])-1;

      connects[curType].pushBack(S, S+nods);
      ++nb_read;
      continue;
    }
    
    entity = words[0];
    curType = NONE;

    if (entity == "End\n")
        break;
    if (entity == "Dimension\n")
    {
      fgets(cline, 512, fp);
      dim = fast_atoindex(cline);
		  continue;
    }
    if ( (entity == "Vertices\n") || (entity == "Edges\n") || (entity == "Triangles\n") || (entity == "Quadrilaterals\n") || (entity == "Tetrahedra\n") || (entity == "Hexahedra\n") )
    {      
      fgets(cline, 512, fp);
      line = cline;
      get_words(line, ' ', words);

      nb_entities = atoi(words[0]);
      nb_read=0;

      if (entity == "Vertices\n")
      {
        pos.reserve(dim, nb_entities);
        curType = VERT;
      }
      else if (entity == "Edges\n")
      {
        connects[EDGE].reserve(2, connects[EDGE].size() + nb_entities);
        curType = EDGE;
      }
      else if (entity == "Triangles\n")
      {
        connects[TRI].reserve(3, connects[TRI].size() + nb_entities);
        curType = TRI;
      }
      else if (entity == "Quadrilaterals\n")
      {
        connects[QUAD].reserve(4, connects[QUAD].size() + nb_entities);
        curType = QUAD;
      }
      else if (entity == "Tetrahedra\n")
      {
        connects[TET].reserve(4, connects[TET].size() + nb_entities);
        curType = TET;
      }
      else if (entity == "Hexahedra\n")
      {
        connects[HEX].reserve(8, connects[HEX].size() + nb_entities);
        curType = HEX;
      }
      continue;
    }
  }

  // Gather the meshes.
  E_Int none = IDX_NONE;
  if (connects[HEX].size() != 0)
  {
    connect = connects[HEX];
    connects[TET].resize((E_Int)8, connects[TET].cols(), &none);
    connect.pushBack(connects[TET]);
  }
  else if (connects[TET].size() != 0)
  {
    connect = connects[TET];
  }
  else if (connects[QUAD].size() != 0)
  {
    connect = connects[QUAD];
    connects[TRI].resize((E_Int)4, connects[TRI].cols(), &none);
    connect.pushBack(connects[TRI]);
  }
  else if (connects[TRI].size() != 0)
  {
    connect = connects[TRI];
  }
  else if (connects[EDGE].size() != 0)
  {
    connect = connects[EDGE];
  }

  return 0;
}

template <typename T1, typename T2>
static E_Int read(const char* filename, T1*& pcrd, T2& dim, T2& npts, bool& calloc1, T2*& pcnt, T2& nnodes, T2& nelts, bool& calloc2)
{
  K_FLD::DynArray<T1> pos;
  K_FLD::DynArray<T2> connect;

  E_Int ret = read(filename, pos, connect);

  pos.relay_mem(pcrd, dim, npts, calloc1);
  connect.relay_mem(pcnt, nnodes, nelts, calloc2);

  return ret;
}

  static E_Int write(const char* fname, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const char* elt_type, const std::vector<bool>* keep = nullptr, const std::vector<E_Int>* colors = nullptr)
  {
    if (crd.cols() == 0)
    return 0;

    std::ostringstream filename;
    filename << wdir << fname << ".mesh";

  FILE * file = fopen(filename.str().c_str(), "w");
  if (file == nullptr) return 1;

  E_Int nb_pts(crd.cols()), COLS(cnt.cols()), dim(crd.rows());

  // Header
  fprintf(file, "MeshVersionFormatted 1\n");
  fprintf(file, "Dimension " SF_D_"\n", dim);
       
  // Points
  fprintf(file, "Vertices\n");
  fprintf(file, SF_D_"\n", nb_pts);
  
  const E_Float* pP;

  if (dim == 3)
  {
    for (E_Int i = 0; i < nb_pts; ++i)
    {
      pP = crd.col(i);
	    fprintf(file, "%.20f %.20f %.20f 0\n", *(pP), *(pP+1), *(pP+2));
    }
  }
  else if (dim ==2)
  {
    for (E_Int i = 0; i < nb_pts; ++i)
    {
      pP = crd.col(i);
      fprintf(file, "%.20f %.20f 0\n", *(pP), *(pP + 1));
    }
  }
  
  fprintf(file, "\n");

  if (COLS==0)
  {
    fprintf(file, "End\n");
    fclose(file);
  }

  // Connectivity.
  E_Int nb_nods(cnt.rows()), nb_elts(COLS);

  std::string et;
  if (elt_type) et = elt_type;

  if (nb_nods == 2) // Edges
    fprintf(file, "Edges\n");
  else if (nb_nods == 3) // Triangles
    fprintf(file, "Triangles\n");
  else if ((nb_nods == 4) && (et.find("QUAD") != std::string::npos)) // Q4
    fprintf(file, "Quadrilaterals\n");
  else if ((nb_nods == 4) && (et.find("TETRA") != std::string::npos)) // TH4
    fprintf(file, "Tetrahedra\n");
  else if (nb_nods == 8) // Hexahedra
    fprintf(file, "Hexaedra\n");
  else //Unknonw
    return 1;

  if (keep)
  {
    for (E_Int i = 0; i < COLS; ++i)
    {
      if ((*keep)[i] == false)
        --nb_elts;
    }
  }

  if (nb_elts == 0)
  {
    fprintf(file, "End\n");
    fclose(file);
    return 0;
  }

	fprintf(file, SF_D_"\n", nb_elts);
  
  bool valid;
  const E_Int* pC;
  for (E_Int i = 0; i < COLS; ++i)
  {
    valid = true;

    if (keep && (*keep)[i] == false)
      continue;

	  pC = cnt.col(i);
	  for (E_Int k = 0; (k < nb_nods) && valid; ++k)
	    valid = (*(pC + k) < nb_pts) && (*(pC + k) >= 0);

	  if (valid)
	  {
		  if (!colors)
		  {
			  switch (nb_nods)
			  {
			  case 2: fprintf(file, SF_D2_" 0\n", *(pC)+1, *(pC + 1) + 1); break;
			  case 3: fprintf(file, SF_D3_" 0\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1); break;
			  case 4: fprintf(file, SF_D4_" 0\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 3) + 1); break;
			  case 8: fprintf(file, SF_D8_" 0\n", *pC + 1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 2) + 1, *(pC + 3) + 1, *(pC + 4) + 1, *(pC + 5) + 1, *(pC + 5) + 1); break;
			  default:break;
			  }
		  }
		  else
		  {
			  switch (nb_nods)
			  {
			    case 2: fprintf(file, SF_D3_"\n", *(pC)+1, *(pC + 1) + 1, (E_Int)(*colors)[i]); break;
			    case 3: fprintf(file, SF_D4_"\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, (E_Int)(*colors)[i]); break;
			    case 4: fprintf(file, SF_D5_"\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 3) + 1, (E_Int)(*colors)[i]); break;
			    case 8: fprintf(file, SF_D9_"\n", *pC + 1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 2) + 1, *(pC + 3) + 1, *(pC + 4) + 1, *(pC + 5) + 1, *(pC + 5) + 1, (E_Int)(*colors)[i]); break;
			    default:break;
			  }
		  }
	  }
  }
    
  fprintf(file, "End\n");
  fclose(file);
  
  return 0;
  }

  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<E_Int>* toprocess = nullptr, E_Int idx_start=0, const std::vector<E_Int>* colors = nullptr)
  {
    std::string stype;
    if      (cnt.rows() == 2) stype="BAR";
    else if (cnt.rows() == 3) stype="TRI";
    else if (cnt.rows() == 4) stype="QUAD"; //fixme : TETRA

    std::vector<bool> keep;
    if (toprocess)
    {
      keep.resize(cnt.cols(), false);

      for (size_t u = 0; u < toprocess->size(); ++u)
        keep[(*toprocess)[u] - idx_start] = true;
    }

    std::vector<E_Int> nids;
    K_FLD::FloatArray tmpcrd(crd);
    K_FLD::IntArray tmpcnt(cnt);
    NUGA::MeshTool::compact_to_mesh(tmpcrd, tmpcnt, nids);

    return write(filename, tmpcrd, tmpcnt, stype.c_str(), keep.empty() ? nullptr : &keep, colors);
  }

  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_unit& pgs, const std::vector<E_Int>* toprocess = nullptr, E_Int idx_start = 0, const std::vector<E_Int>* colors = nullptr)
  {

    K_FLD::IntArray cT;
    std::vector<E_Int> Tcolors;

    E_Int pureBasic = 0;

    if (toprocess)
    {
      for (size_t i = 0; i < toprocess->size(); ++i)
      {
        E_Int PGi = (*toprocess)[i] - idx_start;
        E_Int str = pgs.stride(PGi);
        if (pureBasic == 0)
        {
          if (str != 3 && str != 4) break;
          pureBasic = str;
        }
        else if (pureBasic != str)
        {
          pureBasic = 0;
          break;
        }
      }
    }
    else
    {
      for (E_Int i = 0; i < pgs.size(); ++i)
      {
        E_Int str = pgs.stride(i);
        if (pureBasic == 0)
        {
          if (str != 3 && str != 4) break;
          pureBasic = str;
        }
        else if (pureBasic != str)
        {
          pureBasic = 0;
          break;
        }
      }
    }

    if (!pureBasic)
    {
      if (toprocess)
      {
        for (size_t i = 0; i < toprocess->size(); ++i)
        {
          E_Int PGi = (*toprocess)[i] - idx_start;
          //K_MESH::Polygon::triangulate(crd, pgs.get_facets_ptr(PGi), pgs.stride(PGi), 1, cT);
          if (colors)Tcolors.resize(cT.cols(), (E_Int)(*colors)[PGi]);
        }
      }
      else
      {
        for (E_Int i = 0; i < pgs.size(); ++i)
        {
          //K_MESH::Polygon::triangulate(crd, pgs.get_facets_ptr(i), pgs.stride(i), 1, cT);
          if (colors)Tcolors.resize(cT.cols(), (E_Int)(*colors)[i]);
        }
      }

      std::vector<E_Int> nids;
      K_FLD::FloatArray tmpcrd(crd);
      K_FLD::IntArray tmpcnt(cT);
      NUGA::MeshTool::compact_to_mesh(tmpcrd, tmpcnt, nids);

      return write(filename, tmpcrd, tmpcnt, "TRI", nullptr, colors ? &Tcolors : nullptr);
    }
    else
    {
      std::vector<E_Int> newPGcolor;
      const std::vector<E_Int> *pColors(colors);

      ngon_unit pgs2(pgs);

      if (toprocess)
      {
        Vector_t<E_Int> oids;
        std::vector<E_Int> shiftedToPoc(ALL(*toprocess));
        K_CONNECT::IdTool::shift(shiftedToPoc, -idx_start);
        pgs.extract(shiftedToPoc, pgs2, oids);
        //update PG color
        if (colors)
        {
          newPGcolor.resize(pgs.size(), 0);
          for (size_t k = 0; k < newPGcolor.size(); ++k)
            newPGcolor[k] = (*colors)[oids[k]];
          pColors = &newPGcolor;
        }
      }
      ngon_unit::convert_ngon_unit_to_fixed_stride(pgs2, 1, cT);

      std::vector<E_Int> nids;
      K_FLD::FloatArray tmpcrd(crd);
      K_FLD::IntArray tmpcnt(cT);
      NUGA::MeshTool::compact_to_mesh(tmpcrd, tmpcnt, nids);

      if (cT.rows() == 3) return write(filename, tmpcrd, tmpcnt, "TRI", nullptr, pColors);
      else return write(filename, tmpcrd, tmpcnt, "QUAD", nullptr, pColors);
    }
  }

  ///  
  template <typename aELT>
  static E_Int write(const char* filename, const aELT& ae)
  {
    return write(filename, ae.m_crd, ae.m_pgs);
  }

  template <typename ngo_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngo_t& ng)
  {    
    ngo_t ngt(ng);
    std::vector<E_Int> pgnids, phnids;
    ngt.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdl(crd);
    ngo_t::compact_to_used_nodes(ngt.PGs, crdl);
    
    return write(filename, crdl, ngt.PGs);
  }

  template <typename ngo_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngo_t& ng, E_Int i)
  {
    if (i >= ng.PHs.size())
      return 1;

    ng.PGs.updateFacets();
    ng.PHs.updateFacets();

    ngon_unit ph;
    ph.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));

    ngo_t one_ph(ng.PGs, ph);
    Vector_t<E_Int> pgnids, phnids;
    one_ph.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdl(crd);
    ngo_t::compact_to_used_nodes(one_ph.PGs, crdl);

    write(filename, crdl, one_ph.PGs);

    return 0;
  }

  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Int i)
  {
    PGs.updateFacets();
    
    write(filename, crd, PGs.get_facets_ptr(i), PGs.stride(i), 1);  

    return 0;
  }
  
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int idx_start)
  {
    K_FLD::IntArray cnt(2,nb_nodes);
    for (E_Int i=0; i < nb_nodes; ++i)
    {
      cnt(0,i) = nodes[i] - idx_start;
      cnt(1,i) = nodes[(i+1)%nb_nodes] - idx_start;
    }
    return write(filename, crd, cnt);
  }

  ///
  static void draw_wired_PG(const char* fname, const K_FLD::FloatArray& coord, const ngon_unit& PGs, E_Int ith, E_Float *normal)
  {
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(coord);
    K_FLD::IntArray connectE;
    //E_Int n0, n1;
    E_Float P0[3], P1[3], Lmin(NUGA::FLOAT_MAX), L2;

    E_Int nb_nodes, E[2];
    const E_Int* pNi = PGs.get_facets_ptr(ith);
    nb_nodes = PGs.stride(ith);

    for (E_Int j = 0; j < nb_nodes; ++j)
    {
      E[0] = *(pNi + j) - 1;
      E[1] = *(pNi + (j + 1) % nb_nodes) - 1;
      connectE.pushBack(E, E + 2);
      L2 = NUGA::sqrDistance(coord.col(E[0]), coord.col(E[1]), 3);
      Lmin = (L2 < Lmin) ? L2 : Lmin;
    }

    Lmin = 0.5*::sqrt(Lmin);

    K_MESH::Polygon::iso_barycenter<acrd_t, 3 >(acrd, pNi, nb_nodes, 1, P0);

    K_FLD::FloatArray crd(coord);

    E_Float Norm[3];
    if (normal == nullptr)
    {
      K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNi, nb_nodes, 1, Norm);
      normal = Norm;
    }

    NUGA::sum<3>(1., P0, Lmin, normal, P1);
    crd.pushBack(P0, P0 + 3);
    E[0] = crd.cols() - 1;
    crd.pushBack(P1, P1 + 3);
    E[1] = crd.cols() - 1;
    connectE.pushBack(E, E + 2);

    medith::write(fname, crd, connectE, "BAR");
  }
  
  
  ///
  static void
    get_words
    (const std::string& str_line, char delim, char** oWords)
  {
    //oWords.clear();
    oWords[0] = oWords[1] = oWords[2] = oWords[3] = oWords[4] = oWords[5] = oWords[6] = oWords[7] = oWords[8] = oWords[9] = { 0 };
    
    char* buf = const_cast<char*>(str_line.c_str());
    char* pch = strtok(buf," ");
    int c=0;
    while (pch != NULL)
    {
      //oWords.push_back(pch);
      oWords[c++]=pch;
      pch = strtok (NULL, " ");
    }
  }
  
  #define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')

  static double fast_atof(const char *p)
  {
    int frac;
    double sign, value, scale;

    // Skip leading white space, if any.

    while (white_space(*p)) {
      p += 1;
    }

    // Get sign, if any.

    sign = 1.0;
    if (*p == '-') {
      sign = -1.0;
      p += 1;

    }
    else if (*p == '+') {
      p += 1;
    }

    // Get digits before decimal point or exponent, if any.

    for (value = 0.0; valid_digit(*p); p += 1) {
      value = value * 10.0 + (*p - '0');
    }

    // Get digits after decimal point, if any.

    if (*p == '.') {
      double pow10 = 10.0;
      p += 1;
      while (valid_digit(*p)) {
        value += (*p - '0') / pow10;
        pow10 *= 10.0;
        p += 1;
      }
    }

    // Handle exponent, if any.

    frac = 0;
    scale = 1.0;
    if ((*p == 'e') || (*p == 'E')) {
      unsigned int expon;

      // Get sign of exponent, if any.

      p += 1;
      if (*p == '-') {
        frac = 1;
        p += 1;

      }
      else if (*p == '+') {
        p += 1;
      }

      // Get digits of exponent, if any.

      for (expon = 0; valid_digit(*p); p += 1) {
        expon = expon * 10 + (*p - '0');
      }
      if (expon > 308) expon = 308;

      // Calculate scaling factor.

      while (expon >= 50) { scale *= 1E50; expon -= 50; }
      while (expon >= 8) { scale *= 1E8;  expon -= 8; }
      while (expon >   0) { scale *= 10.0; expon -= 1; }
    }

    // Return signed and scaled floating point result.

    return sign * (frac ? (value / scale) : (value * scale));
  }

  static int fast_atoindex(const char *p)
  {
    int value;

    // Get digits before decimal point or exponent, if any.

    for (value = 0; valid_digit(*p); p += 1) {
      value = value * 10 + (*p - '0');
    }
    return value;
  }
  
private:
  medith(void){}
  ~medith(void){}
};
/*
template <> inline
E_Int medith::write<NUGA::aPolygon>(const char* filename, const NUGA::aPolygon& ae)
{
  return write(filename, ae.m_crd, &ae.m_nodes[0], ae.m_nodes.size(), 0);
}

template <> inline
E_Int medith::write<NUGA::aPolyhedron<0>>(const char* filename, const NUGA::aPolyhedron<0>& ae)
{
  return write(filename, ae.m_crd, ae.m_pgs);

}
*/





#endif
