/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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

#include "Nuga/include/Triangulator.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/polygon.hxx"
#include "Nuga/include/polyhedron.hxx"
#include "Nuga/include/cdatastruct.hxx"
#include "Nuga/include/Tetrahedron.h"
#include "Nuga/include/Hexahedron.h"

using ngon_type = ngon_t<K_FLD::IntArray>;

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
  int                          nb_entities, i, nb_read;
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
      for (i = 0; i < dim; ++i)
        P[i] = fast_atof(words[i]);

      pos.pushBack(P, P+dim);
      ++nb_read;
      continue;
    }

    if ((curType != NONE) && (nb_read < nb_entities))
    {
      nods = nb_nodes[curType];
      for (i = 0; i < nods; ++i)
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
  int none = IDX_NONE;
  if (connects[HEX].size() != 0)
  {
    connect = connects[HEX];
    connects[TET].resize(8, connects[TET].cols(), &none);
    connect.pushBack(connects[TET]);
  }
  else if (connects[TET].size() != 0)
  {
    connect = connects[TET];
  }
  else if (connects[QUAD].size() != 0)
  {
    connect = connects[QUAD];
    connects[TRI].resize(4, connects[TRI].cols(), &none);
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

template <typename phmesh_type>
static E_Int read(const char* filename, phmesh_type& mesh)
{
  using INT_t = typename phmesh_type::INT_t;
  using FLT_t = typename phmesh_type::FLT_t;
  using cnt_t = K_FLD::DynArray<INT_t>;
  using ngon_type = ngon_t<cnt_t>;
  using crd_t = K_FLD::DynArray<FLT_t>;

  crd_t pos(mesh.CALLOC == 1);      // pass the allocation style
  cnt_t connect(mesh.CALLOC == 1);  // pass the allocation style

  E_Int ret = read(filename, pos, connect);

  if (pos.cols() == 0)
  {
    std::cout << "could not open file : " << filename << std::endl;
    return 1;
  }

  ngon_type ng;

  INT_t rows = connect.rows();

  if (rows == 4) // TETRA
    ngon_type::template convert<K_MESH::Tetrahedron>(connect, ng);
  else if (rows == 8) // HEXA
    ngon_type::template convert<K_MESH::Hexahedron>(connect, ng);
  else
  {
    std::cout << "elt type no handled" << std::endl;
    return 1;
  }
  
  ngon_type::clean_connectivity(ng, pos);

  // grabs coordinates
  int dim{ 3 };
  bool alloc;
  pos.relay_mem(mesh.crd.p, dim, mesh.crd.n, alloc);
  assert(int(alloc) == mesh.CALLOC);

  // grabs mesh
  ng.PGs.relay_mem(mesh.pgs);
  ng.PHs.relay_mem(mesh.phs);


  return ret;
}

  template< typename color_t = E_Int>
  static E_Int write(const char* fname, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const char* elt_type, const std::vector<bool>* keep = nullptr, const std::vector<color_t>* colors = nullptr)
  {
    if (crd.cols() == 0)
    return 0;

    std::ostringstream filename;
    filename << wdir << fname << ".mesh";

  FILE * file = fopen(filename.str().c_str(), "w");
  if (file == nullptr) return 1;

  E_Int  nb_pts(crd.cols()), COLS(cnt.cols()), dim(crd.rows());

  // Header
  fprintf(file, "MeshVersionFormatted 1\n");
  fprintf(file, "Dimension %i\n", dim);
       
  // Points
  fprintf(file, "Vertices\n");
  fprintf(file, "%i\n", nb_pts);
  
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

	fprintf(file, "%i\n", nb_elts);
  
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
			  case 2: fprintf(file, "%i %i 0\n", *(pC)+1, *(pC + 1) + 1); break;
			  case 3: fprintf(file, "%i %i %i 0\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1); break;
			  case 4: fprintf(file, "%i %i %i %i 0\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 3) + 1); break;
			  case 8: fprintf(file, "%i %i %i %i %i %i %i %i 0\n", *pC + 1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 2) + 1, *(pC + 3) + 1, *(pC + 4) + 1, *(pC + 5) + 1, *(pC + 5) + 1); break;
			  default:break;
			  }
		  }
		  else
		  {
			  switch (nb_nods)
			  {
			    case 2: fprintf(file, "%i %i %i\n", *(pC)+1, *(pC + 1) + 1, (E_Int)(*colors)[i]); break;
			    case 3: fprintf(file, "%i %i %i %i\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, (E_Int)(*colors)[i]); break;
			    case 4: fprintf(file, "%i %i %i %i %i\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 3) + 1, (E_Int)(*colors)[i]); break;
			    case 8: fprintf(file, "%i %i %i %i %i %i %i %i %i\n", *pC + 1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 2) + 1, *(pC + 3) + 1, *(pC + 4) + 1, *(pC + 5) + 1, *(pC + 5) + 1, (E_Int)(*colors)[i]); break;
			    default:break;
			  }
		  }
	  }
  }
    
  fprintf(file, "End\n");
  fclose(file);
  
  return 0;
  }

  template< typename color_t = E_Int>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<E_Int>* toprocess = nullptr, E_Int idx_start=0, const std::vector<color_t>* colors = nullptr)
  {
    std::string stype;
    if      (cnt.rows() == 2) stype="BAR";
    else if (cnt.rows() == 3) stype="TRI";
    else if (cnt.rows() == 4) stype="QUAD"; //fixme : TETRA

    std::vector<bool> keep;
    if (toprocess)
    {
      keep.resize(cnt.cols(), false);

      for (size_t u = 0; u< toprocess->size(); ++u)
        keep[(*toprocess)[u] - idx_start] = true;
    }

    std::vector<E_Int> nids;
    K_FLD::FloatArray tmpcrd(crd);
    K_FLD::IntArray tmpcnt(cnt);
    NUGA::MeshTool::compact_to_mesh(tmpcrd, tmpcnt, nids);

    return write<color_t>(filename, tmpcrd, tmpcnt, stype.c_str(), keep.empty() ? nullptr : &keep, colors);
  }

  template< typename color_t = E_Int>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_unit& pgs, const std::vector<E_Int>* toprocess = nullptr, E_Int idx_start = 0, const std::vector<color_t>* colors = nullptr)
  {
    
    K_FLD::IntArray cT;
    DELAUNAY::Triangulator dt;
    std::vector<E_Int> Tcolors;

    E_Int pureBasic = 0;

    if (toprocess)
    {
      for (E_Int i = 0; i < toprocess->size(); ++i)
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
          K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, crd, pgs.get_facets_ptr(PGi), pgs.stride(PGi), 1, cT);
          if (colors)Tcolors.resize(cT.cols(), (E_Int)(*colors)[PGi]);
        }
      }
      else
      {
        for (E_Int i = 0; i < pgs.size(); ++i)
        {
          K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, crd, pgs.get_facets_ptr(i), pgs.stride(i), 1, cT);
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
      std::vector<color_t> newPGcolor;
      const std::vector<color_t> *pColors(colors);

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

  template <typename crd3D_t, typename vngon_unit>
  static E_Int write(const char* filename, crd3D_t& crd3D, const vngon_unit& pgs)
  {
    K_FLD::FloatArray crd(crd3D.p, 3, crd3D.n, (crd3D.CALLOC == 1));

    ngon_unit PGS(pgs.elts, pgs.range, pgs.nrange);

    write<E_Int>(filename, crd, PGS);

    //hack because crd grabbed memory
    int dim{ 3 };
    bool calloc{ false };
    crd.relay_mem(crd3D.p, dim, crd3D.n, calloc);

    return 0;
  }

  template< typename color_t = E_Int>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_type& ng, const std::vector<E_Int>* toprocess = nullptr, E_Int idx_start = 0, const std::vector<color_t>* colors = nullptr)
  {
    K_FLD::IntArray cT3;
    DELAUNAY::Triangulator dt;
    std::vector<color_t> PHcolor;
    if (colors) PHcolor.insert(PHcolor.end(), ALL(*colors));

    ngon_type ngt(ng);

    if (toprocess)
    {
      ngon_unit phs;
      Vector_t<E_Int> oids;
      std::vector<E_Int> shiftedToPoc(ALL(*toprocess));
      K_CONNECT::IdTool::shift(shiftedToPoc, -idx_start);
      ng.PHs.extract(shiftedToPoc, phs, oids);
      ngt.PHs = phs;
      //update PH color
      if (colors)
      {
        std::vector<color_t> newPHcolor(ngt.PHs.size(), 0);
        for (size_t k = 0; k < newPHcolor.size(); ++k)
          newPHcolor[k] = PHcolor[oids[k]];
        PHcolor = newPHcolor;
      }
    }
    
    // clean up with unrelevant pgs
    std::vector<E_Int> pgnids, phnids;
    ngt.remove_unreferenced_pgs(pgnids, phnids);
    K_CONNECT::IdTool::compact(PHcolor, phnids);

    std::vector<E_Int> PGcolor;
    if (colors)
    {
      PGcolor.resize(ngt.PGs.size(), 0);
      for (size_t i = 0; i < ngt.PHs.size(); ++i)
      {
        const E_Int* PGi = ngt.PHs.get_facets_ptr(i);
        E_Int stride = ngt.PHs.stride(i);
        for (size_t j = 0; j < stride; ++j)
          if (PGcolor[PGi[j] - 1]==0) PGcolor[PGi[j] - 1] = PHcolor[i];
      }
    }
    
    return medith::write(filename, crd, ngt.PGs, nullptr, 0,colors ? &PGcolor : nullptr);
  }
  
  template <typename aELT>
  static E_Int write(const char* filename, const aELT& ae);

  template <typename ngo_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngo_t& ng)
  {    
    ngo_t ngt(ng);
    std::vector<E_Int> pgnids, phnids;
    ngt.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdl(crd);
    ngo_t::compact_to_used_nodes(ngt.PGs, crdl);
    
    write(filename, crdl, ngt.PGs);
    
    return 0;
  }

  template <typename ngo_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngo_t& ng, E_Int i)
  {
    ng.PGs.updateFacets();
    ng.PHs.updateFacets();

    ngon_unit ph;
    ph.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));

    ngo_t one_ph(ng.PGs, ph);
    Vector_t<E_Int> pgnids, phnids;
    one_ph.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdl(crd);
    ngo_t::compact_to_used_nodes(one_ph.PGs, crdl);

    write<E_Int>(filename, crdl, one_ph.PGs);

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



#endif
