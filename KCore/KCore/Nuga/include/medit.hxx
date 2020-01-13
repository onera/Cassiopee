

#ifndef MEDITH_HXX
#define MEDITH_HXX

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <cstring>

#include "Nuga/Delaunay/Triangulator.h"
#include "Fld/ngon_t.hxx"
#include "MeshElement/Polygon.h"

using ngon_type = ngon_t<K_FLD::IntArray>;

class medith
{
public:
  enum eType { NONE = -2, VERT = -1, EDGE, TRI, QUAD, TET, HEX};
  enum eColor { COL_NONE = 0, RED = 1, GREEN = 2, YELLOW = 3};

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
  int none = E_IDX_NONE;
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
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const char* elt_type, const std::vector<bool>* keep = 0, const std::vector<E_Int>* colors = 0)
  {
    if (crd.cols() == 0)
    return 0;
  FILE * file = fopen(filename, "w");
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
			    case 2: fprintf(file, "%i %i %i\n", *(pC)+1, *(pC + 1) + 1, (*colors)[i]); break;
			    case 3: fprintf(file, "%i %i %i %i\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, (*colors)[i]); break;
			    case 4: fprintf(file, "%i %i %i %i %i\n", *(pC)+1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 3) + 1, (*colors)[i]); break;
			    case 8: fprintf(file, "%i %i %i %i %i %i %i %i %i\n", *pC + 1, *(pC + 1) + 1, *(pC + 2) + 1, *(pC + 2) + 1, *(pC + 3) + 1, *(pC + 4) + 1, *(pC + 5) + 1, *(pC + 5) + 1, (*colors)[i]); break;
			    default:break;
			  }
		  }
	  }
  }
    
  fprintf(file, "End\n");
  fclose(file);
  
  return 0;
  }
  
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_unit& pgs)
  {
    K_FLD::IntArray cT3;
    DELAUNAY::Triangulator dt;
    for (E_Int i=0; i < pgs.size(); ++i)
    {
      K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, crd, pgs.get_facets_ptr(i), pgs.stride(i), 1, cT3);
    }
    
    write(filename, crd, cT3, "TRI");
    
    return 0;
  }
  
  template <typename ngon_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngon_t& ng, const std::vector<E_Int>& toprocess)
  {
    ngon_unit phs;
    Vector_t<E_Int> oids;
    ng.PHs.extract(toprocess, phs, oids);
    
    ngon_t ngt(ng.PGs, phs);
    std::vector<E_Int> pgnids, phnids;
    ngt.remove_unreferenced_pgs(pgnids, phnids);
    
    write(filename, crd, ngt.PGs);
    
    return 0;
  }
  
  template <typename ngo_t>
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const ngo_t& ng)
  {    
    ngo_t ngt(ng);
    std::vector<E_Int> pgnids, phnids;
    ngt.remove_unreferenced_pgs(pgnids, phnids);
    
    write(filename, crd, ngt.PGs);
    
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

    write(filename, crdl, one_ph.PGs);

    return 0;
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


#endif
