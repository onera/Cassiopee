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
//Authors : Sam Landier (sam.landier@onera.fr)

#include "Nuga/include/iodata.h"
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace DELAUNAY
{

iodata::iodata(void)
{
}

iodata::~iodata(void)
{
}

void
    __get_words
    (const std::string& str_line, char delim, std::vector<std::string>& oWords)
  {
    std::stringstream    sline (str_line);
    std::string          word;

    oWords.clear();

    do
    {
      word.clear();
      std::getline (sline, word, sline.widen(delim)); 
      if (!word.empty())
        oWords.push_back (word);
    }
    while (!word.empty());
  }

///
template <typename T>
void write_vector(std::ofstream& file, const std::vector<T>& vals)
{
  if (vals.size() == 0)
    return;

  file << vals.size() << std::endl;
  for (size_t i = 0; i < vals.size(); ++i)
    file << vals[i] << std::endl;
  file << std::endl;
}

///
template <typename T>
void write_matrix(std::ofstream& file, const K_FLD::DynArray<T>& vals)
{
  if (vals.cols() == 0)
    return;

  file << vals.rows() << " " << vals.cols() << std::endl;
  for (E_Int j = 0; j < vals.cols(); ++j)
  {
    for (E_Int i = 0; i < vals.rows(); ++i)
    {
      file << vals(i,j) << " ";
    }
    file << std::endl;
  }
  file << std::endl;
}

E_Int iodata::write (const char* filename, const DELAUNAY::MeshData& data)
{

  if (data.pos->cols() == 0)
    return 0;
  std::ofstream                file (filename);

  // Header
  file << "Mesh Data Format"  << std::endl;
       
  // POS
  file.flags (std::ios::scientific);
  file.precision (10);
  file << "POS"  << std::endl;
  write_matrix(file, *data.pos);

  // CONNECTB
  if (data.connectB->cols())
  {
    file << "CONNECTB"  << std::endl;
    write_matrix(file, *data.connectB);
  }

  // HARD NODES
  if (!data.hardNodes.empty())
  {
    file << "HARDNODES"  << std::endl;
    write_vector(file, data.hardNodes);
  }

  // CONNECTM
  if (data.connectM.cols())
  {
    file << "CONNECTM"  << std::endl;
    write_matrix(file, data.connectM);
  }

  // NEIGHBORS
  if (data.neighbors.cols())
  {
    file << "NEIGHBORS"  << std::endl;
    write_matrix(file, data.neighbors);
  }

  // ANCESTORS
  if (!data.ancestors.empty())
  {
    file << "ANCESTORS"  << std::endl;
    write_vector(file, data.ancestors);
  }

  // COLORS
  if (!data.colors.empty())
  {
    file << "COLORS"  << std::endl;
    write_vector(file, data.colors);
  }

  // NEIGHBORS
  if (data.metrics.cols())
  {
    file << "METRICS"  << std::endl;
    write_matrix(file, data.metrics);
  }

  // MASK
  if (!data.mask.empty())
  {
    file << "MASK"  << std::endl;
    write_vector(file, data.mask);
  }
    
  file << "End" << std::endl;
  file.close();
 
  return 0;
}

///
E_Int iodata::read(const char* filename, DELAUNAY::MeshData& data)
{
  std::string                  line, entity;
  std::vector<std::string>     words;
  std::ifstream                file (filename);
  std::set<std::string>        keys;
  E_Int                        rows{ 0 }, cols{ 0 };
  E_Float                      P[3];
  E_Int                        K[3], i;

  if (!file.is_open())
    return 1;

  keys.insert("POS");
  keys.insert("CONNECTB");
  keys.insert("HARDNODES");
  keys.insert("CONNECTM");
  keys.insert("NEIGHBORS");
  keys.insert("ANCESTORS");
  keys.insert("COLORS");
  keys.insert("METRICS");
  keys.insert("MASK");
  keys.insert("End");
  
  K_FLD::IntArray* cb = const_cast<K_FLD::IntArray*>(data.connectB);
  
  while (std::getline (file, line))
  {
    if (line.empty())
      continue;
    
    __get_words(line, ' ', words);

    if (keys.find(words[0]) != keys.end())
    {
      if (words[0] == "End")
        break;
	  
      entity = words[0];
      std::getline (file, line);
      __get_words(line, ' ', words);
      rows = atoi(words[0].c_str());
      cols = 0;
      if (words.size() > 1)
        cols = atoi(words[1].c_str());

      if (entity == "POS")
      {
        data.pos->clear();
        data.pos->reserve(rows, cols);
      }
      else if (entity == "CONNECTB")
      {
        cb->clear();
        cb->reserve(rows, cols);
      }
      else if (entity == "HARDNODES")
      {
        data.hardNodes.clear();
        data.hardNodes.reserve(rows);
      }
      else if (entity == "CONNECTM")
      {
        data.connectM.clear();
        data.connectM.reserve(rows, cols);
      }
      else if (entity == "NEIGHBORS")
      {
        data.neighbors.clear();
        data.neighbors.reserve(rows, cols);
      }
      else if (entity == "ANCESTORS")
      {
        data.ancestors.clear();
        data.ancestors.reserve(rows);
      }
      else if (entity == "COLORS")
      {
        data.colors.clear();
        data.colors.reserve(rows);
      }
      else if (entity == "METRICS")
      {
        data.metrics.clear();
        data.metrics.reserve(rows, cols);
      }
      else if (entity == "MASK")
      {
        data.mask.clear();
        data.mask.reserve(rows);
      }
      continue;
    }
    else if (entity == "POS")
    {
      for (i = 0; i < rows; ++i)
        P[i] = atof(words[i].c_str());
      data.pos->pushBack(P, P+rows);
    }
    else if (entity == "CONNECTB")
    {
      for (i = 0; i < 2; ++i)
        P[i] = atoi(words[i].c_str());
      
      cb->pushBack(P, P+2);
    }
    else if (entity == "HARDNODES")
    {
       data.hardNodes.push_back(atoi(words[0].c_str()));
    }
    else if (entity == "CONNECTM")
    {
      for (i = 0; i < rows; ++i)
        K[i] = atoi(words[i].c_str());
      data.connectM.pushBack(K, K+rows);
    }
    else if (entity == "NEIGHBORS")
    {
      for (i = 0; i < rows; ++i)
        K[i] = atoi(words[i].c_str());
      data.neighbors.pushBack(K, K+rows);
    }
    else if (entity == "ANCESTORS")
    {
      data.ancestors.push_back(atoi(words[0].c_str()));
    }
    else if (entity == "COLORS")
    {
       data.colors.push_back(atoi(words[0].c_str()));
    }
    else if (entity == "METRICS")
    {
      for (i = 0; i < rows; ++i)
        P[i] = atof(words[i].c_str());
      data.metrics.pushBack(P, P+rows);
    }
    else if (entity == "MASK")
    {
       data.mask.push_back(atoi(words[0].c_str()));
    }
  }

  return 0;
}




}
