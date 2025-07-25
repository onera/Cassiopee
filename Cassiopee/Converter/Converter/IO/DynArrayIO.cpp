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

#ifndef __DYNARRAYIO_CPP__
#define	__DYNARRAYIO_CPP__

#include "IO/DynArrayIO.h"
#include <vector>
#include "Fld/FldArray.h"
# include "Nuga/include/DynArray.h"
#include <string.h>
#include "IO/GenIO.h"
#include <iostream>
#include <sstream>
#include "Connect/connect.h"
#include "kcore.h"
#include "converter.h"
#include "Nuga/include/merge.h"

std::string K_CONVERTER::DynArrayIO::rdir = "";
std::string K_CONVERTER::DynArrayIO::wdir = "";

namespace K_CONVERTER
{

// Writes a point cloud into a file.
E_Int DynArrayIO::write(const char* filename, const K_FLD::FloatArray& coord)
{
  return write(filename, coord, K_FLD::IntArray());
}
  
E_Int DynArrayIO::read
(const char* fileName,
 std::vector<K_FLD::FloatArray>& coords, std::vector<K_FLD::IntArray>& connects)
{
  std::vector<K_FLD::FldArrayF*> field;// field read for each zone
  std::vector<E_Int> im,jm,km;
  char* varString = NULL;
  std::vector<K_FLD::FldArrayI*> c, cc;
  std::vector<K_FLD::FldArrayF*> ufield;
  std::vector<E_Int> et;
  std::vector<char*> zoneNames;
  char* varStringc = NULL; // added for center fields read
  std::vector<K_FLD::FldArrayF*> fieldc;
  std::vector<K_FLD::FldArrayF*> ufieldc;
  
  E_Int ret = 1;
  
  E_Int l = strlen(fileName);
  char* fname=0;
  
  if (rdir == "")
  {
    //const char* -> char*
    fname = new char [l+1];
    strcpy(fname, fileName); 
  }
  else
  {
    //const char* -> char*
    fname = new char [l+rdir.size()+1];
    strcpy(fname, rdir.c_str()); 
    strcat(fname, fileName); 
  }
  
  const char* fileFmt = get_fmt(fileName);
  if (!fileFmt) return 1;
  
  if (strcmp(fileFmt, "bin_tp") == 0)
  {
    // Binary tecplot read
    ret = K_IO::GenIO::getInstance()->tecread
            (fname, varString, field, im, jm, km, ufield, c, et, zoneNames, varStringc, fieldc, ufieldc);
  }
  else if (strcmp(fileFmt, "fmt_tp") == 0)
  {
    // Formatted tecplot read
    ret = K_IO::GenIO::getInstance()->tpread
            (fname, varString, field, im, jm, km, ufield, c, et, zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_mesh") == 0)
  {
    // Formatted mesh read
    ret = K_IO::GenIO::getInstance()->meshread(fname, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  
  delete [] fname;
  
  if (ret) return ret;
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  //E_Float eps = 1.e-12;
  
  size_t sz = zoneNames.size();
  connects.resize(sz);
  coords.resize(sz);
  cc.resize(sz);
   
  for (size_t i = 0; i < sz ; ++i)
  {
    bool structured = (field.size() > i && field[i]);
    K_FLD::FldArrayF* cur_field = (structured ? field[i] : ufield[i]);  
    
    E_Int shift = (et[i] == 8) ? 0 : -1;
        
    //connectivity
    K_FLD::IntArray& cont = connects[i];
    if (!structured)
      cont = K_FLD::IntArray(*c[i], shift);
    //else
      //build_connectivity(im[i], jm[i], km[i], cont);
    
    //coordinates
    K_FLD::FloatArray& coord = coords[i];
    
    //E_Int nb_pts = cur_field->getSize();
    //E_Int dim = cur_field->getNfld();
    coord = K_FLD::FloatArray(*cur_field);
  }
  
  // clean
  //field c ufield zoneNames
  for (size_t i = 0; i < field.size(); ++i)
    delete field[i];
  for (size_t i = 0; i < ufield.size(); ++i)
    delete ufield[i];
  for (size_t i = 0; i < c.size(); ++i)
    delete c[i];
  for (size_t i = 0; i < zoneNames.size(); ++i)
    delete zoneNames[i];
  for (size_t i = 0; i < cc.size(); ++i)
    delete cc[i]; 
  
  return ret;
  
}

E_Int DynArrayIO::read
(const char* fileName,
 std::vector<K_FLD::FldArrayF>& coords, std::vector<K_FLD::FldArrayI>& connects)
{
  std::vector<K_FLD::FldArrayF*> field;            // field read for each zone
  std::vector<E_Int> im,jm,km;
  char* varString = NULL;
  std::vector<K_FLD::FldArrayI*> c, cc;
  std::vector<K_FLD::FldArrayF*> ufield;
  std::vector<E_Int> et;
  std::vector<char*> zoneNames;
  char* varStringc = NULL; // centers
  std::vector<K_FLD::FldArrayF*> fieldc; 
  std::vector<K_FLD::FldArrayF*> ufieldc; 
  E_Int ret = 1;
  
  E_Int l = strlen(fileName);
  char* fname=0;
  
  if (rdir == "")
  {
    //const char* -> char*
    fname = new char [l+1];
    strcpy(fname, fileName); 
  }
  else
  {
    //const char* -> char*
    fname = new char [l+rdir.size()+1];
    strcpy(fname, rdir.c_str()); 
    strcat(fname, fileName);
  }
  
  const char* fileFmt = get_fmt(fileName);
  if (!fileFmt) return 1;
  
  if (strcmp(fileFmt, "bin_tp") == 0)
  {
    // Binary tecplot read
    ret = K_IO::GenIO::getInstance()->tecread
            (fname, varString, field, im, jm, km, ufield, c, et, zoneNames, varStringc, fieldc, ufieldc);
  }
  else if (strcmp(fileFmt, "fmt_tp") == 0)
  {
    // Formatted tecplot read
    ret = K_IO::GenIO::getInstance()->tpread
            (fname, varString, field, im, jm, km, ufield, c, et, zoneNames);
  }
   else if (strcmp(fileFmt, "fmt_mesh") == 0)
  {
    // Formatted mesh read
    ret = K_IO::GenIO::getInstance()->meshread(fname, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  
  delete [] fname;
  
  if (ret) return ret;
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  //E_Float eps = 1.e-12;
  
  size_t sz = zoneNames.size();
  connects.resize(sz);
  coords.resize(sz);
  cc.resize(sz);
   
  for (size_t i = 0; i < sz ; ++i)
  {
    bool structured = (field.size() > i && field[i]);
    K_FLD::FldArrayF* cur_field = (structured ? field[i] : ufield[i]);
    
    //connectivity
    K_FLD::FldArrayI& cont = connects[i];
    if (!structured) cont = *c[i];
    
    //coordinates
    K_FLD::FldArrayF& coord = coords[i];
    
    //E_Int nb_pts = cur_field->getSize();
    //E_Int dim = cur_field->getNfld();
    coord = *cur_field;
  }
  
  // clean
  //field c ufield zoneNames
  for (size_t i = 0; i < field.size(); ++i)
    delete field[i];
  for (size_t i = 0; i < ufield.size(); ++i)
    delete ufield[i];
  for (size_t i = 0; i < c.size(); ++i)
    delete c[i];
  for (size_t i = 0; i < zoneNames.size(); ++i)
    delete zoneNames[i];
  for (size_t i = 0; i < cc.size(); ++i)
    delete cc[i]; 
  
  return ret;
  
}

E_Int DynArrayIO::write
(const char* fileName,
 const std::vector<K_FLD::FloatArray>& coords,
const std::vector<K_FLD::IntArray>& connects,
 const std::vector<std::string>& elt_type)
{
  if (coords.empty() || connects.empty()) return 0;
  
  const char* fileFmt = get_fmt(fileName);
  if (!fileFmt) return 1;
  
  std::vector<char*> zoneNames;
  std::vector<K_FLD::FldArrayF*> dummyfield;
  std::vector<E_Int> im_dum, jm_dum, km_dum;
  std::vector<K_FLD::FloatArray> coordz(coords);
  char datfmt[20]; strcpy(datfmt, "%.9e ");
  char varString[128]; strcpy(varString, "x,y,z");
  if (coordz[0].rows() == 2)
  {
    E_Float zero = 0.;
    coordz[0].resize(3, coordz[0].cols(), &zero);
  }
  std::vector<K_FLD::FldArrayI*> c;
  std::vector<K_FLD::FldArrayF*> ufield;
  std::vector<E_Int> et;
  E_Int ret = 1;
  size_t sz = connects.size();
  
  ufield.resize(sz);
  c.resize(sz);
  
  std::ostringstream o;
  for (size_t i = 0; i < sz; ++i)
  {
    o.str("");
    o << "zone_" << i;
    zoneNames.push_back(const_cast<char*>(o.str().c_str()));
    
    ufield[i] = new K_FLD::FldArrayF;
    coordz[i].convert(*ufield[i], 0);   
    
    c[i] = new K_FLD::FldArrayI;
    connects[i].convert(*c[i], 1);
    
    et.push_back(getElementTypeId(elt_type[i].c_str()));
  }
  
  std::vector< std::vector<E_Int> > ets;
  for (size_t i = 0; i < sz; ++i)
  {
    std::vector<E_Int> v(1); v[0] = et[i];
    ets.push_back(v);
  }

  E_Int l = strlen(fileName);
  char* fname=0;
  
  if (wdir == "")
  {
    //const char* -> char*
    fname = new char [l+1];
    strcpy(fname, fileName);
  }
  else
  {
    //const char* -> char*
    fname = new char [l+wdir.size()+1];
    strcpy(fname, wdir.c_str()); 
    strcat(fname, fileName); 
  }
    
  if (strcmp(fileFmt, "bin_tp") == 0) // binary tecplot
  {
    ret = 
      K_IO::GenIO::getInstance()->tecwrite(fname, datfmt, varString,
                                           im_dum, jm_dum, km_dum,
                                           dummyfield, ufield, c, ets,
                                           zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_tp") == 0) // fmt tecplot
  { 
    ret = K_IO::GenIO::getInstance()->tpwrite(fname, datfmt, varString,
                                              im_dum, jm_dum, km_dum,
                                              dummyfield, ufield, c, ets,
                                              zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_mesh") == 0) // fmt mesh
  {
    if (dummyfield.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    ret = K_IO::GenIO::getInstance()->meshwrite(fname, datfmt, varString,
                                                im_dum, jm_dum, km_dum,
                                                dummyfield, ufield, c, ets,
                                                zoneNames);
  }
  
  // clean
  delete [] fname;
  //field c ufield zoneNames
  for (size_t i = 0; i < ufield.size(); ++i) delete ufield[i];
  for (size_t i = 0; i < c.size(); ++i) delete c[i];
  
  return ret;
}

E_Int DynArrayIO::write
(const char* fileName, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect,
 const char* elt_type, const std::vector<bool> * keep, const std::vector<E_Int>* colors)
{
  if (coord.cols() == 0) return 0;
  
  const K_FLD::FloatArray *pcrd(&coord);
  K_FLD::FloatArray crd;
  if (coord.rows() == 2){
    crd=coord;
    crd.resize(3, crd.cols(), 0.);
    pcrd=&crd;
  }
    
  
  const char* fileFmt = get_fmt(fileName);
  if (!fileFmt) return 1;
  
  std::vector<char*> zoneNames;
  std::vector<K_FLD::FldArrayF*> dummyfield;
  std::vector<E_Int> im_dum, jm_dum, km_dum;
  char datfmt[20]; strcpy(datfmt, "%.9e ");
  char varString[128]; strcpy(varString, "x,y,z");
  if (pcrd->rows() == 2) strcpy(varString, "x,y");
  std::vector<K_FLD::FldArrayI*> c;
  std::vector<K_FLD::FldArrayF*> ufield;
  std::vector<E_Int> et;
  E_Int ret = 1;
    
  ufield.resize(1);
  c.resize(1);
  char* zoneName = new char [128];
  strcpy(zoneName, "zone1");
  zoneNames.push_back(zoneName);
  K_FLD::IntArray connect1;
  std::vector<E_Int> new_cols;
  if (keep)
  {
    E_Int ROWS = connect.rows();
    for (E_Int i = 0; i < connect.cols(); ++i)
    {
      if ((*keep)[i])
      {
        connect1.pushBack(connect.col(i), connect.col(i)+ROWS);
        if (colors)new_cols.push_back((*colors)[i]);
      }
    }
  }
  else
  {
    connect1 = connect;
    if (colors) new_cols = *colors;
  }
    
  ufield[0] = new K_FLD::FldArrayF;
  pcrd->convert(*ufield[0], 0);
  
  E_Int id;
  if (!elt_type)
    id = connect1.rows()-1;
  else
    id = getElementTypeId(elt_type);
  
  et.push_back(id);
  std::vector< std::vector<E_Int> > ets(1); ets[0] = et;

  c[0] = new K_FLD::FldArrayI;
  E_Int shift = (et[0] == 8) ? 0 : 1;
  connect1.convert(*c[0], shift);
      
  E_Int l = strlen(fileName);
  char* fname=0;
  
  if (wdir == "")
  {
    //const char* -> char*
    fname = new char [l+1];
    strcpy(fname, fileName);
  }
  else
  {
    //const char* -> char*
    fname = new char [l+wdir.size()+1];
    strcpy(fname, wdir.c_str()); 
    strcat(fname, fileName); 
  }
  
  if (strcmp(fileFmt, "bin_tp") == 0) // binary tecplot
  {
    ret = K_IO::GenIO::getInstance()->tecwrite(fname, datfmt, varString,
                                               im_dum, jm_dum, km_dum,
                                               dummyfield, ufield, c, ets,
                                               zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_tp") == 0) // fmt tecplot
  { 
    ret = K_IO::GenIO::getInstance()->tpwrite(fname, datfmt, varString,
                                              im_dum, jm_dum, km_dum,
                                              dummyfield, ufield, c, ets,
                                              zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_mesh") == 0) // fmt mesh
  {
    if (dummyfield.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    std::vector<std::vector<E_Int> > cols;
    if (colors)
      cols.push_back(new_cols);
    
    ret = K_IO::GenIO::getInstance()->meshwrite(fname, datfmt, varString,
                                                im_dum, jm_dum, km_dum,
                                                dummyfield, ufield, c, ets,
                                                zoneNames, (colors)? &cols : 0);
  } 
  
  // clean
  delete [] fname;
  
  //field c ufield zoneNames
  for (size_t i = 0; i < ufield.size(); ++i)
    delete ufield[i];
  for (size_t i = 0; i < c.size(); ++i)
    delete c[i];
  
  return ret;
}

E_Int DynArrayIO::write
  (const char* fileName, const K_FLD::FldArrayF& coord, const K_FLD::FldArrayI& connect,
   const char* elt_type, const std::vector<bool>* mask, const std::vector<E_Int>* colors)
{
  if (coord.getSize()*connect.getSize() == 0) return 0;
  
  const char* fileFmt = get_fmt(fileName);
  if (!fileFmt) return 1;
  
  std::vector<char*> zoneNames;
  std::vector<K_FLD::FldArrayF*> dummyfield;
  std::vector<E_Int> im_dum, jm_dum, km_dum;
  char datfmt[20]; strcpy(datfmt, "%.9e ");
  char varString[128]; strcpy(varString, "x,y,z");
  if (coord.getNfld() == 2) strcpy(varString, "x,y");
  std::vector<K_FLD::FldArrayI*> c;
  std::vector<K_FLD::FldArrayF*> ufield;
  std::vector<E_Int> et;
  E_Int ret = 1;
    
  ufield.resize(1);
  c.resize(1);
  
  char* zoneName = new char [128];
  strcpy(zoneName, "zone1");
  zoneNames.push_back(zoneName);
  K_FLD::FldArrayI connect1;
  if (mask)
  {
    E_Int ROWS = connect.getNfld();
    E_Int* entry = new E_Int[ROWS];
    for (E_Int i = 0; i < connect.getSize(); ++i)
    {
      if ((*mask)[i])
      {
        for (E_Int j = 0; j < ROWS; ++j)
          entry[j] = connect(i,j+1);
        connect1.pushBack(entry, entry+ROWS);
      }
    }
    delete [] entry;
  }
  else
    connect1=connect;
    
  ufield[0] = const_cast<K_FLD::FldArrayF*>(&coord);
    
  et.push_back(getElementTypeId(elt_type));
  std::vector< std::vector<E_Int> > ets(1); ets[0] = et;

  c[0] = &connect1;

  int l = strlen(fileName);
  char* fname=0;
  
  if (wdir == "")
  {
    //const char* -> char*
    fname = new char [l+1];
    strcpy(fname, fileName); 
  }
  else
  {
    //const char* -> char*
    fname = new char [l+wdir.size()+1];
    strcpy(fname, wdir.c_str()); 
    strcat(fname, fileName); 
  }
  
  if (strcmp(fileFmt, "bin_tp") == 0) // binary tecplot
  {
    ret = K_IO::GenIO::getInstance()->tecwrite(fname, datfmt, varString,
                                               im_dum, jm_dum, km_dum,
                                               dummyfield, ufield, c, ets,
                                               zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_tp") == 0) // fmt tecplot
  { 
    ret = K_IO::GenIO::getInstance()->tpwrite(fname, datfmt, varString,
                                              im_dum, jm_dum, km_dum,
                                              dummyfield, ufield, c, ets,
                                              zoneNames);
  }
  else if (strcmp(fileFmt, "fmt_mesh") == 0) // fmt mesh
  {
    if (dummyfield.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    std::vector<std::vector<E_Int> > cols;
    if (colors)
      cols.push_back(*colors);
    
    ret = K_IO::GenIO::getInstance()->meshwrite(fname, datfmt, varString,
                                                im_dum, jm_dum, km_dum,
                                                dummyfield, ufield, c, ets,
                                                zoneNames, (colors)? &cols : 0);
  } 
  
  // clean
  delete [] fname;
    
  return ret;
}

E_Int DynArrayIO::getElementTypeId(const char* eltType)
{
  if (strcmp(eltType, "NODE") == 0) // NODE-> structure
    return 0;
  if (strcmp(eltType, "BAR") == 0)
    return 1;
  if (strcmp(eltType, "TRI") == 0)
    return 2;
  if (strcmp(eltType, "QUAD") == 0)
    return 3;
  if (strcmp(eltType, "TETRA") == 0)
    return 4;
  if (strcmp(eltType, "PYRA") == 0)
    return 5;
  if (strcmp(eltType, "PENTA") == 0)
    return 6;
  if (strcmp(eltType, "HEXA") == 0)
    return 7;
  if (strcmp(eltType, "NGON") == 0)
    return 8;
  return -1;//unknown
}

const char* DynArrayIO::getElementType(E_Int id)
{
  if (id == 0)
    return "NODE";
  if (id == 1)
    return "BAR";
  if (id == 2)
    return "TRI";
  if (id == 3)
    return "QUAD";
  if (id == 4)
    return "TETRA";
  if (id == 5)
    return "PYRA";
  if (id == 6)
    return "PENTA";
  if (id == 7)
    return "HEXA";
  if (id == 8)
    return "NGON";
  return ""; 
}

const char* DynArrayIO::get_fmt(const char* fname)
{
  std::string tmp(fname);
  if (tmp.find(".plt") != std::string::npos)
    return "bin_tp";
  if (tmp.find(".tp") != std::string::npos)
    return "fmt_tp";
  if (tmp.find(".mesh") != std::string::npos)
    return "fmt_mesh";
  std::cout << fname << " is not a supported format" << std::endl;
  return 0;
}
/*
void DynArrayIO::build_connectivity(E_Int im, E_Int jm, E_Int km, K_FLD::IntArray& connect)
{
  E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  for (E_Int k = 0; k < km; k++)
  {
    for (E_Int j = 0; j < jm; j++)
    {
      for (E_Int i = 0; i < im; i++)
      {
        ind1 = 1 + i + j*nil + k*ninj; //A(  i,  j,k)
        ind2 = ind1 + 1;              //B(i+1,  j,k)
        ind3 = ind2 + nil;            //C(i+1,j+1,k)
        ind4 = ind3 - 1;              //D(  i,j+1,k)
        ind5 = ind1 + ninj;           //E(  i,  j,k+1)
        ind6 = ind2 + ninj;           //F(i+1,  j,k+1)
        ind7 = ind3 + ninj;           //G(i+1,j+1,k+1)
        ind8 = ind4 + ninj;           //H(  i,j+1,k+1) 
        ind = i+j*ni1+k*ni1*nj1;
        cn1[ind] = ind1;
        cn2[ind] = ind2;
        cn3[ind] = ind3;
        cn4[ind] = ind4;
        cn5[ind] = ind5;
        cn6[ind] = ind6;
        cn7[ind] = ind7;
        cn8[ind] = ind8;
      }
    }
  }
}
*/
}

#endif	/* __DYNARRAYIO_CPP__ */

