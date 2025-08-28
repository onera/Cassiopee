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
#ifndef _CONVERTER_GENIO_H_
#define _CONVERTER_GENIO_H_

#include <cctype>
#include <locale>
# include <vector>
# include <list>
# include <stdio.h>
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefCplusPlusConst.h"
# include "String/kstring.h"
# include "kPython.h"
# include "Nuga/include/ngon_t.hxx"

// Define ftell and fseek (long return)
#ifdef _WIN64
#define KFTELL _ftelli64
#define KFSEEK _fseeki64
#else
#define KFTELL ftell
#define KFSEEK fseek
#endif

// Define conversion from int big endian to int little endian
#define SBE(x) *((short*)K_IO::GenIO::conv2((char*)&x))
// Define conversion from int big endian to int little endian
#define IBE(x) *((int*)K_IO::GenIO::conv4((char*)&x))
// Define conversion from unsigned int big endian to unsigned int little endian
#define UIBE(x) *((unsigned int*)K_IO::GenIO::conv4((char*)&x))
// Define conversion from long big endian to long little endian
#define LBE(x) *((E_LONG*)K_IO::GenIO::conv8((char*)&x))
// Define conversion between float big endian to double little endian
#define FBE(x) *((float*)K_IO::GenIO::conv4((char*)&x))
// Define conversion between double big endian to double little endian
#define DBE(x) *((double*)K_IO::GenIO::conv8((char*)&x))

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

namespace K_IO
{
  const E_Int BUFSIZE = 8192;
// ============================================================================
// @Name GenIO
// @Memo Input/output routines
// @See
/* @Text
*/
//=============================================================================
class GenIO
{
  public:

    ///+ 1- Access to the singleton.
    static GenIO* getInstance();
    ///-

    ///+ 2- Misc Functions
    /** Wait s milliseconds. */
    E_Int wait_micsec(E_Int s);
    /** Return the convert endian variable. */
    E_Bool getConvertEndian();
    /** Return endianess of the current machine:
        1: big endian, 0: little endian. */
    E_Int machineEndianess();
    /** Convert endians on a field. */
    void convertEndianField(FldArrayF& field);
    /** Convert endians on a field. */
    void convertEndianField(FldArrayI& field);
    /** Convert two bytes Little<->Big endian. */
    static char* conv2(char* x);
    /** Convert four bytes Little<->Big endian. */
    static char* conv4(char* x);
    /** Convert eight bytes Little<->Big endian. */
    static char* conv8(char* x);
    /** Retourne 1 si word match dans buf */
    E_Int matchInString(char* buf, const char* word);
    ///-

    ///+ 3- Binary tecplot file functions
    /** All routines return 1 if file doesnt exist,
        return 2 if zone doesnt exist,
        and return 0 if ok. */
    /** Return nfield.
        Fill the varString string with the variables name. */
    E_Int tecstat(char* file, char* varString,
                  E_Int& nfield);
    /** Read all zones in file.
        varString must be allocated by the calling routine. */
    E_Int tecread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames,
      char*& centerVarString,
      std::vector<FldArrayF*>& centerStructField,
      std::vector<FldArrayF*>& centerUnstructField);
    /** Write structured field in binary tec format. One zone.
        return 1 if failed */
    E_Int tecwrite(char* file, char* dataFmt, char* varstring,
                   E_Int ni, E_Int nj, E_Int nk,
                   const FldArrayF& coord, const FldArrayF& field,
                   std::vector<char*>& zoneNames);
    /** Write structured field in binary tec format. Multiple zones. */
    E_Int tecwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<char*>& zoneNames);
    /** Write unstructured field in tec format. Multiple zones.
        element type: 2 (TRI), 3 (QUAD), 4 (TETRA), 7 (HEXA).
        return 1 if failed */
    E_Int tecwrite(char* file, char* dataFmt, char* varString,
                   std::vector<FldArrayF*>& unstructField,
                   std::vector<FldArrayI*>& connect,
                   std::vector< std::vector<E_Int> >& eltTypes,
                   std::vector<char*>& zoneNames);
    /** Write structured and unstructured field in tec format.
        Multiple zones. */
     E_Int tecwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    E_Int tecwrite75(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    E_Int tecwrite108(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);

    /** Write triangles (tp format) described by field and connect.
	Field contains coord and other possible fields. */
    E_Int tpwriteTriangles(char* file, FldArrayF& field, FldArrayI& connect,
                           E_Bool add=false);
    /** Write quadrangles (tp format) described by field and connect.
	Field contains coord and other possible fields. */
    E_Int tpwriteQuads(char* file, FldArrayF& field, FldArrayI& connect,
                       E_Bool add=false);
    /** Write text "number" near each point */
    E_Int tpwriteText(char* file, FldArrayF& field, FldArrayI& number);

    ///+ bin v3d file functions
    /** Read */
    E_Int v3dread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    /** Write */
    E_Int v3dwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field,
      std::vector<char*>& zoneNames, E_Int isize=8,
      E_Int rsize=8, E_Bool convertEndian=false);
    ///-

    ///+ bin plot3d file functions
    /** Read */
    E_Int plot3dread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    /** Write */
    E_Int plot3dwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field,
      std::vector<char*>& zoneNames, E_Int isize=8,
      E_Int rsize=8, E_Bool convertEndian=false);
    ///-

    ///+ fmt plot3d file functions
    /** Read */
    E_Int fp3dread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    /** Write */
    E_Int fp3dwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    ///-

    ///+ fmt tecplot file functions
    /** Read */
    E_Int tpread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int tpread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<std::vector<E_Int> >& eltType, std::vector<char*>& zoneNames,
      E_Int api=1);
    E_Int tpwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames);
    E_Int tpwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ fmt v3d file functions
    /** Read */
    E_Int fv3dread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    /** Write */
    E_Int fv3dwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, std::vector<char*>& zoneNames);
    ///-

    ///+ fmt su2 file functions
    /** Read */
    E_Int su2read(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames,
      std::vector<FldArrayI*>& BCFaces, std::vector<char*>& BCNames);
    E_Int su2read(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<std::vector<E_Int> >& eltType, std::vector<char*>& zoneNames,
      std::vector<FldArrayI*>& BCFaces, std::vector<char*>& BCNames,
      E_Int api=1);
    /** Write */
    E_Int su2write(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames,
      PyObject* BCFaces);
    E_Int su2write(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames,
      PyObject* BCFaces);
    ///-

    ///+ fmt foam file functions
    /** Read */
    E_Int foamread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames,
      std::vector<FldArrayI*>& BCFaces, std::vector<char*>& BCNames,
      std::vector<FldArrayF*>& BCFields,
      char*& varStringc,
      std::vector<FldArrayF*>& centerStructField,
      std::vector<FldArrayF*>& centerUnstructField);
    /** Write */
    E_Int foamwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames,
      PyObject* BCFaces);
    E_Int foamWritePoints(char* file, FldArrayF& f);
    E_Int foamReadPoints(char* file, FldArrayF& f);
    E_Int foamReadFields(char* file,
      std::vector<FldArrayF*>& centerUnstructField, E_Int ncells,
      char*& varStringc,
      const std::vector<char*> &BCNames,
      const std::vector<FldArrayI*> &BCFaces,
      std::vector<FldArrayF*> &BCFields, E_Int *owner,
      const std::vector<E_Float> &delta, E_Int nifaces,
      const std::vector<E_Int> &indir);
    E_Int foamWriteFaces(char* file, FldArrayI& cn,
      const std::vector<E_Int>& faces);
    E_Int foamReadFaces(char* file, E_Int& nfaces, FldArrayI& cn);
    E_Int foamWriteOwner(char* file, const std::vector<E_Int> &owner,
      const std::vector<E_Int> &faces);
    E_Int foamReadOwner(char* file, FldArrayI& PE);
    E_Int foamWriteNeighbour(char* file, const std::vector<E_Int>& neigh,
      const std::vector<E_Int>& faces);
    E_Int foamReadNeighbour(char* file, FldArrayI& PE);
    E_Int foamWriteBoundary(char* file, const std::vector<char*>& bc_names, 
      const std::vector<E_Int>& bc_nfaces,
      const std::vector<E_Int>& bc_startfaces);
    E_Int foamReadBoundary(char* file, std::vector<FldArrayI*>& BCFaces,
      std::vector<char*>& BCNames, std::vector<E_Int> &indir);
    E_Int readScalarField(char *file, FldArrayF& f, E_Int idx, E_Int* owner,
      const std::vector<FldArrayI*>& BCFaces,
      std::vector<FldArrayF*> &BCFields, const std::vector<E_Float>& delta,
      E_Int nifaces, const std::vector<E_Int>& indir);
    E_Int readVectorField(char* file, FldArrayF& f, E_Int idx, E_Int* owner,
      const std::vector<FldArrayI*>& BCFaces,
      std::vector<FldArrayF*>& BCFields, const std::vector<E_Float>& delta,
      E_Int nifaces, const std::vector<E_Int>& indir);
    E_Int readTensorField(char* file, FldArrayF& f, E_Int idx);
    ///-

    ///+ Povray functions
    /** Mesh2 Read */
    E_Int povread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** Mesh2 write */
    E_Int povwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames, E_Int colormap);
    /** Density df3 write */
    E_Int df3write(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    /* Create colormaps */
    void createColormap(E_Int colormap, FldArrayF& rgb);
    ///-

    ///+ Mesh (INRIA) functions
    /** Mesh read */
    E_Int meshread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int meshread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<std::vector<E_Int> >& eltType, std::vector<char*>& zoneNames,
      E_Int api=1);
    /** Mesh write */
    E_Int meshwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames,
      std::vector<std::vector<E_Int> > * colors=0);
    E_Int meshwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames,
      std::vector<std::vector<E_Int> > * colors=0);
    ///-

    ///+ Gmsh (Louvain) functions
    /** Formatted Gmsh read */
    E_Int gmshread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, 
      std::vector<char*>& zoneNames);
    E_Int gmshread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<std::vector<E_Int> >& eltType, std::vector<char*>& zoneNames,
      E_Int api=1);
    /** Formatted Gmsh write */
    E_Int gmshwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames);
    E_Int gmshwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    /** Binary Gmsh read */
    E_Int bingmshread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** Binary Gmsh write */
    E_Int bingmshwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ Obj functions
    /** Obj read */
    E_Int objread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** Obj write */
    E_Int objwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ STL functions
    /** bin STL */
    E_Int stlread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int stlwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    /** fmt STL */
    E_Int fstlread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int fstlwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ fmt selig functions
    E_Int seligread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int seligwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);

    /** PLY (Stanford) read */
    E_Int plyread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int plywrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    /** 3DS (3D studio) read */
    E_Int f3dsread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    E_Int f3dswrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ GLTF functions
    /** gltf read */
    E_Int gltfread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** gltf write */
    E_Int gltfwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ VTK legacy functions
    /** vtk read */
    E_Int binvtkread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** gltf write */
    E_Int binvtkwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ WAV functions
    /** wavwrite */
    E_Int wavwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ PNG functions
    /** pngread */
    E_Int pngread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** pngwrite */
    E_Int pngwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ JPG functions
    /** jpgread */
    E_Int jpgread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** jpgwrite */
    E_Int jpgwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ Xfig functions
    /** Xfig read */
    E_Int xfigread(
      char* file, char*& varString, E_Int NptsCurve, E_Int NptsLine,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    /** XFig write */
    E_Int xfigwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ SVG functions
    /** svg read */
    E_Int svgread(
      char* file, char*& varString, E_Float density, E_Int NptsCurve,
      E_Int NptsLine,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames);
    
    /** svg write */
    E_Int svgwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    /* Lecture des coordonnees x,y, ds le fichier
       retourne 0 si fin du fichier, 1 sinon */
    E_Int readTwoCoordinates(FILE* ptrFile, E_Float* pt);
    E_Int readTwoCoordinates(char* buf, FILE* ptrFile, E_Float* pt);
    ///-

    ///+ GTS functions
    /** gts read */
    E_Int gtsread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, 
      std::vector<char*>& zoneNames);
    /** gts write */
    E_Int gtswrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ CEDRE functions
    /** cedre read */
    E_Int cedreread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames,
      std::vector<FldArrayI*>& BCFaces,
      std::vector<char*>& BCNames);

    /** cedre write */
    E_Int cedrewrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames,
      PyObject* BCFaces);
    ///-

    ///+ CEDRE archive functions
    /** archive cedre read */
    E_Int arcread(
      char* file, char*& varString,
      std::vector<FldArrayF*>& structField,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connectivity,
      std::vector<E_Int>& eltType, std::vector<char*>& zoneNames,
      char*& centerVarString,
      std::vector<FldArrayF*>& centerStructField,
      std::vector<FldArrayF*>& centerUntructField);

    /** archive cedre write */
    E_Int arcwrite(
      char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< std::vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames);
    ///-

    ///+ ADF/CGNS functions
    /* Lecture dans un arbre */
    E_Int adfcgnsread(char* file, PyObject*& tree, int skeleton=0,
                      int maxFloatSize=5, int maxDepth=-1);
    /* Ecriture d'un arbre */
    E_Int adfcgnswrite(char* file, PyObject* tree);
    /* Lecture a partir de chemins donnes (base/zone) */
    PyObject* adfcgnsReadFromPaths(char* file, PyObject* paths, 
                                   E_Int maxFloatSize=1.e6, E_Int maxDepth=-1);
    /* Ecrit des parties d'un arbre python */
    E_Int adfcgnsWritePaths(char* file, PyObject* nodeList, PyObject* paths,
      E_Int maxDepth=-1, E_Int mode=0);
    void getABFromPath(char* path, std::vector<char*>& pelts);
    /* Efface un chemin d'un fichier */
    E_Int adfcgnsDeletePaths(char* file, PyObject* paths);
    ///-

    ///+ HDF/CGNS functions
    /* Lecture dans un arbre */
    E_Int hdfcgnsread(char* file, PyObject*& tree, PyObject* dataShape,
                      PyObject* links, 
                      int skeleton=0, int maxFloatSize=5, int maxDepth=-1,
                      int readIntMode=0, PyObject* skipTypes=NULL);
    /* Ecriture d'un arbre */
    E_Int hdfcgnswrite(char* file, PyObject* tree, PyObject* links=NULL,
                       int writeIntMode=0, int writeRealMode=0);

    /* Lecture a partir de chemins donnes */
    PyObject* hdfcgnsReadFromPaths(char* file, PyObject* paths,
                                   E_Int maxFloatSize=1.e6, E_Int maxDepth=-1,
                                   E_Int readIntMode=0, PyObject* dataShape=NULL,
                                   PyObject* skipTypes=NULL, 
                                   PyObject* mpi4pyCom=NULL);
    PyObject* hdfcgnsReadFromPathsPartial(char* file, E_Int readIntMode,
                                          PyObject* Filters,
                                          PyObject* mpi4pyCom=NULL);
    /* Ecrit des parties d'un arbre python */
    E_Int hdfcgnsWritePaths(char* file, PyObject* nodeList, PyObject* paths,
                            PyObject* links=NULL, E_Int maxDepth=-1,
                            E_Int mode=0, E_Int isize=8, E_Int rsize=8);
    E_Int hdfcgnsWritePathsPartial(char* file, PyObject* tree,
                                   PyObject* Filter,
                                   int skeleton=0,
                                   PyObject* mpi4pyCom=NULL);
    E_Int hdfcgnsDeletePaths(char* file, PyObject* paths);
    void ripEndOfPath(char* path, char*& startPath);
    void getEndOfPath(char* path, char*& startPath);
    ///-

    ///+ NETCDF3/TAU functions
    /* Lecture dans un arbre */
    E_Int tauread(char* file, PyObject*& tree);
    /* Ecriture d'un arbre */
    E_Int tauwrite(char* file, PyObject* tree);

    ///+ HDF/FSDM functions
    /* Lecture dans un arbre */
    E_Int hdffsdmread(char* file, PyObject*& tree);
    /* Ecriture d'un arbre */
    E_Int hdffsdmwrite(char* file, PyObject* tree);

    ///+ CPlot functions
    /* Create the socket for communications */
    E_Int cplotClient(char* machine, char* service, int *mySocket);
    /* Send an integer to cplot server */
    void cplotSendInteger(int* mySocket, int value);
    /* Send a double to cplot server */
    void cplotSendDouble(int* mySocket, double value);
    /* Send the number of variables and grid dimensions to cplot server */
    void cplotSendHeader(int* mySocket, int nvar, int nfield,
			 int zonenumber, int ni, int nj, int nk);
    /** Send zone data to cplot server */
    void cplotSend(int* mySocket, int size, double* array);
    ///-

  private:
    /* Return 0: same endian for file and machine,
       1: different endians,
      -1: file doesn't exist. */
    E_Int tecCheckEndian(char* file);
    E_Int v3dCheckEndian(char* file);

    /* bin_tp: read the header,
       return nfield and fill varString and version. */
    E_Int readHeader(FILE *ptr_file, E_Int& nfield,
                     char*& varString, char* version);
    E_Int readHeaderCE(FILE *ptr_file, E_Int& nfield,
                       char*& varString, char* version);
    /* bin_tp: read the zone header. */
    E_Int readZoneHeader75(FILE* ptr_file, E_Int& dim,
                           E_Int& ni, E_Int& nj, E_Int& nk,
                           E_Int& npts, E_Int& nelts, E_Int& eltType,
                           char* zoneName, E_Int& dataPacking,
                           std::vector<FldArrayF*>& geom);
    E_Int readZoneHeader75CE(FILE* ptr_file, E_Int& dim,
                             E_Int& ni, E_Int& nj, E_Int& nk,
                             E_Int& npts, E_Int& nelts, E_Int& eltType,
                             char* zoneName, E_Int& dataPacking,
                             std::vector<FldArrayF*>& geom);
    E_Int readZoneHeader108(
      E_Int version, E_Int nfield,
      FILE* ptr_file, E_Int& dim,
      E_Int& ni, E_Int& nj, E_Int& nk,
      E_Int& npts, E_Int& nelts,
      E_Int& numFaces, E_Int& numFaceNodes,
      E_Int& numBoundaryFaces, E_Int& numBoundaryConnections,
      E_Int& eltType, E_Int& rawlocal,
      char* zoneName, E_Int& dataPacking,
      E_Int& strand, E_Float& time,
      std::vector<E_Int>& loc,
      std::vector<FldArrayF*>& geom);
    E_Int readZoneHeader108CE(
      E_Int version, E_Int nfield,
      FILE* ptr_file, E_Int& dim,
      E_Int& ni, E_Int& nj, E_Int& nk,
      E_Int& npts, E_Int& nelts,
      E_Int& numFaces, E_Int& numFaceNodes,
      E_Int& numBoundaryFaces, E_Int& numBoundaryConnections,
      E_Int& eltType, E_Int& rawlocal,
      char* zoneName, E_Int& dataPacking,
      E_Int& strand, E_Float& time, 
      std::vector<E_Int>& loc,
      std::vector<FldArrayF*>& geom);
    /* bin_tp: read data */
    E_Int readData108(E_Int version, FILE* ptrFile, E_Int ni, E_Int nj,
                      E_Int nk,
                      E_Int dataPacking, std::vector<E_Int>& loc,
                      FldArrayF* f, FldArrayF* fc);
    E_Int readData108(E_Int version, FILE* ptrFile,
                      E_Int dataPacking, std::vector<E_Int>& loc, E_Int et,
                      E_Int numFaces, E_Int numFaceNodes,
                      E_Int numBoundaryFaces, E_Int numBoundaryConnections,
                      E_Int ne, E_Int rawlocal,
                      FldArrayF* f, FldArrayI& c, FldArrayF* fc);
    E_Int readData108CE(E_Int version, FILE* ptrFile, E_Int ni, E_Int nj,
                        E_Int nk,
                        E_Int dataPacking, std::vector<E_Int>& loc,
                        FldArrayF* f, FldArrayF* fc);
    E_Int readData108CE(E_Int version, FILE* ptrFile,
                        E_Int dataPacking, std::vector<E_Int>& loc, E_Int et,
                        E_Int numFaces, E_Int numFaceNodes,
                        E_Int numBoundaryFaces, E_Int numBoundaryConnections,
                        E_Int ne, E_Int rawlocal,
                        FldArrayF* f, FldArrayI& c, FldArrayF* fc);
    E_Int readData75(FILE* ptrFile, E_Int ni, E_Int nj, E_Int nk,
                     E_Int dataPacking,
                     FldArrayF& f);
    E_Int readData75(FILE* ptrFile,
                     E_Int dataPacking,
                     FldArrayF& f, FldArrayI& c);
    E_Int readData75CE(FILE* ptrFile, E_Int ni, E_Int nj, E_Int nk,
                       E_Int dataPacking,
                       FldArrayF& f);
    E_Int readData75CE(FILE* ptrFile,
                       E_Int dataPacking,
                       FldArrayF& f, FldArrayI& c);
    E_Int readGeom108(FILE* ptrFile,
                      FldArrayF*& f);
    E_Int readGeom108CE(FILE* ptrFile,
                        FldArrayF*& f);
    void writeInt(E_Int value, E_Int isize,
                  FILE* ptrFile,
                  E_Bool convertEndian,
                  E_Int sizeInt, E_Int sizeLong);
    E_Int readInt(E_Int& value,
                  E_Int size,
                  FILE* ptrFile,
                  E_Bool convertEndian,
                  E_Int sizeInt, E_Int sizeLong);
    E_Int readIntTuple(FILE* ptrFile, E_Int& value);
    E_Int readIntTuple2(FILE* ptrFile, E_Int& value1, E_Int& value2);
    E_Int readIntTuple3(FILE* ptrFile, E_Int& value1, E_Int& value2,
      E_Int& value3);

    E_Int convertString2Int(char* str);

    E_Int findOpenablePlot3dFiles(
      char* file,
      char*& xyzFile, char*& qFile, char*& fFile);
    E_Int xyzread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field);
    E_Int p3dqread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field);
    E_Int p3dfread(
      char* file, char*& varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field);
    E_Int findWritablePlot3dFiles(
      char* file, char* varString,
      char*& xyzFile, char*& qFile, char*& fFile);
    E_Int xyzwrite(
      char* file, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, E_Int isize, E_Int rsize,
      E_Bool convertEndian);
    E_Int qwrite(
      char* file, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, E_Int isize, E_Int rsize,
      E_Bool convertEndian);
    E_Int ffwrite(
      char* file, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector<FldArrayF*>& field, E_Int isize, E_Int rsize,
      E_Bool convertEndian, E_Bool writeConservative);
    E_Int readWord(FILE* ptrFile, char* buf);
    E_Int readGivenKeyword(FILE* ptrFile, const char* keyword);
    E_Int readGivenKeyword(FILE* ptrFile, const char* keyword1,
      const char* keyword2);
    E_Int readGivenKeyword(FILE* ptrFile, const char* keyword1,
      const char* keyword2, const char* keyword3);
    E_Int readKeyword(FILE* ptrFile, char* buf);
    E_Int readDataAndKeyword(
      FILE* ptrFile, char* buf,
      std::list<const char*>& knownKeywords,
      char* prevData, char* keyword);
    E_Int readDouble(FILE* ptrFile, E_Float& value, E_Int formatLength=-1);
    E_Int readDouble(char* buf, E_Int size, E_Int& pos, E_Float& value);
    E_Int readInt(FILE* ptrFile, E_Int& value, E_Int formatLength=-1);
    E_Int readInt(char* buf, E_Int size, E_Int& pos, E_Int& value);
    E_Int readHexaInt(FILE* ptrFile, E_Int& value, E_Int formatLength=-1);
    void compressString(char* buf);
    E_Int getFormatLength(char* formatString);
    E_Int numeralVersion(char* version);
    E_Int skipComment(FILE*& ptrFile, char commentChar);
    E_Int skipLine(FILE*& ptrFile);
    E_Int readline(FILE*& ptrFile, char* buf, E_Int size);

    // Creation d'elements discrets
    void createPoint(std::vector<FldArrayF*>& unstructField,
                     std::vector<FldArrayI*>& connect,
                     std::vector<E_Int>& eltType,
                     E_Float density,
                     E_Float x, E_Float y, E_Float z);
    void createLine(std::vector<FldArrayF*>& structField,
                    std::vector<E_Int>& ni, std::vector<E_Int>& nj,
                    std::vector<E_Int>& nk,
                    E_Float density,
                    E_Float x1, E_Float y1, E_Float z1,
                    E_Float x2, E_Float y2, E_Float z2);

    // Traitement BCFaces
    E_Int getSizeOfBCFaces(PyObject* BCFaces);
    E_Int getBCFaces(PyObject* BCFaces, E_Int i,
                     char* name, FldArrayI& faces);

  private:
    /* Unique constructor. */
    GenIO();
    /* The destructor. */
    ~GenIO();

    // Should I convert endians?
    E_Bool _convertEndian;

    // Integers length (4 or 8)
    E_Int _intLength;

    // Real lentgh (4 or 8)
    E_Int _realLength;

    // The singleton instance.
    static GenIO* _instance;
};

}

#undef FldArrayF
#undef FldArrayI

#endif
