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

// Binary vtu file support
# include <inttypes.h>
# include <string.h>
# include <stdio.h>
# include <zlib.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

// Read a line: return 0 if OK, 1 otherwise
E_Int readLine(FILE* ptrFile, char* buf, const E_Int bufsize)
{
  E_Int i = 0;
  int c;

  while ((c = fgetc(ptrFile)) != EOF)
  {
    if (c == '\n') break;
    if (i < bufsize-1) buf[i++] = (char)c;
  }

  buf[i] = '\0';
  if (c == EOF && i == 0) return 1;
  return 0;
}

//===========================================================================
static const unsigned char b64Table[256] = {
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,62,64,64,64,63,
  52,53,54,55,56,57,58,59,60,61,64,64,64,64,64,64,
  64, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
  15,16,17,18,19,20,21,22,23,24,25,64,64,64,64,64,
  64,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
  41,42,43,44,45,46,47,48,49,50,51,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
  64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64
};

//===========================================================================
E_Int base64Decode(const char* in, E_Int inLen, unsigned char*& out)
{
  // skip leading whitespace
  while (inLen > 0 && (*in == ' ' || *in == '\t' || *in == '\r' || *in == '\n'))
  { in++; inLen--; }

  // trim trailing whitespace
  while (inLen > 0 && (in[inLen-1] == ' ' || in[inLen-1] == '\t' ||
                        in[inLen-1] == '\r' || in[inLen-1] == '\n'))
    inLen--;

  // printf("base64Decode inLen=" SF_D_ " first4='%.4s'\n", inLen, in);

  if (inLen <= 0 || inLen % 4 != 0)
  {
    printf("Warning: base64Decode: invalid input length=" SF_D_ "\n", inLen);
    out = NULL;
    return 0;
  }

  E_Int outLen = (inLen / 4) * 3;
  if (in[inLen-1] == '=') outLen--;
  if (in[inLen-2] == '=') outLen--;

  if (outLen <= 0)
  {
    printf("Warning: base64Decode: invalid output length=" SF_D_ "\n", outLen);
    out = NULL;
    return 0;
  }

  out = new unsigned char[(size_t)outLen];

  E_Int j = 0;
  for (E_Int i = 0; i < inLen; i += 4)
  {
    unsigned char b0 = b64Table[(unsigned char)in[i]];
    unsigned char b1 = b64Table[(unsigned char)in[i+1]];
    unsigned char b2 = b64Table[(unsigned char)in[i+2]];
    unsigned char b3 = b64Table[(unsigned char)in[i+3]];
    out[j++] = (b0 << 2) | (b1 >> 4);
    if (j < outLen) out[j++] = (b1 << 4) | (b2 >> 2);
    if (j < outLen) out[j++] = (b2 << 6) | b3;
  }
  return outLen;
}

//===========================================================================
// Read appended data block, handling both compressed and uncompressed cases.
// Returns a newly allocated buffer in data and its size in bytes in dataSize.
// Caller is responsible for freeing data.
void readAppendedBlock(FILE* ptrFile, const off_t base, const off_t offset,
  const E_Bool compressed, const E_Bool inlineBinary, const char* headerType,
  void*& data, E_Int& dataSize)
{
  data = NULL; dataSize = 0;
  printf("DEBUG: base=%jd offset=%jd seek=%jd\n",
    (intmax_t)base,
    (intmax_t)offset,
    (intmax_t)(base + offset));

  if (inlineBinary)
  {
    // read full line dynamically
    std::vector<char> lineBufVec;
    lineBufVec.reserve(1024*1024);
    int c;
    while ((c = fgetc(ptrFile)) != EOF && c != '\n')
      lineBufVec.push_back((char)c);
    lineBufVec.push_back('\0');
    E_Int len = (E_Int)lineBufVec.size() - 1;

    // strip trailing spaces and spaces
    while (len > 0 && (lineBufVec[len-1] == '\r' || lineBufVec[len-1] == ' '))
    { lineBufVec[len-1] = '\0'; len--; }

    // strip trailing XML tag
    for (E_Int i = 0; i < len; i++)
    {
      if (lineBufVec[i] == '<') { lineBufVec[i] = '\0'; len = i; break; }
    }

    // strip all whitespace from line to get a clean base64 string
    std::vector<char> cleanBuf;
    cleanBuf.reserve(len);
    for (E_Int i = 0; i < len; i++)
    {
      char ch = lineBufVec[i];
      if (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n')
        cleanBuf.push_back(ch);
    }
    cleanBuf.push_back('\0');
    const char* lineBuf = cleanBuf.data();
    len = (E_Int)cleanBuf.size() - 1;
    // printf("clean base64 line length=" SF_D_ " first8='%.8s'\n", len, lineBuf);

    if (!compressed)
    {
      // header = single uint32/uint64 = uncompressed byte count
      E_Int headerB64Len = (strcmp(headerType, "UInt64") == 0) ? 12 : 8;
      unsigned char* decodedHeader = NULL;
      base64Decode(lineBuf, headerB64Len, decodedHeader);
      if (decodedHeader == NULL) { printf("Warning: vturead: failed to decode header.\n"); return; }
      if (strcmp(headerType, "UInt64") == 0) dataSize = (E_Int)(*(uint64_t*)decodedHeader);
      else dataSize = (E_Int)(*(uint32_t*)decodedHeader);
      delete [] decodedHeader;
      printf("uncompressed dataSize=" SF_D_ "\n", dataSize);

      unsigned char* decodedData = NULL;
      E_Int dataBlockLen = len - headerB64Len;
      while (dataBlockLen % 4 != 0 && dataBlockLen > 0) dataBlockLen--;
      base64Decode(lineBuf + headerB64Len, dataBlockLen, decodedData);
      if (decodedData == NULL) { printf("Warning: vturead: failed to decode data.\n"); return; }
      data = decodedData;
    }
    else
    {
      // header = [nblocks][uncompBlockSize][lastBlockSize][compSize0]...[compSizeN]
      // decode 3 words to get nblocks
      E_Int headerFixedB64Len = (strcmp(headerType, "UInt64") == 0) ? 32 : 16;
      unsigned char* decodedFixed = NULL;
      base64Decode(lineBuf, headerFixedB64Len, decodedFixed);
      if (decodedFixed == NULL) { printf("Warning: vturead: failed to decode fixed header.\n"); return; }
      uint32_t nblocks = *(uint32_t*)decodedFixed;
      delete [] decodedFixed;

      // decode header: (3 + nblocks) * sizeof(uint32_t) bytes
      E_Int headerBytes = (3 + nblocks) * (E_Int)sizeof(uint32_t);
      E_Int headerB64Len = ((headerBytes + 2) / 3) * 4;
      unsigned char* decodedHeader = NULL;
      base64Decode(lineBuf, headerB64Len, decodedHeader);
      if (decodedHeader == NULL) { printf("Warning: vturead: failed to decode header.\n"); return; }

      // decode data block
      unsigned char* decodedData = NULL;
      E_Int dataBlockLen = len - headerB64Len;
      while (dataBlockLen % 4 != 0 && dataBlockLen > 0) dataBlockLen--;
      base64Decode(lineBuf + headerB64Len, dataBlockLen, decodedData);
      if (decodedData == NULL) { printf("Warning: vturead: failed to decode data.\n"); delete [] decodedHeader; return; }

      uint32_t* h = (uint32_t*)decodedHeader;
      uint32_t uncompBlockSize = h[1];
      uint32_t lastBlockSize = h[2];
      printf("nblocks=%" PRIu32 " uncompBlockSize=%" PRIu32 " lastBlockSize=%" PRIu32 "\n",
             nblocks, uncompBlockSize, lastBlockSize);

      dataSize = (nblocks - 1) * (E_Int)uncompBlockSize + (E_Int)lastBlockSize;
      data = new unsigned char[dataSize];
      unsigned char* outPtr = (unsigned char*)data;
      unsigned char* inPtr = decodedData;

      for (uint32_t b = 0; b < nblocks; b++)
      {
        uint32_t compSize = h[3 + b];
        uLongf destLen = (b < nblocks - 1) ? uncompBlockSize : lastBlockSize;
        int ierr = uncompress(outPtr, &destLen, inPtr, compSize);
        if (ierr != Z_OK)
          printf("Warning: vturead: zlib uncompress failed (block=%" PRIu32 " err=%d).\n", b, ierr);
        outPtr += destLen;
        inPtr += compSize;
      }
      delete [] decodedHeader;
      delete [] decodedData;
    }
  }
  else  // appended binary
  {
    fseeko(ptrFile, base + offset, SEEK_SET);

    if (!compressed)
    {
      if (strcmp(headerType, "Int64") == 0)
      {
        uint64_t byteCount;
        fread(&byteCount, sizeof(uint64_t), 1, ptrFile);
        dataSize = (E_Int)byteCount;
      }
      else
      {
        uint32_t byteCount;
        fread(&byteCount, sizeof(uint32_t), 1, ptrFile);
        dataSize = (E_Int)byteCount;
      }
      data = new unsigned char[dataSize];
      fread(data, 1, dataSize, ptrFile);
    }
    else
    {
      // [nblocks][uncompBlockSize][lastBlockSize][compSize0]...[data0]...
      uint32_t nblocks, uncompBlockSize, lastBlockSize;
      fread(&nblocks, sizeof(uint32_t), 1, ptrFile);
      fread(&uncompBlockSize, sizeof(uint32_t), 1, ptrFile);
      fread(&lastBlockSize, sizeof(uint32_t), 1, ptrFile);
      printf("nblocks%" PRIu32 " uncompBlockSize=%" PRIu32 " lastBlockSize=%" PRIu32 "\n",
             nblocks, uncompBlockSize, lastBlockSize);

      uint32_t* compSizes = new uint32_t[nblocks];
      fread(compSizes, sizeof(uint32_t), nblocks, ptrFile);

      dataSize = (nblocks - 1) * (E_Int)uncompBlockSize + (E_Int)lastBlockSize;
      data = new unsigned char[dataSize];
      unsigned char* outPtr = (unsigned char*)data;

      for (uint32_t b = 0; b < nblocks; b++)
      {
        unsigned char* compBuf = new unsigned char[compSizes[b]];
        fread(compBuf, 1, compSizes[b], ptrFile);
        uLongf destLen = (b < nblocks - 1) ? uncompBlockSize : lastBlockSize;
        int ierr = uncompress(outPtr, &destLen, compBuf, compSizes[b]);
        if (ierr != Z_OK)
          printf("Warning: vturead: zlib uncompress failed (block=%" PRIu32 " err=%d).\n", b, ierr);
        outPtr += destLen;
        delete [] compBuf;
      }
      delete [] compSizes;
    }
  }
}

//========================================================================
// read POINTS
void readPointsVTU(
  FILE* ptrFile, const E_Bool changeEndian, const E_Bool compressed,
  const char* headerType, const E_Int npts, E_Bool& formatted, double*& coords
)
{
  const E_Int bufsize = 1024;
  char buf[bufsize]; char type[64] = "";
  char format[64] = "";
  off_t offset = 0;
  char* p;

  // find <Points>
  while (readLine(ptrFile, buf, bufsize) == 0)
    if (strstr(buf, "<Points>") != NULL) break;

  // read DataArray line
  readLine(ptrFile, buf, bufsize);
  printf("DataArray line: %s\n", buf);
  p = strstr(buf, "type=\"");
  if (p != NULL) sscanf(p, "type=\"%[^\"]\"", type);
  p = strstr(buf, "format=\"");
  if (p != NULL) sscanf(p, "format=\"%[^\"]\"", format);
  p = strstr(buf, "offset=\"");
  if (p != NULL) { intmax_t tmp; sscanf(p, "offset=\"%jd\"", &tmp); offset = (off_t)tmp; }
  E_Bool isDoublePrecision = strcmp(type, "Float64") == 0;
  E_Bool inlineBinary = strcmp(format, "binary") == 0;
  printf("VTU type=%s format=%s\n", type, format);

  if (strcmp(format, "ascii") == 0)
  {
    formatted = true;
    for (E_Int i = 0; i < npts; i++)
      fscanf(ptrFile, "%lf %lf %lf", &coords[3*i], &coords[3*i+1], &coords[3*i+2]);
    // skip until </DataArray>
    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "</DataArray>") != NULL) break;
    // skip </Points>
    readLine(ptrFile, buf, bufsize);
  }
  else if (strcmp(format, "appended") == 0 || strcmp(format, "binary") == 0)
  {
    formatted = false;
    off_t base = 0;
    if (!inlineBinary)
    {
      // find appended section and base
      while (readLine(ptrFile, buf, bufsize) == 0)
        if (strstr(buf, "<AppendedData")) break;
      while (fgetc(ptrFile) != '_');
      base = ftello(ptrFile);
    }
    // read raw, appended, or inline binary block
    void* data; E_Int dataSize;
    readAppendedBlock(ptrFile, base, offset,
                      compressed, inlineBinary, headerType,
                      data, dataSize);

    // fill coords from data buffer
    if (isDoublePrecision)
    {
      double* dbuf = (double*)data;
      if (changeEndian)
        for (E_Int i = 0; i < 3*npts; i++) coords[i] = DBE(dbuf[i]);
      else
        for (E_Int i = 0; i < 3*npts; i++) coords[i] = dbuf[i];
    }
    else
    {
      float* fbuf = (float*)data;
      if (changeEndian)
        for (E_Int i = 0; i < 3*npts; i++) coords[i] = FBE(fbuf[i]);
      else
        for (E_Int i = 0; i < 3*npts; i++) coords[i] = (double)fbuf[i];
    }
    delete [] (unsigned char*)data;
  }
  printf("point0 %lf %lf %lf\n", coords[0], coords[1], coords[2]);
  // printf("point1 %lf %lf %lf\n", coords[3], coords[4], coords[5]);
  // printf("point2 %lf %lf %lf\n", coords[6], coords[7], coords[8]);
}

//===========================================================================
// read CELLS
void readCellsVTU(
  FILE* ptrFile,
  const E_Bool changeEndian, const E_Bool compressed, const E_Bool formatted,
  const char* headerType, const E_Int ntotElts,
  E_Int*& cells, E_Int*& vtuOffsets, std::vector<E_Int>& nepc, char*& eltType,
  std::vector<E_Int>& eltMap, std::vector<E_Int>& offsets
)
{
  const E_Int bufsize = 1024;
  char buf[bufsize];
  off_t offConn = 0, offOffsets = 0, offTypes = 0;
  char formatConn[64] = "";
  char* p;
  cells = NULL; vtuOffsets = NULL;

  const E_Int nmaxVTKTypes = 32;
  std::vector<E_Int> tmp_nepc(nmaxVTKTypes, 0);
  const char* vtkTypeNames[nmaxVTKTypes] = {
    NULL,       // 0  unused
    NULL,       // 1  unused
    NULL,       // 2  unused
    "BAR",      // 3  VTK_LINE
    NULL,       // 4  unused
    "TRI",      // 5  VTK_TRIANGLE
    NULL,       // 6  unused
    NULL,       // 7  unused
    NULL,       // 8  unused
    "QUAD",     // 9  VTK_QUAD
    "TETRA",    // 10 VTK_TETRA
    NULL,       // 11 unused
    "HEXA",     // 12 VTK_HEXAHEDRON
    "PENTA",    // 13 VTK_WEDGE
    "PYRA",     // 14 VTK_PYRAMID
    NULL,       // 15 unused
    NULL,       // 16 unused
    NULL,       // 17 unused
    NULL,       // 18 unused
    NULL,       // 19 unused
    NULL,       // 20 unused
    NULL,       // 21 unused
    "TRI_6",    // 22 VTK_QUADRATIC_TRIANGLE
    "QUAD_8",   // 23 VTK_QUADRATIC_QUAD
    "TETRA_10", // 24 VTK_QUADRATIC_TETRA
    "HEXA_20",  // 25 VTK_QUADRATIC_HEXAHEDRON
    NULL,       // 26 unused
    NULL,       // 27 unused
    NULL,       // 28 unused
    NULL,       // 29 unused
    NULL,       // 30 unused
    NULL        // 31 unused
  };

  // find <Cells>
  rewind(ptrFile);
  while (readLine(ptrFile, buf, bufsize) == 0)
    if (strstr(buf, "<Cells")) break;

  // read connectivity DataArray header line
  while (readLine(ptrFile, buf, bufsize) == 0)
    if (strstr(buf, "Name=\"connectivity\"")) break;
  printf("conn header: %s\n", buf);
  p = strstr(buf, "format=\"");
  if (p != NULL) sscanf(p, "format=\"%[^\"]\"", formatConn);
  p = strstr(buf, "offset=\"");
  if (p != NULL) { intmax_t tmp; sscanf(p, "offset=\"%jd\"", &tmp); offConn = (off_t)tmp; }
  char connType[64] = "";
  p = strstr(buf, "type=\"");
  if (p != NULL) sscanf(p, "type=\"%[^\"]\"", connType);
  E_Bool inlineBinary = strcmp(formatConn, "binary") == 0;
  printf("inlineBinary=%d connType=%s\n", inlineBinary, connType);

  void* data; E_Int dataSize;
  if (inlineBinary)
  {
    // read connectivity data line
    readAppendedBlock(ptrFile, 0, 0, compressed, true, headerType, data, dataSize);
    printf("connectivity dataSize=" SF_D_ "\n", dataSize);
    if (strcmp(connType, "Int64") == 0)
    {
      E_Int cnsize = dataSize/sizeof(int64_t);
      if (sizeof(E_Int) == sizeof(int64_t))
      {
        cells = (E_Int*)data;
        if (changeEndian)
          for (E_Int i = 0; i < cnsize; i++) cells[i] = IBE(cells[i]);
      }
      else  // Int64 to Int32
      {
        int64_t* tmp64 = (int64_t*)data;
        cells = new E_Int[cnsize];
        for (E_Int i = 0; i < cnsize; i++)
        {
          int64_t v = tmp64[i];
          if (changeEndian) v = IBE(v);
          cells[i] = (E_Int)v;
        }
        delete [] (unsigned char*)data;
      }
    }
    else
    {
      E_Int cnsize = dataSize/sizeof(int32_t);
      if (sizeof(E_Int) == sizeof(int32_t))
      {
        cells = (E_Int*)data;
        if (changeEndian)
          for (E_Int i = 0; i < cnsize; i++) cells[i] = IBE(cells[i]);
      }
      else  // Int32 to Int64
      {
        int32_t* tmp32 = (int32_t*)data;
        cells = new E_Int[cnsize];
        for (E_Int i = 0; i < cnsize; i++)
        {
          int32_t v = tmp32[i];
          if (changeEndian) v = IBE(v);
          cells[i] = (E_Int)v;
        }
        delete [] (unsigned char*)data;
      }
    }
    printf("first connectivity=" SF_D4_ "\n", cells[0], cells[1], cells[2], cells[4]);

    // read offsets DataArray header line
    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "Name=\"offsets\"")) break;
    printf("offsets header: %s\n", buf);
    char offType[64] = "";
    p = strstr(buf, "type=\"");
    if (p != NULL) sscanf(p, "type=\"%[^\"]\"", offType);
    printf("offsets type=%s\n", offType);

    // read offsets data line
    void* offData; E_Int offDataSize;
    readAppendedBlock(ptrFile, 0, 0, compressed, true, headerType, offData, offDataSize);
    printf("offsets dataSize=" SF_D_ "\n", offDataSize);
    vtuOffsets = new E_Int[ntotElts];
    if (strcmp(offType, "Int64") == 0)
    {
      if (sizeof(E_Int) == sizeof(int64_t))
      {
        vtuOffsets = (E_Int*)offData;
        if (changeEndian)
          for (E_Int i = 0; i < ntotElts; i++) vtuOffsets[i] = IBE(vtuOffsets[i]);
      }
      else  // Int64 to Int32
      {
        int64_t* offTmp64 = (int64_t*)offData;
        vtuOffsets = new E_Int[ntotElts];
        for (E_Int i = 0; i < ntotElts; i++)
        {
          int64_t v = offTmp64[i];
          if (changeEndian) v = IBE(v);
          vtuOffsets[i] = (E_Int)v;
        }
        delete [] (unsigned char*)offData;
      }
    }
    else
    {
      if (sizeof(E_Int) == sizeof(int32_t))
      {
        vtuOffsets = (E_Int*)offData;
        if (changeEndian)
          for (E_Int i = 0; i < ntotElts; i++) vtuOffsets[i] = IBE(vtuOffsets[i]);
      }
      else  // Int32 to Int64
      {
        int32_t* offTmp32 = (int32_t*)offData;
        vtuOffsets = new E_Int[ntotElts];
        for (E_Int i = 0; i < ntotElts; i++)
        {
          int32_t v = offTmp32[i];
          if (changeEndian) v = IBE(v);
          vtuOffsets[i] = (E_Int)v;
        }
        delete [] (unsigned char*)offData;
      }
    }

    // read types DataArray header line
    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "Name=\"types\"")) break;
    printf("types header: %s\n", buf);

    // read types data line
    readAppendedBlock(ptrFile, 0, 0, compressed, true, headerType, data, dataSize);
    printf("types dataSize=" SF_D_ "\n", dataSize);
  }
  else
  {
    // appended: read offsets and types header lines for their offsets
    intmax_t tmp;
    p = strstr(buf, "offset=\"");
    if (p != NULL) { sscanf(p, "offset=\"%jd\"", &tmp); offConn = (off_t)tmp; }

    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "Name=\"offsets\"")) break;
    printf("offsets header: %s\n", buf);
    char offType[64] = "";
    p = strstr(buf, "type=\"");
    if (p != NULL) sscanf(p, "type=\"%[^\"]\"", offType);
    p = strstr(buf, "offset=\"");
    if (p != NULL) { sscanf(p, "offset=\"%jd\"", &tmp); offOffsets = (off_t)tmp; }

    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "Name=\"types\"")) break;
    printf("types header: %s\n", buf);
    p = strstr(buf, "offset=\"");
    if (p != NULL) { sscanf(p, "offset=\"%jd\"", &tmp); offTypes = (off_t)tmp; }

    // find appended data section
    while (readLine(ptrFile, buf, bufsize) == 0)
      if (strstr(buf, "<AppendedData")) break;
    while (fgetc(ptrFile) != '_');
    off_t base = ftello(ptrFile);
    printf("base=%jd\n", (intmax_t)base);

    // read connectivity
    readAppendedBlock(ptrFile, base, offConn, compressed, false, headerType, data, dataSize);
    printf("connectivity dataSize=" SF_D_ "\n", dataSize);
    if (strcmp(connType, "Int64") == 0)
    {
      E_Int cnsize = dataSize/sizeof(int64_t);
      if (sizeof(E_Int) == sizeof(int64_t))
      {
        cells = (E_Int*)data;
        if (changeEndian)
          for (E_Int i = 0; i < cnsize; i++) cells[i] = IBE(cells[i]);
      }
      else  // Int64 to Int32
      {
        int64_t* tmp64 = (int64_t*)data;
        cells = new E_Int[cnsize];
        for (E_Int i = 0; i < cnsize; i++)
        {
          int64_t v = tmp64[i];
          if (changeEndian) v = IBE(v);
          cells[i] = (E_Int)v;
        }
        delete [] (unsigned char*)data;
      }
    }
    else
    {
      E_Int cnsize = dataSize/sizeof(int32_t);
      if (sizeof(E_Int) == sizeof(int32_t))
      {
        cells = (E_Int*)data;
        if (changeEndian)
          for (E_Int i = 0; i < cnsize; i++) cells[i] = IBE(cells[i]);
      }
      else  // Int32 to Int64
      {
        // convert 
        int32_t* tmp32 = (int32_t*)data;
        cells = new E_Int[cnsize];
        for (E_Int i = 0; i < cnsize; i++)
        {
          int32_t v = tmp32[i];
          if (changeEndian) v = IBE(v);
          cells[i] = (E_Int)v;
        }
        delete [] (unsigned char*)data;
      }
    }
    printf("first connectivity=" SF_D4_ "\n", cells[0], cells[1], cells[2], cells[4]);

    // read offsets
    void* offData; E_Int offDataSize;
    readAppendedBlock(ptrFile, base, offOffsets, compressed, false, headerType, offData, offDataSize);
    printf("offsets dataSize=" SF_D_ "\n", offDataSize);
    vtuOffsets = new E_Int[ntotElts];
    if (strcmp(offType, "Int64") == 0)
    {
      if (sizeof(E_Int) == sizeof(int64_t))
      {
        vtuOffsets = (E_Int*)offData;
        if (changeEndian)
          for (E_Int i = 0; i < ntotElts; i++) vtuOffsets[i] = IBE(vtuOffsets[i]);
      }
      else  // Int64 to Int32
      {
        int64_t* offTmp64 = (int64_t*)offData;
        vtuOffsets = new E_Int[ntotElts];
        for (E_Int i = 0; i < ntotElts; i++)
        {
          int64_t v = offTmp64[i];
          if (changeEndian) v = IBE(v);
          vtuOffsets[i] = (E_Int)v;
        }
        delete [] (unsigned char*)offData;
      }
    }
    else
    {
      if (sizeof(E_Int) == sizeof(int32_t))
      {
        vtuOffsets = (E_Int*)offData;
        if (changeEndian)
          for (E_Int i = 0; i < ntotElts; i++) vtuOffsets[i] = IBE(vtuOffsets[i]);
      }
      else  // Int32 to Int64
      {
        int32_t* offTmp32 = (int32_t*)offData;
        vtuOffsets = new E_Int[ntotElts];
        for (E_Int i = 0; i < ntotElts; i++)
        {
          int32_t v = offTmp32[i];
          if (changeEndian) v = IBE(v);
          vtuOffsets[i] = (E_Int)v;
        }
        delete [] (unsigned char*)offData;
      }
    }

    // read types
    readAppendedBlock(ptrFile, base, offTypes, compressed, false, headerType, data, dataSize);
    printf("types dataSize=" SF_D_ "\n", dataSize);
  }

  // process types (common to both paths)
  unsigned char* tmp = (unsigned char*)data;
  std::vector<E_Int> cellTypes(ntotElts);
  for (E_Int i = 0; i < ntotElts; i++)
  {
    cellTypes[i] = (E_Int)tmp[i];
    tmp_nepc[cellTypes[i]]++;
  }
  delete [] (unsigned char*)data;
  printf("cellTypes[0]=" SF_D_ "\n", cellTypes[0]);

  // set eltType and compress nepc
  nepc.clear(); nepc.reserve(nmaxVTKTypes);
  eltType[0] = '\0';
  for (E_Int i = 0; i < nmaxVTKTypes; i++)
  {
    if (tmp_nepc[i] > 0)
    {
      nepc.push_back(tmp_nepc[i]);
      if (eltType[0] == '\0') strcpy(eltType, vtkTypeNames[i]);
      else
      {
        strcat(eltType, ",");
        strcat(eltType, vtkTypeNames[i]);
      }
    }
  }
  printf("eltType=%s\n", eltType);

  // build starting offsets for each type
  E_Int tmp_offsets[nmaxVTKTypes]; tmp_offsets[0] = 0;
  for (E_Int i = 1; i < nmaxVTKTypes; i++)
    tmp_offsets[i] = tmp_offsets[i-1] + tmp_nepc[i-1];

  // fill element map: eltMap[sorted_position] = original_index
  eltMap.clear(); eltMap.resize(ntotElts);
  for (E_Int i = 0; i < ntotElts; i++)
  {
    E_Int ct = cellTypes[i];
    eltMap[tmp_offsets[ct]++] = i;
  }

  // compress offsets
  offsets.clear(); offsets.resize(nepc.size() + 1); offsets[0] = 0;
  E_Int nc = 0;
  for (E_Int i = 0; i < nmaxVTKTypes; i++)
  {
    if (tmp_nepc[i] > 0)
    {
      offsets[nc+1] = offsets[nc] + tmp_nepc[i];
      nc++;
    }
  }
  printf("readCellsVTU done\n");
}

//========================================================================
void readScalarVTU(FILE* ptrFile, const E_Bool changeEndian, const E_Bool formatted,
                   const E_Int npts, char* varName, float*& fieldf, double*& fieldd)
{ 
  // fieldf = NULL; fieldd = NULL;

  // char buf[256];

  // // recupere dataName, dataType, numComp
  // char dataType[256]; E_Int numComp;
  // fscanf(ptrFile, " %s %s " SF_D_ "\n", varName, dataType, &numComp);
  // printf("varName=%s, dataType=%s, numCom=" SF_D_ "\n", varName, dataType, numComp);

  // // Passe lookup table
  // readLine(ptrFile, buf);

  // if (strcmp(dataType, "float") == 0)
  // {
  //   fieldf = new float [npts];
  //   if (formatted == false)
  //   {
  //       fread(fieldf, sizeof(float), npts, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%f", &fieldf[i]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < npts; i++) fieldf[i] = FBE(fieldf[i]);
  //   }
  // }
  // else if (strcmp(dataType, "double") == 0)
  // {
  //   fieldd = new double [npts];
  //   if (formatted == false)
  //   {
  //       fread(fieldd, sizeof(double), npts, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%lf", &fieldd[i]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < npts; i++) fieldd[i] = DBE(fieldd[i]);
  //   }
  // }
}

//==========================================================================
void readVectorVTU(FILE* ptrFile, const E_Bool changeEndian, const E_Bool formatted,
                   const E_Int npts, char* varName, float*& fieldf, double*& fieldd)
{ 
  // fieldf = NULL; fieldd = NULL;

  // // recupere dataName, dataType
  // char dataType[256];
  // fscanf(ptrFile, " %s %s\n", varName, dataType);
  // printf("varName=%s, dataType=%s\n", varName, dataType);

  // if (strcmp(dataType, "float") == 0)
  // {
  //   fieldf = new float [3*npts];
  //   if (formatted == false)
  //   {
  //       fread(fieldf, sizeof(float), 3*npts, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%f %f %f", &fieldf[3*i], &fieldf[3*i+1], &fieldf[3*i+2]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < 3*npts; i++) fieldf[i] = FBE(fieldf[i]);
  //   }
  // }
  // else if (strcmp(dataType, "double") == 0)
  // {
  //   fieldd = new double [3*npts];
  //   if (formatted == false)
  //   {
  //       fread(fieldd, sizeof(double), 3*npts, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%lf %lf %lf", &fieldd[3*i], &fieldd[3*i+1], &fieldd[3*i+2]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < 3*npts; i++) fieldd[i] = DBE(fieldd[i]);
  //   }
  // }
}

//========================================================================
void readFieldVTU(FILE* ptrFile, E_Bool changeEndian, E_Bool formatted,
                  E_Int& npts, E_Int& nfields, char* varName, float*& fieldf, double*& fieldd)
{ 
  // fieldf = NULL; fieldd = NULL;

  // // recupere dataName, numArrays
  // char dataName[256]; E_Int numArrays;
  // fscanf(ptrFile, " %s " SF_D_ "\n", dataName, &numArrays);
  // printf("dataName=%s, numArrays=" SF_D_ "\n", dataName, numArrays);
  
  // // Lit uniquement la premiere variable
  // E_Int numComponents; E_Int numTuples; char dataType[256];
  // fscanf(ptrFile, "%s " SF_D2_ " %s\n", varName, &numComponents, &numTuples, dataType);
  // npts = numTuples;
  // nfields = numComponents;
  // E_Int size = npts*nfields;
  // printf("varName=%s, dataType=%s, npts=" SF_D_ "\n", varName, dataType, npts);

  // if (strcmp(dataType, "float") == 0)
  // {
  //   fieldf = new float [size];
  //   if (formatted == false)
  //   {
  //       fread(fieldf, sizeof(float), size, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < size; i++) fscanf(ptrFile, "%f", &fieldf[i]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < size; i++) fieldf[i] = FBE(fieldf[i]);
  //   }
  // }
  // else if (strcmp(dataType, "double") == 0)
  // {
  //   fieldd = new double [size];
  //   if (formatted == false)
  //   {
  //       fread(fieldd, sizeof(double), size, ptrFile);
  //       fgetc(ptrFile); // avoid \n
  //   }
  //   else
  //   {
  //       for (E_Int i = 0; i < size; i++) fscanf(ptrFile, "%lf", &fieldd[i]);
  //   }
  //   if (changeEndian)
  //   {
  //       for (E_Int i = 0; i < size; i++) fieldd[i] = DBE(fieldd[i]);
  //   }
  // }
}

//=============================================================================
// binvturead
//=============================================================================
E_Int K_IO::GenIO::binvturead(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType,
  vector<char*>& zoneNames,
  E_Int api
)
{
  // Open file
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");
  if (ptrFile == NULL)
  {
    printf("Warning: vturead: cannot open file %s.\n", file);
    return 1;
  }
  const E_Int bufsize = 1024;
  char buf[bufsize];

  // Read XML declaration
  readLine(ptrFile, buf, bufsize);
  printf("XML: %s\n", buf);

  // Read VTKFile line
  readLine(ptrFile, buf, bufsize);
  printf("VTKFile: %s\n", buf);

  // parse byte_order and compressor independently to handle arbitrary attribute ordering
  char byteOrder[64] = "LittleEndian";
  char fileType[64]  = "";
  char headerType[64] = "UInt32";
  char* p;
  p = strstr(buf, "type=\"");
  if (p != NULL) sscanf(p, "type=\"%[^\"]\"", fileType);
  p = strstr(buf, "byte_order=\"");
  if (p != NULL) sscanf(p, "byte_order=\"%[^\"]\"", byteOrder);
  p = strstr(buf, "header_type=\"");
  if (p != NULL) sscanf(p, "header_type=\"%[^\"]\"", headerType);
  
  printf("type=%s\n", fileType);
  printf("byte_order=%s\n", byteOrder);
  printf("headerType=%s\n", headerType);

  // Endianness
  E_Bool changeEndian = strcmp(byteOrder, "BigEndian") == 0;
  printf("changeEndian=%d\n", changeEndian);

  // Compression
  E_Bool compressed = strstr(buf, "compressor=") != NULL;
  if (compressed &&
      strstr(buf, "ZLib") == NULL &&
      strstr(buf, "zlib") == NULL &&
      strstr(buf, "ZLIB") == NULL)
  {
    printf("Warning: vturead: unsupported compressor, only zlib is supported.\n");
    fclose(ptrFile);
    return 1;
  }
  printf("compressed=%d\n", compressed);

  // Find Piece line, then number of vertices and elements
  while (readLine(ptrFile, buf, bufsize) == 0)
    if (strstr(buf, "<Piece")) break;
  // printf("Piece: %s\n", buf);
  E_Int npts = 0;
  E_Int ntotElts = 0;
  p = strstr(buf, "NumberOfPoints=\"");
  if (p != NULL) sscanf(p, "NumberOfPoints=\"" SF_D_ "\"", &npts);
  p = strstr(buf, "NumberOfCells=\"");
  if (p != NULL) sscanf(p, "NumberOfCells=\"" SF_D_ "\"", &ntotElts);
  printf("npts=" SF_D_ ", ncells=" SF_D_ "\n", npts, ntotElts);

  // Read vertices and whether this is formatted or binary
  E_Bool formatted = true;  // default unless appended/binary detected
  double* coords = new double[3*npts];
  readPointsVTU(ptrFile, changeEndian, compressed, headerType, npts, formatted, coords);

  // Read elements
  E_Int* cells; E_Int* vtuOffsets;
  char* eltString = new char[K_ARRAY::VARSTRINGLENGTH];
  std::vector<E_Int> nepc, eltMap, offsets;
  readCellsVTU(ptrFile, changeEndian, compressed, formatted, headerType,
               ntotElts, cells, vtuOffsets, nepc, eltString, eltMap, offsets);

  // Close file
  fclose(ptrFile);

  // // print coords
  // printf("coords:\n");
  // for (E_Int i = 0; i < npts; i++)
  //   printf("  point%d: %lf %lf %lf\n", i, coords[3*i], coords[3*i+1], coords[3*i+2]);

  // // print cells
  // E_Int cnsize = 0;
  // for (E_Int i = 0; i < (E_Int)nepc.size(); i++) cnsize += nepc[i] * /* nvpe */ 1;
  // printf("ntotElts=%d\n", ntotElts);
  // printf("eltString=%s\n", eltString);

  // // print nepc
  // printf("nepc (number of elements per connectivity, size=%d):\n", (int)nepc.size());
  // for (E_Int i = 0; i < (E_Int)nepc.size(); i++)
  //   printf("  nepc[%d]=%d\n", i, nepc[i]);

  // // print offsets
  // printf("offsets (size=%d):\n", (int)offsets.size());
  // for (E_Int i = 0; i < (E_Int)offsets.size(); i++)
  //   printf("  offsets[%d]=%d\n", i, offsets[i]);

  // // print eltMap
  // printf("eltMap (size=%d):\n", (int)eltMap.size());
  // for (E_Int i = 0; i < (E_Int)eltMap.size(); i++)
  //   printf("  eltMap[%d]=%d\n", i, eltMap[i]);

  // // print cells — we need total connectivity size
  // printf("cells:\n");
  // E_Int cellIdx = 0;
  // for (E_Int ic = 0; ic < (E_Int)nepc.size(); ic++)
  // {
  //   // find the vtk type for this connectivity block
  //   // eltString contains comma-separated names, nepc[ic] elements
  //   printf("  connectivity block %d (%d elements):\n", ic, nepc[ic]);
  //   for (E_Int i = 0; i < nepc[ic]; i++)
  //   {
  //     printf("    cell%d: ", i);
  //     // we don't know nvpe here without parsing eltString
  //     // print raw cells starting from cellIdx
  //     printf("cells[%d]=%d\n", cellIdx, cells[cellIdx]);
  //     cellIdx++;
  //   }
  // }

  // Build ME connectivity
  E_Int nc = nepc.size();
  std::vector<E_Int> nvpe(nc), eltIds(nc);
  std::vector<char*> eltStrings;
  K_ARRAY::extractVars(eltString, eltStrings);
  char dummyStr[128]; E_Int dummy;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_ARRAY::eltString2TypeId(eltStrings[ic], dummyStr, nvpe[ic],
                              dummy, eltIds[ic]);
  }

  E_Int nfld = 3;
  varString = new char[K_ARRAY::VARSTRINGLENGTH];
  strcpy(varString, "x,y,z");

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nepc,
                                       eltString, 0, api); 
  FldArrayI* cn2; FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  #pragma omp parallel
  {
    E_Int elt, ind, off;
    for (E_Int n = 0; n < nfld; n++)
    {
      E_Float* f2p = f2->begin(n+1);
      #pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = coords[3*i+n];
    }

    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm2 = *(cn2->getConnect(ic));
      off = offsets[ic];
      #pragma omp for
      for (E_Int i = 0; i < nepc[ic]; i++)
      {
        elt = eltMap[off+i];  // original element index
        if (elt == 0) ind = 0;
        else ind = vtuOffsets[elt-1];
        // std::cout << ic << ", " << i << ", " << elt << ", " << ind << std::endl;
        for (E_Int j = 0; j < nvpe[ic]; j++) cm2(i,j+1) = (E_Int)cells[ind+j] + 1;  // starts at 1
      }
    }
  }

  // for (E_Int ic = 0; ic < nc; ic++)
  // {
  //   FldArrayI& cm2 = *(cn2->getConnect(ic));
  //   E_Int off = offsets[ic];
  //   std::cout << "offsets[ic] = " << offsets[ic] << std::endl;
  //   // for (E_Int i = 0; i < nepc[ic]; i++)
  //   // {
  //   //   for (E_Int j = 0; j < nvpe[ic]; j++) std::cout << "cm2("<<i<<","<<j+1<<")" << cm2(i,j+1) << std::endl;
  //   // }
  //   // std::cout << "eltIds[ic] = " << eltIds[ic] << std::endl;
  // }

  eltType.clear();
  if (api == 3)
  {
    eltType.resize(1);
    for (E_Int ic = 0; ic < nc; ic++) eltType[0].push_back(eltIds[ic]);
    for (E_Int ic = 0; ic < nc; ic++) std::cout << eltType[0][ic] << std::endl;
    unstructField.push_back(f2);
    connect.push_back(cn2);
  }
  else
  {
    eltType.resize(nc);
    for (E_Int ic = 0; ic < nc; ic++)
    {
      unstructField.push_back(f2);
      connect.push_back(cn2->getConnect(ic));
      eltType[ic].push_back({eltIds[ic]});
    }
  }

  // Free memory
  delete [] coords; delete [] cells; delete [] eltString;
  delete [] vtuOffsets;
  for (size_t ic = 0; ic < eltStrings.size(); ic++) delete [] eltStrings[ic];
  // RELEASESHAREDU(tpl, f2, cn2);

  // Cree le nom de zone
  for (size_t i = 0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }

  return 0;
}

//=============================================================================
// Only write BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA.
//=============================================================================
E_Int K_IO::GenIO::binvtuwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA supported
    if (eltTypes[zone][0] > 0 && eltTypes[zone][0] < 8) nvalidZones++;
    else
      printf("Warning: binvtuwrite: zone " SF_D_ " not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: binvtuwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  printf("Warning: binvtuwrite: not implemented.");
  return 1;
}
