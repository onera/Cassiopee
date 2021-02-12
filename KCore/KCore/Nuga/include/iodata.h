/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __DynArrayIO_MDATA_H__
#define __DynArrayIO_MDATA_H__

#include "MeshData.h"

namespace DELAUNAY
{
class iodata
{
public:
  static E_Int read (const char* filename, DELAUNAY::MeshData& data);
  static E_Int write (const char* filename, const DELAUNAY::MeshData& data);
private:
  iodata(void);
  ~iodata(void);
};
}

#endif
