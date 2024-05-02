# include <cstdlib>
# include "DataDL.h"
//# include "DataVBO.h"
// ########################################################################
Data* Data::getInstance()
{
  switch(Data::_renderID) {
  case Data::DL:
    return DataDL::getInstance();
  case Data::VBO:
  default:
    exit(EXIT_FAILURE);
    return NULL;
  }
}
