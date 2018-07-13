#include "Def/Global/DefInclude.h"
#include "Def/Global/DefCompiler.h"
#include "Fld/Base/FldArray.h"
#include "Gen/IO/GenIO.h"

/* Add iblank array for current block in generic.ib file */
void writeIBlankArray(FldArrayI& nit, FldArrayI& njt, 
                      vector<FldArrayF*>& field);
