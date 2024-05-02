#include "Proto.h"

void compute_principal_vecs(AMesh *M, std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    switch (M->cellTree->type(i)) {
      case TETRA:
        make_pdirs_tetra(i, M, Dirs[i]);
        break;
      case PENTA:
        make_pdirs_penta(i, M, Dirs[i]);
        break;
      case PYRA:
        make_pdirs_pyra(i, M, Dirs[i]);
        break;
      case HEXA:
        make_pdirs_hexa(i, M, Dirs[i]);
        break;
      default:
        assert(0);
        break;
    }
  }
}