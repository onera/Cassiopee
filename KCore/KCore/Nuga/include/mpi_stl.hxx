/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_MPI_STL_HXX
#define NUGA_MPI_STL_HXX

#include "mpi.h"

namespace NUGA
{

  namespace MPI
  {

    template <typename T>
    struct mpi_type_trait;

    template <>
    struct mpi_type_trait<E_Int>
    {
      static constexpr MPI_Datatype type = MPI_INT;
    };


    template <typename T>
    int Isend(const std::vector<T>& data, int target_rank_id, int TAG_1_BASED, MPI_Comm COM, MPI_Request* req)
    {
      int data_sz = data.size();
      int err = MPI_Send(&data_sz, 1, mpi_type_trait<T>::type, target_rank_id, TAG_1_BASED, COM);
      if (err) return 1;

      err = MPI_Isend(&data[0], data_sz, mpi_type_trait<T>::type, target_rank_id, TAG_1_BASED+1, COM, req);
      return err;
    }

    template <typename T>
    int Irecv(int sz, int source_rank_id, int TAG_1_BASED, MPI_Comm COM, std::vector<T>& data, MPI_Request* req)
    {
      data.resize(sz);
      return MPI_Irecv(&data[0], sz, mpi_type_trait<T>::type, source_rank_id, TAG_1_BASED, COM, req);
    }

  }
}



#endif
