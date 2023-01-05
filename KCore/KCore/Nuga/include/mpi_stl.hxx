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
    static inline MPI_Datatype mpi_get_type() noexcept {

      if (std::is_same<T, signed int>::value) return MPI_INT;
      if (std::is_same<T, signed long int>::value) return MPI_LONG;
      if (std::is_same<T, signed long long int>::value) return MPI_LONG_LONG;

      if (std::is_same<T, float>::value) return MPI_FLOAT;
      if (std::is_same<T, double>::value) return MPI_DOUBLE;
      if (std::is_same<T, long double>::value) return MPI_LONG_DOUBLE;

      if (std::is_same<T, bool>::value) return MPI_C_BOOL;
      
      if (std::is_same<T, char>::value) return MPI_CHAR;
      if (std::is_same<T, signed char>::value) return MPI_SIGNED_CHAR;
      if (std::is_same<T, unsigned char>::value) return MPI_UNSIGNED_CHAR;
      if (std::is_same<T, wchar_t>::value) return MPI_WCHAR;

      if (std::is_same<T, signed short>::value) return MPI_SHORT;
      if (std::is_same<T, unsigned short>::value) return MPI_UNSIGNED_SHORT;
      if (std::is_same<T, unsigned int>::value) return MPI_UNSIGNED;
      if (std::is_same<T, unsigned long int>::value) return MPI_UNSIGNED_LONG;
      if (std::is_same<T, unsigned long long int>::value) return MPI_UNSIGNED_LONG_LONG;

      if (std::is_same<T, std::int8_t>::value) return MPI_INT8_T;
      if (std::is_same<T, std::int16_t>::value) return MPI_INT16_T;
      if (std::is_same<T, std::int32_t>::value) return MPI_INT32_T;
      if (std::is_same<T, std::int64_t>::value) return MPI_INT64_T;
      if (std::is_same<T, std::uint8_t>::value) return MPI_UINT8_T;
      if (std::is_same<T, std::uint16_t>::value) return MPI_UINT16_T;
      if (std::is_same<T, std::uint32_t>::value) return MPI_UINT32_T;
      if (std::is_same<T, std::uint64_t>::value) return MPI_UINT64_T;

      /*
      if (std::is_same<T, std::complex<float>>::value) return MPI_C_COMPLEX;
      if (std::is_same<T, std::complex<double>>::value) return MPI_C_DOUBLE_COMPLEX;
      if (std::is_same<T, std::complex<long double>>::value) return MPI_C_LONG_DOUBLE_COMPLEX;
      */
      return MPI_DATATYPE_NULL;
    }


    template <typename T>
    int Isend(const std::vector<T>& data, int target_rank_id, int TAG_1_BASED, MPI_Comm COM, MPI_Request* req)
    {
      int data_sz = data.size();
      int err = MPI_Send(&data_sz, 1, mpi_get_type<T>(), target_rank_id, TAG_1_BASED, COM);
      if (err) return 1;

      err = MPI_Isend(&data[0], data_sz, mpi_get_type<T>(), target_rank_id, TAG_1_BASED+1, COM, req);
      return err;
    }

    template <typename T>
    int Irecv(int sz, int source_rank_id, int TAG_1_BASED, MPI_Comm COM, std::vector<T>& data, MPI_Request* req)
    {
      data.resize(sz);
      return MPI_Irecv(&data[0], sz, mpi_get_type<T>(), source_rank_id, TAG_1_BASED, COM, req);
    }

  }
}



#endif
