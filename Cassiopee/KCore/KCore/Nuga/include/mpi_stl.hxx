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

//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_MPI_STL_HXX
#define NUGA_MPI_STL_HXX

#include <map>
#include <vector>
#include "kcore.h"
#include "Nuga/include/macros.h"

#ifdef _MPI
#include "mpi.h"

namespace NUGA
{
  namespace MPI
  {

    template <typename T>
    static inline MPI_Datatype mpi_get_type() noexcept 
    {
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

#else
#define MPI_INT int
#define MPI_Datatype int
#define MPI_LONG long
#define MPI_LONG_LONG long
#define MPI_FLOAT float
#define MPI_DOUBLE double
#define MPI_LONG_DOUBLE double
#define MPI_C_BOOL int
#define MPI_CHAR char
#define MPI_SIGNED_CHAR char
#define MPI_UNSIGNED_CHAR char
#define MPI_WCHAR char
#define MPI_SHORT short
#define MPI_UNSIGNED_SHORT short
#define MPI_UNSIGNED int
#define MPI_UNSIGNED_LONG long
#define MPI_UNSIGNED_LONG_LONG long
#define MPI_INT8_T int
#define MPI_INT16_T int
#define MPI_INT32_T int
#define MPI_INT64_T int
#define MPI_UINT8_T int
#define MPI_UINT16_T int
#define MPI_UINT32_T int
#define MPI_UINT64_T int
#define MPI_Barrier(X)
#define MPI_Comm_rank(X1,X2)
#define MPI_Comm_size(X1,X2)
#define MPI_Recv(X1,X2,X3,X4,X5,X6,X7)
#define MPI_Waitall(X1,X2,X3)
#define MPI_Request int
#define MPI_Isend(X1,X2,X3,X4,X5,X6,X7)
#define MPI_Send(X1,X2,X3,X4,X5,X6)
#define MPI_Irecv(X1,X2,X3,X4,X5,X6,X7) 0
#define MPI_Comm int
#define MPI_Status int
#define MPI_DATATYPE_NULL int
#endif

#endif