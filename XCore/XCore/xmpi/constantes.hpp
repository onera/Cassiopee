#ifndef _XCORE_XMPI_CONSTANTES_HPP_
#define _XCORE_XMPI_CONSTANTES_HPP_

#include <complex>

# include "xmpi/x_mpi.h"

namespace xcore
{
#if defined( _MPI )
        const int any_tag = MPI_ANY_TAG;
        const int any_source = MPI_ANY_SOURCE;
        const int undefined = MPI_UNDEFINED;

        enum error {
            success = MPI_SUCCESS,
            count   = MPI_ERR_COUNT,
            rank    = MPI_ERR_RANK,
            tag     = MPI_ERR_TAG,
            buffer  = MPI_ERR_BUFFER,
            unknown = MPI_ERR_UNKNOWN
        };

        typedef MPI_Op Operation;

        const Operation null = MPI_OP_NULL;
        const Operation max = MPI_MAX;
        const Operation min = MPI_MIN;
        const Operation sum = MPI_SUM;
        const Operation prod = MPI_PROD;
        const Operation logical_and = MPI_LAND;
        const Operation binary_and = MPI_BAND;
        const Operation logical_or = MPI_LOR;
        const Operation binary_or = MPI_BOR;
        const Operation logical_xor = MPI_LXOR;
        const Operation binary_xor = MPI_BXOR;
        const Operation minloc = MPI_MINLOC;
        const Operation maxloc = MPI_MAXLOC;
        const Operation replace = MPI_REPLACE;

        typedef MPI_Comm Ext_Communicator;

        template <typename K>
        struct Type_MPI {
            static bool must_be_packed() { return true; }
            static MPI_Datatype mpi_type() { return MPI_PACKED; }
        };
        //
        template <>
        struct Type_MPI<short> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_SHORT; }
        };
        //
        template <>
        struct Type_MPI<char> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_CHAR; }
        };
        //
        template <>
        struct Type_MPI<int> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_INT; }
        };
        //
        template <>
        struct Type_MPI<long> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_LONG; }
        };
        //
        template <>
        struct Type_MPI<float> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_FLOAT; }
        };
        //
        template <>
        struct Type_MPI<double> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_DOUBLE; }
        };
        //
        template <>
        struct Type_MPI<unsigned char> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_BYTE; }
        };
        //
        template <>
        struct Type_MPI<unsigned short> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_UNSIGNED_SHORT; }
        };
        //
        template <>
        struct Type_MPI<unsigned> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_UNSIGNED; }
        };
        //
        template <>
        struct Type_MPI<unsigned long> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_UNSIGNED_LONG; }
        };
        template <>
        struct Type_MPI<std::complex<float>> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_COMPLEX; }
        };

        template <>
        struct Type_MPI<std::complex<double>> {
            static bool must_be_packed() { return false; }
            static MPI_Datatype mpi_type() { return MPI_DOUBLE_COMPLEX; }
        };
#else
        const int any_tag = -1;    /*!< Constant to receive from any tag */
        const int any_source = -1; /*!< Constant to receive from any source */
        const int undefined = -1;  /*!< Constant for undefined parameter */

        typedef int Ext_Communicator;

        typedef int Operation;

        enum operation {
            null = 0,
            max,
            min,
            sum,
            prod,
            logical_and,
            binary_and,
            logical_or,
            binary_or,
            logical_xor,
            binary_xor,
            minloc,
            maxloc,
            replace
        };
        /*!
     * \enum error
     * \brief To manage error coming for parallel calls
     */
        enum error {
            success = 0, /*!< Success. No error */
            count,       /*!< Error on the size of receiving data */
            rank,        /*!< Wrong rank. Doesn't exist */
            tag,         /*!< Wrong tag. */
            buffer,      /*!< Problems with the used buffer (null ? )*/
            unknown      /*!< Unknown problem */
        };
#endif
}  // namespace xcore
#endif
