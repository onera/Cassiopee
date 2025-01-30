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
#ifndef _XCORE_XMPI_STATUS_HPP_
#define _XCORE_XMPI_STATUS_HPP_
#include <array>

#include "xmpi/constantes.hpp"

namespace xcore
{
#      if defined(_MPI)
        struct status
        {
        public:
            MPI_Status mpi_status;

            status() = default;
            status( const MPI_Status& st ) : mpi_status(st)
            {}
            status( const status& ) = default;

            status& operator = ( const status& ) = default;
            /**
             * @brief Retourne le nombre d'elements recus par une reception
             * @details Retourne via le status de la reception le nombre d'elements recus
             *          par la reception associee au status
             *          Usage : 
             *             status st;
             *             ...
             *             st = xcore::recv(double_array, expediteur);
             *             st.count<double>();// Nombre de reels double precision recus
             * @return Le nombre d'elements recus
             */
            template<typename K> int count() const
            {
                int cnt; MPI_Get_count(&mpi_status, Type_MPI<K>::mpi_type(), &cnt );
                return cnt;
            }
            /**
             * @brief Retourne le numero de l'expediteur du message
             */
            int source() const { return mpi_status.MPI_SOURCE; }
            /**
             * @brief Retourne l'identificateur associe au message
             */
            int tag()    const { return mpi_status.MPI_TAG;    }
            /**
             * @brief Retourne le numero d'erreur associe a la reception
             * @details Retourne un numero associe a une reussite ou une erreur durant
             *          la transmission du message recu. Si l'erreur renvoye est egal a
             *          xcore::success, cela veut dire que la reception du message s'est
             *          deroule correctement, sinon, les erreurs possibles sont :
             *               - xmpi::count dans le cas ou le nombre d'elements attendus etait inferieur au nombre d'elements recus
             *               - xmpi::rank  si le numero de l'expediteur est un numero invalide ( negatif ou superieur ou egal au nombre de processus )
             *               - xmpi::tag   si l'identificateur du message est un numero invalide ( ? )
             *               - xmpi::buffer si mpi n'a pas pu allouer un buffer temporaire pour recevoir le message
             * @return le numero d'erreur
             */
            int error()  const { return mpi_status.MPI_ERROR;  }
        };
#      else
        struct status
        {
            status() = default;
            status( int ct, int tg, int err ) : m_count(ct), m_tag(tg), m_error(err)
            {}
            status( const std::array<int,3>& st ) : m_count(st[0]), m_tag(st[1]), m_error(st[2])
            {}

            template<typename K> int count() const { return m_count; }
            int tag() const { return m_tag; }
            int source() const { return 0; }
            int error() const { return m_error; }
        private:
            int m_count, m_tag, m_error;
        };
#      endif
}

#endif
