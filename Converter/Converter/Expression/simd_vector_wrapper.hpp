#ifndef _CORE_HPC_SIMD_VECTOR_WRAPPER_HPP_
#define _CORE_HPC_SIMD_VECTOR_WRAPPER_HPP_
#include "Memory/vector_view.hpp"
#include <array>
#include <cassert>
#include <iostream>
#include <cmath>

namespace Expression {
    /**
     * @brief A wrapper on data of an array
     * @details simd_vector_wrapper instance wraps some existing data or
     * contains some local temporary data. Usage :
     * -------
     *     std::vector<double> u(256);
     *     simd_vector_wrapper u(3.);
     *     simd_vector_wrapper v({1.,2.,3.,4.,6.});
     *     simd_vector_wrapper w(u.data());
     *
     */
    class simd_vector_wrapper {
      public:
        static constexpr std::size_t max_size = 512;
        using iterator       = K_MEMORY::vector_view<double>::iterator;
        using const_iterator = K_MEMORY::vector_view<double>::const_iterator;

        /** @name Constructors and destructor
         *
         * Constructors to initialize a simd_vector_wrapper instance.
         */
        ///@{
        /**
         * @brief Default constructor
         * @details Build a new simd vector which has owned not yet initialized
         * data.
         *
         */
        simd_vector_wrapper()
            : m_own_data(), m_pt_data(m_own_data.data(), max_size) {
            assert(has_owned_data() &&
                   "Error, simd vector must own here his data !");
        }
        /**
         * @brief Build a simd vector wrapping data given by a pointer
         * @details Build a new simd vector which not owns data but wraps
         *          preexisting data.
         *
         * @param pt_data The pointer on the first element of the data to wrap.
         */
        simd_vector_wrapper(double *pt_data)
            : m_own_data(), m_pt_data(pt_data, max_size) {}
        /**
         * @brief build a simd vector with owned data initialized from arr
         * @details build a new simd vector with data in arr array. The data are
         *          owned by the simd vector.
         *
         * @param arr data to copy in simd vector
         */
        simd_vector_wrapper(const std::array<double, max_size> &arr)
            : m_own_data(arr), m_pt_data(m_own_data.data(), max_size) {
            assert(has_owned_data() &&
                   "Error, simd vector must own here his data !");
        }

        /**
         * @brief Broadcast a scalar in simd vector
         * @details Convert a simple double as simd vector with each element
         *          equal to the scalar.
         *
         * @param scal The scalar to broadcast
         */
        simd_vector_wrapper(double scal)
            : m_own_data(), m_pt_data(m_own_data.data(), max_size) {
            assert(has_owned_data() &&
                   "Error, simd vector must own here his data !");
            double *pt = m_own_data.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt[i] = scal;
        }

        /**
         * @brief Copy constructor
         * @details Copy data from u in new simd vector. If the data is owned by
         * u, the constructor effectively copy the data from u in the new
         * vector, else the new simd vector wraps and shares the same
         * data than u
         *
         * @param u The simd vector to copy
         */
        simd_vector_wrapper(const simd_vector_wrapper &u)
            : m_own_data(), m_pt_data() {
            if (u.has_owned_data()) {
                m_own_data = u.m_own_data;
                m_pt_data  = {m_own_data.data(), max_size};
            } else {
                m_pt_data = u.m_pt_data;
            }
        }

        ///! Default move constructor
        simd_vector_wrapper(simd_vector_wrapper &&u)
            : m_own_data(), m_pt_data() {
            if (u.has_owned_data()) {
                m_own_data = std::move(u.m_own_data);
                m_pt_data  = {m_own_data.data(), max_size};
                assert(has_owned_data() &&
                       "Error, simd vector must own here his data !");
            } else {
                m_pt_data = std::move(u.m_pt_data);
            }
        }

        ///! Default destructor
        ~simd_vector_wrapper() = default;
        ///@}

        /** @name Getters and setters
         */
        ///@{

        /**
         * @brief Tell if the simd vector owns his data
         * @details Return true if the simd vector owns his data, false else.
         * @return True if the data is owned by the simd vector
         */
        bool has_owned_data() const {
            return m_own_data.data() == m_pt_data.data();
        }

        /**
         * @brief Return the size of the simd vector
         * @details Return the size of the simd vector. All simd vector have
         * same size, a.k.a max_size constant
         * @return The size of the simd vector
         */
        std::size_t size() const { return max_size; }
        /**
         * @brief Return a pod pointer on the first element of the simd vector
         * @details Return a pod pointer on the first element of the simd
         * vector, owned or not.
         * @return The address of the first element.
         */
        ///@{
        double *      data() { return m_pt_data.data(); }
        const double *data() const { return m_pt_data.data(); }
        ///@}

        /**
         * @brief Return an iterator on the beginning of the simd vector
         * @details Return a writable iterator on the beginning of the simd
         * vector
         * @return An iterator
         */
        iterator begin() { return m_pt_data.begin(); }

        /**
         * @brief Return a const_iterator on the beginning of the simd vector
         * @details Return a only readable iterator on the beginning of the simd
         * vector
         * @return A const_iterator
         */
        ///@{
        const_iterator cbegin() const { return m_pt_data.cbegin(); }
        const_iterator begin() const { return m_pt_data.cbegin(); }
        ///@}

        /**
         * @brief Return an iterator on the end of the simd vector
         * @details Return a writable iterator on the end of the simd vector
         * @return An iterator
         */
        iterator end() { return m_pt_data.end(); }
        /**
         * @brief Return a const_iterator on the end of the simd vector
         * @details Return a only readable iterator on the end of the simd
         * vector
         * @return A const_iterator
         */
        ///@{
        const_iterator cend() const { return m_pt_data.cend(); }
        const_iterator end() const { return m_pt_data.cend(); }

        /**
         * @brief Copy data in current memory from another simd_vector_wrapper
         * @details A faire
         * 
         * @param v the input vector
         */
        void copy_in( const simd_vector_wrapper& v )
        {
            const double * __restrict__ dv = v.data();
            double * __restrict__     du = data();
            #pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i) {
                du[i] = dv[i];
            }
        }     
        ///@}

        /** @name Operators
         */
        ///@{
        simd_vector_wrapper &operator=(const simd_vector_wrapper &v) {
            if (this->data() != v.data()) {
                if (has_owned_data()) {
                    const double *dv = v.data();
                    double *      du = data();
#pragma omp simd
                    for (std::size_t i = 0; i < max_size; ++i) {
                        du[i] = dv[i];
                    }
                    m_pt_data = {m_own_data.data(), max_size};
                } else if (v.has_owned_data()) {
                    m_own_data = v.m_own_data;
                    m_pt_data  = {m_own_data.data(), max_size};
                } else {
                    m_pt_data = v.m_pt_data;
                }
            }
            return *this;
        }

        simd_vector_wrapper &operator=(simd_vector_wrapper &&v) {
            if (this->data() != v.data()) {
                if (v.has_owned_data()) {
                    m_own_data = std::move(v.m_own_data);
                    m_pt_data  = {m_own_data.data(), max_size};
                } else if (has_owned_data()) {
                    const double *dv = v.data();
                    double *      du = data();
#pragma omp simd
                    for (std::size_t i = 0; i < max_size; ++i) {
                        du[i] = dv[i];
                    }
                } else
                    m_pt_data = std::move(v.m_pt_data);
            }
            return *this;
        }

        simd_vector_wrapper &operator+=(const simd_vector_wrapper &v)  {
            double *      pt_u = data();
            const double *pt_v = v.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i) {
                pt_u[i] += pt_v[i];
            }
            return *this;
        }

        simd_vector_wrapper operator+(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *__restrict__ pt_u = data();
            const double *__restrict__ pt_v = v.data();
            double *__restrict__ pt_w       = w.m_own_data.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = pt_u[i] + pt_v[i];
            return w;
        }

        simd_vector_wrapper &operator-=(const simd_vector_wrapper &v)  {
            double *__restrict__ pt_u       = data();
            const double *__restrict__ pt_v = v.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i) {
                pt_u[i] -= pt_v[i];
            }
            return *this;
        }

        simd_vector_wrapper operator-(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double * __restrict__ pt_u = data();
            const double * __restrict__ pt_v = v.data();
            double *__restrict__ pt_w       = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = pt_u[i] - pt_v[i];
            return w;
        }

        simd_vector_wrapper operator-() const {
            simd_vector_wrapper w;
            const double * __restrict__     pt_u = data();
            double * __restrict__           pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = -pt_u[i];
            return w;
        }

        simd_vector_wrapper &operator*=(const simd_vector_wrapper &v)  {
            double * __restrict__     pt_u = data();
            const double * __restrict__ pt_v = v.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i) {
                pt_u[i] *= pt_v[i];
            }
            return *this;
        }

        simd_vector_wrapper operator*(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double * __restrict__     pt_u = data();
            const double * __restrict__ pt_v = v.data();
            double * __restrict__            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = pt_u[i] * pt_v[i];
            return w;
        }

        simd_vector_wrapper &operator/=(const simd_vector_wrapper &v)  {
            double * __restrict__      pt_u = data();
            const double * __restrict__ pt_v = v.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i) {
                pt_u[i] /= pt_v[i];
            }
            return *this;
        }

        simd_vector_wrapper operator/(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double * __restrict__      pt_u = data();
            const double * __restrict__ pt_v = v.data();
            double * __restrict__            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = pt_u[i] / pt_v[i];
            return w;
        }

        simd_vector_wrapper operator==(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double * __restrict__      pt_u = data();
            const double * __restrict__ pt_v = v.data();
            double * __restrict__             pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] == pt_v[i]);
            return w;
        }

        simd_vector_wrapper operator!=(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *      pt_u = data(), *pt_v = v.data();
            double *            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] != pt_v[i]);
            return w;
        }

        simd_vector_wrapper operator>(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *      pt_u = data(), *pt_v = v.data();
            double *            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] > pt_v[i]);
            return w;
        }

        simd_vector_wrapper operator>=(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *      pt_u = data(), *pt_v = v.data();
            double *            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] >= pt_v[i]);
            return w;
        }

        simd_vector_wrapper operator<(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *      pt_u = data(), *pt_v = v.data();
            double *            pt_w = w.data();
#pragma omp simd
            for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] < pt_v[i]);
            return w;
        }

        simd_vector_wrapper operator<=(const simd_vector_wrapper &v) const
             {
            simd_vector_wrapper w;
            const double *__restrict__ pt_u = data();
            const double *__restrict__ pt_v = v.data();
            double *__restrict__ pt_w       = w.data();
#pragma omp simd
              for (std::size_t i = 0; i < max_size; ++i)
                pt_w[i] = double(pt_u[i] <= pt_v[i]);
            return w;
        }

        ///@}

      private:
        alignas(64) std::array<double, max_size>  m_own_data;
        K_MEMORY::vector_view<double> m_pt_data;
    };
} // namespace Expression

#endif
