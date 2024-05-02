#ifndef _KCORE_NUMPY_VECTOR_HPP_
#define _KCORE_NUMPY_VECTOR_HPP_
#include "Numpy/Numpy.h"
#include "Numpy/types.hpp"
#include <algorithm>
#include <vector>

namespace KCore {
    namespace Numpy {
        template <typename K, class Allocator = std::allocator<K>>
        class vector {
          public:
            typedef K                                     value_type;
            typedef Allocator                             allocator_type;
            typedef K *                                   pointer;
            typedef const K *                             const_pointer;
            typedef K &                                   reference;
            typedef const K &                             const_reference;
            typedef pointer                               iterator;
            typedef const_pointer                         const_iterator;
            typedef std::reverse_iterator<iterator>       reverse_iterator;
            typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
            typedef std::size_t                           size_type;
            typedef std::ptrdiff_t                        difference_type;
            typedef PyArrayObject *                       numpy_pointer;
            typedef const PyArrayObject *                 const_numpy_pointer;

            // =================================================================
            // Constructeurs et destructeur
            /**
             * @brief Construit un vecteur vide
             * @details Construit un vecteur non alloué et non associé à un tableau numpy
             *
             */
            vector() : m_size(0), m_pt_np_array(nullptr) {}

            /**
             * @brief Construit un vecteur monodimensionnel
             * @details Construit un vecteur monodimensionnel de taille sz sans initialisation des valeurs
             *
             * @param sz Taille du tableau
             */
            vector(const size_type &sz) : m_size(sz) {
                npy_intp dim[1];
                dim[0]        = sz;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(1, dim, KCore::Numpy::type<K>::value(), 0);
            }

            /**
             * @brief Construit un tableau monodimensionnel de taille sz
             * @details Construit un tableau monodimensionnel de taille sz remplit avec une valeur val
             *
             * @param sz  La taille du tableau
             * @param val La valeur de remplissage du tableau
             */
            explicit vector(size_type sz, const K &val) : m_size(sz) {
                npy_intp dim[1];
                dim[0]        = sz;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(1, dim, KCore::Numpy::type<K>::value(), 0);
                std::fill_n((K *)PyArray_DATA(m_pt_np_array), sz, val);
            }

            /**
             * @brief Construit un tableau multi dimensionnel donné par la liste d'initialisation shape
             * @details Construit un tableau multi dimensionnel donné par la liste shape sous la forme :
             *              vector<double> np_vect({3,5,7}); // Construit un tableau de dimension 3 x 5 x 7
             *
             * @param shape Dimensions du vecteur
             */
            explicit vector(std::initializer_list<size_type> shape) : m_size(1) {
                std::vector<npy_intp> dim(shape.begin(), shape.end());
                for (auto d : dim)
                    m_size *= d;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(shape.size(), dim.data(),
                                                               KCore::Numpy::type<K>::value(), 0);
            }

            /**
             * @brief Construit un tableau multi dimensionnel initialisé donné par la liste d'initialisation shape
             * @details Construit un tableau multi-dimensionnel initialisé par le paramètre val et dont la forme
             * multidimensionnel est donné par le paramèter shape
             *          Usage : vector<double> np_vect({3,5,7}, 3.14);
             *
             * @param dim [description]
             * @param e [description]
             */
            explicit vector(std::initializer_list<size_type> shape, const K &val) : m_size(1) {
                std::vector<npy_intp> dim(shape.begin(), shape.end());
                for (auto d : dim)
                    m_size *= d;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(shape.size(), dim.data(),
                                                               KCore::Numpy::type<K>::value(), 0);

                std::fill_n((K *)PyArray_DATA(m_pt_np_array), m_size, val);
            }

            /**
             * @brief Construit un vecteur initialisé par les données comprises entre beg (inclu) et end (exclu)
             * @details Construit un vecteur initialisé par des données comprises entre deux itérateurs ( qui
             * peuvent être des pointeurs ), beg (inclu) pour le début des données à copier et end(exclu) la fin
             * des données à copier.
             *
             * @param beg Itérateur sur le premier élément à copier
             * @param end Itérateur sur l'élément suivant le dernier élément à copier
             */
            template <typename InputRandomIterator>
            explicit vector(const InputRandomIterator &beg, const InputRandomIterator &end) : m_size(end - beg) {
                npy_intp dim[1];
                dim[0]        = m_size;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(1, dim, KCore::Numpy::type<K>::value(), 0);
                std::copy(beg, end, (K *)PyArray_DATA(m_pt_np_array));
            }

            /**
             * @brief Converti un tableau de type C++ (std::vector) en tableau numpy
             * @details Converti un tableau std::vector en un tableau numpy mono dimensionnel
             *
             * @param arr Le tableau des éléments à copier
             */
            explicit vector(const std::vector<K> &arr) : m_size(arr.size()) {
                npy_intp dim[1];
                dim[0]        = m_size;
                m_pt_np_array = (PyArrayObject *)PyArray_EMPTY(1, dim, KCore::Numpy::type<K>::value(), 0);
                std::copy(arr.begin(), arr.end(), (K *)PyArray_DATA(m_pt_np_array));
            }

            /**
             * @brief Encapsule un tableau numpy dans un vecteur de type numpy
             * @details Encapsule sans copie un tableau numpy dans un vecteur de type numpy.
             * Le compteur de référence du tableau numpy est incrémenté de un pour assurer l'existence
             * des données du tableau de type vector si le tableau numpy du côté python aurait dû être détruit
             * par Python.
             *
             * @param pt_numpy_array Le pointeur sur le tableau numpy
             */
            vector(numpy_pointer pt_numpy_array) : m_size(0), m_pt_np_array(pt_numpy_array) {
                if (pt_numpy_array != nullptr) {
                    Py_XINCREF(pt_numpy_array);
                    m_size = PyArray_Size((PyObject *)pt_numpy_array);
                }
            }

            /*
            vector( PyObject* obj ) : m_size( 0 ), m_pt_np_array( nullptr )
            {
                if ( PyArray_Check( obj ) ) {
                    m_pt_np_array = (PyArrayObject*)obj;
                    Py_XINCREF( m_pt_np_array );
                    m_size = PyArray_Size( (PyObject*)m_pt_np_array );
                } else {
                    assert( PySequence_Check( obj ) );
                    m_size = Bemuse::Python::sequence_size( obj );
                    npy_intp dim[ 1 ];
                    dim[ 0 ] = npy_intp( m_size );
                    m_pt_np_array =
                        (PyArrayObject*)PyArray_EMPTY( 1, dim, Bemuse::Python::numpy::type<K>::value(), 0 );
                    Bemuse::Python::sequence2array( obj, m_pt_np_array->data );
                }
            }*/

            /**
             * @brief Constructeur de copie
             * @details Copie un tableau numpy encapsulé dans un vecteur dans un nouveau vecteur
             * contenant un nouveau tableau numpy
             *
             * @param other Le vecteur à copier
             */
            explicit vector(const vector &other) : m_size(other.m_size) {
                m_pt_np_array = (PyArrayObject *)PyArray_NewCopy(other.m_pt_np_array, NPY_ANYORDER);
            }

            /**
             * @brief Constructeur de déplacement
             */
            explicit vector(vector &&other) : m_size(other.m_size), m_pt_np_array(other.m_pt_np_array) {
                other.m_size        = 0;
                other.m_pt_np_array = nullptr;
            }

            /**
             * @brief Destructeur
             * @details Décrémente le compteur de un du tableau numpy encapsulé.
             */
            ~vector() { Py_XDECREF(m_pt_np_array); }

            // Opérations sur les vector :
            // ---------------------------
            vector &operator=(const vector &other) {
                if (this != &other) {
                    Py_XDECREF(m_pt_np_array);
                    m_size        = other.m_size;
                    m_pt_np_array = (PyArrayObject *)PyArray_NewCopy(other.m_pt_np_array, NPY_ANYORDER);
                }
                return *this;
            }
            // -----------------------------------------------------
            vector &operator=(vector &&other) {
                if (this != &other) {
                    Py_XDECREF(m_pt_np_array);
                    m_pt_np_array       = other.m_pt_np_array;
                    m_size              = other.m_size;
                    other.m_pt_np_array = nullptr;
                }
                return *this;
            }
            // -----------------------------------------------------
            iterator               begin() { return data(); }
            iterator               end() { return data() + m_size; }
            const_iterator         begin() const { return cbegin(); }
            const_iterator         end() const { return cend(); }
            const_iterator         cbegin() const { return data(); }
            const_iterator         cend() const { return data() + m_size; }
            reverse_iterator       rbegin() { return end() - 1; }
            reverse_iterator       rend() { return begin() - 1; }
            const_reverse_iterator crbegin() { return cend() - 1; }
            const_reverse_iterator crend() { return cbegin() - 1; }
            const_reverse_iterator rbegin() const { return crbegin(); }
            const_reverse_iterator rend() const { return crend(); }

            size_type size() const { return m_size; }
            size_type length() const { return m_size; }
            size_type max_size() const { return m_size; }
            bool      empty() const { return m_size == 0; }

            reference operator[](const size_type pos) {
                assert(pos < size());
                return data()[pos];
            }
            const_reference operator[](const size_type pos) const {
                assert(pos < size());
                return data()[pos];
            }
            reference at(size_type pos) {
                if (pos >= size()) throw std::out_of_range("Wrong index");
                return data()[pos];
            }
            const_reference at(size_type pos) const {
                if (pos >= size()) throw std::out_of_range("Wrong index");
                return data()[pos];
            }
            reference       front() { return data()[0]; }
            reference       back() { return data()[size() - 1]; }
            const_reference front() const { return data()[0]; }
            const_reference back() const { return data()[size() - 1]; }
            pointer         data() { return (pointer)PyArray_DATA(m_pt_np_array); }
            const_pointer   data() const { return (const_pointer)PyArray_DATA(m_pt_np_array); }

            void clear() { *this = vector(); }
            void swap(vector &s) {
                std::swap(m_pt_np_array, s.m_pt_np_array);
                std::swap(m_size, s.m_size);
            }

            int ref_counter() const { return Py_REFCNT(m_pt_np_array); }

            explicit operator std::vector<K>() const { return std::vector<K>(begin(), end()); }
                     operator numpy_pointer() { return m_pt_np_array; }
                     operator const_numpy_pointer() const { return m_pt_np_array; }

          private:
            size_type     m_size;
            numpy_pointer m_pt_np_array;
        };

        /**
         * @brief Return a numpy array to python.
         * @details Return a numpy array to python. This array is shared
         *          between C++ and Python, so we must increment the reference
         *          counter to be sure to destroy the array when neither C++,
         *          neither Python has a reference on the numpy array.
         *
         * @param arr The vector containing the numpy array
         * @tparam K The kind of element in the numpy array
         * @tparam Allocator The allocator. Only to be compatible with std::vector.
         */
        template <typename K, class Allocator>
        PyArrayObject *shared_with_python(vector<K, Allocator> &arr) {
            PyArrayObject *py_array = static_cast<PyArrayObject *>(arr);
            Py_XINCREF(py_array);
            return py_array;
        }
    } // namespace Numpy
} // namespace KCore

#endif
