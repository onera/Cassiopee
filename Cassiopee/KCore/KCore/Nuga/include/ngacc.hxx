/*    
    Copyright 2013-2026 ONERA.

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

//Authors : Sam Landier (sam.landier@onera.fr)

#include "Nuga/include/ngon_t.hxx"

namespace K_FLD
{
template <typename cnt_t>
  class ArrayAccessor<ngon_t<cnt_t> >
  {

  public: /** Typedefs */
    using ngon_type = ngon_t<cnt_t>;
    typedef           E_Int   size_type;

  public:
    /// Constuctor
    explicit ArrayAccessor(const ngon_type& ng, E_Int shift = 0):_ng(ng)
    {}

    /// Destructor
    ~ArrayAccessor(){};

    /// Returns the total number of entries (points/elements).
    inline E_Int size() const {return _ng.PHs.size();}
    
    /// Returns the number of used rows for the i-th element (NGON specific)
    inline E_Int stride(E_Int i) const {return _ng.PHs.stride(i);}
            

    /// Returns the j-th entry (j-th column).
    template <typename ELT>
    inline void getEntry(const E_Int& j, ELT& PHj) const
    {
      PHj._pgs=const_cast<ngon_unit*>(&_ng.PGs);
      PHj._faces=const_cast<E_Int*>(_ng.PHs.get_facets_ptr(j));
      PHj._nb_faces=_ng.PHs.stride(j);
    }
    
    /// Returns the j-th entry's pointer to the first field.
    //inline const E_Int* getEntry(const E_Int& j) const { }
    
    private:
      const ngon_type& _ng;
  };
}