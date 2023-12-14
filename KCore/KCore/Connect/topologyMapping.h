/*    
    Copyright 2013-2024 Onera.

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
#include <vector>
#include <algorithm>

/* Topology structs */
// Topology object designed to work with std::map and std::unordered_map
// No explicit limit on the number of items contained in the Topology
struct Topology
{
  std::size_t n_ = 0; // actual number of items contained in the Topology
  std::vector<E_Int> p_;

  Topology() {}
  Topology(const E_Int* p) { set(p); }
  Topology(const std::vector<E_Int>& p) { set(p); }
  
  void set(const E_Int* p)
  {
    n_ = 0;
    for (const E_Int* tmp_p = p; *tmp_p != 0; ++tmp_p) n_++;
    p_.assign(p, p+n_);
    std::sort(p_.begin(), p_.end());
  }
  
  void set(const std::vector<E_Int>& p)
  {
    n_ = p.size(); p_ = p;
    std::sort(p_.begin(), p_.end());
  }

  E_Bool operator<(const Topology& other) const
  {
    std::size_t minsize = std::min(n_, other.n_);
    for (std::size_t i = 0; i < minsize; i++)
    {
      if (p_[i] < other.p_[i]) return true;
      else if (p_[i] > other.p_[i]) return false;
    }

    // If all common items are equal, the shorter vector is considered
    // smaller
    return n_ < other.n_;
  }

  E_Bool operator==(const Topology& other) const
  {
    if (n_ != other.n_) return false;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (p_[i] != other.p_[i]) return false;
    }
    return true;
  }
};

// Same as Topology but optimised (no dynamic allocation)
const std::size_t nmaxitems = 8;
struct TopologyOpt
{
  std::size_t n_ = 0; // actual number of items contained in the Topology
  E_Int p_[nmaxitems];

  TopologyOpt() {}
  TopologyOpt(const E_Int* p, const std::size_t n) { set(p, n); }
  TopologyOpt(const std::vector<E_Int>& p, const std::size_t n) { set(p, n); }
  
  void set(const E_Int* p, const std::size_t n)
  {
    n_ = n; assert(n_ <= nmaxitems);
    for (std::size_t i = 0; i < n_; i++) p_[i] = p[i];
    std::sort(p_, p_ + n_);
  }
  
  void set(const std::vector<E_Int>& p, const std::size_t n)
  {
    n_ = n; assert(n_ <= nmaxitems);
    for (std::size_t i = 0; i < n_; i++) p_[i] = p[i];
    std::sort(p_, p_ + n_);
  }

  E_Bool operator<(const TopologyOpt& other) const
  {
    std::size_t minsize = std::min(n_, other.n_);
    for (std::size_t i = 0; i < minsize; i++)
    {
      if (p_[i] < other.p_[i]) return true;
      else if (p_[i] > other.p_[i]) return false;
    }

    // If all common items are equal, the shorter array is considered smaller
    return n_ < other.n_;
  }

  E_Bool operator==(const TopologyOpt& other) const
  {
    if (n_ != other.n_) return false;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (p_[i] != other.p_[i]) return false;
    }
    return true;
  }
};